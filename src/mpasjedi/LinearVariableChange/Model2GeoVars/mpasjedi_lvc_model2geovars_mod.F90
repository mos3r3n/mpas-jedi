! (C) Copyright 2020 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

module mpasjedi_lvc_model2geovars_mod

use iso_c_binding

use fckit_configuration_module, only: fckit_configuration
use fckit_log_module, only: fckit_log

!oops
use kinds, only : kind_real

!ufo
use gnssro_mod_transform, only: geometric2geop
use ufo_vars_mod

!MPAS-Model
use mpas_derived_types
use mpas_field_routines
use mpas_kind_types, only: RKIND
use mpas_pool_routines

!mpas-jedi
use mpas_constants_mod
use mpas_fields_mod, only: mpas_fields
use mpas_geom_mod, only: mpas_geom
use mpas4da_mod, only: da_template_pool
use mpas_hydro_pt_mod, only: nHydroPT, hydroBaseNames, is_hydro_base, hydro_pt_name, pt_slope

!MPAS-Model physical constants
use mpas_constants, only: gravity

!TODO: package the variable changes somewhere
use mpas2ufo_vars_mod
!use mpas_vc_utils

implicit none

private
public :: mpasjedi_lvc_model2geovars

type :: mpasjedi_lvc_model2geovars
  type(mpas_pool_type), pointer :: trajectory => null()
  contains
    procedure, public :: create
    procedure, public :: delete
    procedure, public :: multiply
    procedure, public :: multiplyadjoint
end type mpasjedi_lvc_model2geovars

! --------------------------------------------------------------------------------------------------

character(len=1024) :: message

contains

! --------------------------------------------------------------------------------------------------

subroutine create(self, geom, bg, fg, conf)

  class(mpasjedi_lvc_model2geovars), intent(inout) :: self
  type(mpas_geom),                   intent(in)    :: geom
  type(mpas_fields),                 intent(in)    :: bg
  type(mpas_fields),                 intent(in)    :: fg
  type(fckit_configuration),         intent(in)    :: conf

  integer :: iVar

  character(len=MAXVARLEN), parameter :: &
    trajFieldNames(4 + nHydroPT) = &
      [ character(len=MAXVARLEN) :: &
        'air_temperature', &
        'water_vapor_mixing_ratio_wrt_moist_air', &
        'air_pressure', &
        'air_pressure_at_surface', &
        hydroBaseNames & ! hydrometeors: linearization state for <name>_pt increments
      ]
  character(len=MAXVARLEN), allocatable :: trajHasNames(:)

  allocate(trajHasNames(0))

  !Check if the given trajectory variable exists in the background variables
  do iVar = 1, size(trajFieldNames)
    if (bg%has(trajFieldNames(iVar))) then
      trajHasNames = [ character(len=MAXVARLEN) :: trajHasNames, trajFieldNames(iVar) ]
    endif
  end do

  call da_template_pool(geom, self%trajectory, size(trajHasNames), trajHasNames)

  do iVar = 1, size(trajHasNames)
    call bg%copy_to(trajHasNames(iVar), self%trajectory)
  end do

  deallocate(trajHasNames)

end subroutine create

! --------------------------------------------------------------------------------------------------

subroutine delete(self)
  class(mpasjedi_lvc_model2geovars), intent(inout) :: self
  if (associated(self%trajectory)) then
    call mpas_pool_destroy_pool(self%trajectory)
  end if
end subroutine delete

! --------------------------------------------------------------------------------------------------

subroutine multiply(self, geom, dxm, dxg)
  class(mpasjedi_lvc_model2geovars), intent(inout) :: self
  type(mpas_geom),                   intent(inout) :: geom !< mpas mesh descriptors
  class(mpas_fields),                intent(in)    :: dxm  !< model increment fields
  class(mpas_fields),                intent(inout) :: dxg  !< increment containing linear geovar fields

  ! pool-related pointers
  type(mpas_pool_type), pointer :: mFields_tl, gFields_tl
  type(mpas_pool_data_type), pointer :: gdata

  ! reusable arrays
  real(kind=RKIND), dimension(:), pointer :: ptrr1_a
  real(kind=RKIND), dimension(:,:), pointer :: ptrr2_a, ptrr2_b, traj_ptrr2_a, traj_ptrr2_b
  real(kind=RKIND), dimension(:,:), allocatable :: r2, trajr2

  ! iteration-specific variables
  character(len=MAXVARLEN) :: geovar
  integer :: nCells, nVertLevels, nVertLevelsP1
  integer :: iVar

  ! air pressure on w levels
  real(kind=RKIND), allocatable :: plevels(:,:)

  ! convenient local variables
  mFields_tl => dxm % subFields
  gFields_tl => dxg % subFields
  nCells = geom%nCellsSolve
  nVertLevels = geom%nVertLevels
  nVertLevelsP1 = geom%nVertLevelsP1

  ! pre-calculate pressure on w levels
  allocate(plevels(1:nVertLevelsP1,1:nCells))
  call mpas_pool_get_array(self%trajectory, 'air_pressure', ptrr2_a)
  call mpas_pool_get_array(self%trajectory, 'air_pressure_at_surface', ptrr1_a)
  call pressure_half_to_full(ptrr2_a(:,1:nCells), geom%zgrid(:,1:nCells), ptrr1_a(1:nCells), &
                             nCells, nVertLevels, plevels)

  ! populate dxg
  do iVar = 1, dxg % nf

    geovar = trim(dxg % fldnames(iVar))

    if (geom%has_identity(geovar)) then
      ! "identity" variable changes are controlled in an external configuration file
      ! through the geom object. Refer to the mpas_geom type for more information.
      if (dxm%has(geom%identity(geovar))) then
        call dxm%copy_to(geom%identity(geovar), dxg, geovar)
      else if (is_hydro_base(geom%identity(geovar)) .and. &
               dxm%has(hydro_pt_name(geom%identity(geovar)))) then
        ! hydrometeor analyzed in power-transformed space: dq = (dq/dz)|traj * dz
        call dxg%get(geovar, gdata)
        call hydro_q_tl(geom%identity(geovar), gdata%r2%array(:,1:nCells))
      else
        !note: this warning is specifically for pressure => air_pressure in GNSSRO
        !write(message,'(2A)') &
        !   'WARNING: mpasjedi_lvc_model2geovars::multiply: '&
        !  &'increment missing identity field for geovar => ', trim(geovar)
        !call fckit_log%warning(message)
        continue
      end if
    else

      call dxg%get(geovar, gdata)

      select case (trim(geovar))

        case ( var_tv ) !-virtual_temperature
          ! get TL variables
          call dxm%get('air_temperature', ptrr2_a)
          call dxm%get('water_vapor_mixing_ratio_wrt_moist_air', ptrr2_b)

          ! temporary work fields
          allocate(r2(1:nVertLevels,1:nCells))     ! TL mixing_ratio
          allocate(trajr2(1:nVertLevels,1:nCells)) ! NL mixing_ratio

          ! get linearization state
          call mpas_pool_get_array(self%trajectory, 'air_temperature', traj_ptrr2_a)
          call mpas_pool_get_array(self%trajectory, 'water_vapor_mixing_ratio_wrt_moist_air', traj_ptrr2_b)
          call q_to_w(traj_ptrr2_b(:,1:nCells), trajr2(:,1:nCells)) !NL coeff.

          ! calculations
          call q_to_w_tl(ptrr2_b(:,1:nCells), traj_ptrr2_b(:,1:nCells), r2(:,1:nCells))
          call tw_to_tv_tl(ptrr2_a(:,1:nCells), r2(:,1:nCells), &
                           traj_ptrr2_a(:,1:nCells), trajr2(:,1:nCells), &
                           gdata%r2%array(:,1:nCells))

          ! cleanup
          deallocate(r2, trajr2)

        case ( var_mixr ) !-water_vapor_mixing_ratio_wrt_dry_air
          ! get TL variables
          call dxm%get('water_vapor_mixing_ratio_wrt_moist_air', ptrr2_a)

          ! get linearization state
          call mpas_pool_get_array(self%trajectory, 'water_vapor_mixing_ratio_wrt_moist_air', traj_ptrr2_a)

          ! calculations
          call q_to_w_tl(ptrr2_a(:,1:nCells), traj_ptrr2_a(:,1:nCells), gdata%r2%array(:,1:nCells))

          ! Ensure positive-definite mixing ratios
          !  with respect to precision of crtm::CRTM_Parameters::ZERO.
          ! TODO: this should be moved to the mpas_fields%read step
          !       only other place it should be needed is add_incr (already there)
          where (traj_ptrr2_a(:,1:nCells) <= MPAS_JEDI_ZERO_kr)
            gdata%r2%array(:,1:nCells) = MPAS_JEDI_ZERO_kr
          end where

        case ( var_clw_wp ) !-mass_content_of_cloud_liquid_water_in_atmosphere_layer
          call hydro_wp_tl('cloud_liquid_water', gdata%r2)

        case ( var_cli_wp ) !-mass_content_of_cloud_ice_in_atmosphere_layer
          call hydro_wp_tl('cloud_liquid_ice', gdata%r2)

        case ( var_clr_wp ) !-mass_content_of_rain_in_atmosphere_layer
          call hydro_wp_tl('rain_water', gdata%r2)

        case ( var_cls_wp ) !-mass_content_of_snow_in_atmosphere_layer
          call hydro_wp_tl('snow_water', gdata%r2)

        case ( var_clg_wp ) !-mass_content_of_graupel_in_atmosphere_layer
          call hydro_wp_tl('graupel', gdata%r2)

        case ( var_clh_wp ) !-mass_content_of_hail_in_atmosphere_layer
          call hydro_wp_tl('hail', gdata%r2)

      end select

    end if

  end do !iVar

  deallocate(plevels)

contains

  !> TL hydrometeor mixing-ratio increment dq [kg/kg], taken from the mixing-ratio
  !> increment itself or mapped from its power-transformed counterpart <name>_pt
  subroutine hydro_q_tl(qName, dq)
    character(len=*), intent(in)  :: qName
    real(kind=RKIND), intent(out) :: dq(:,:)   ! (nVertLevels, nCells)
    real(kind=RKIND), dimension(:,:), pointer :: dqm, dz, qTraj
    if (dxm%has(qName)) then
      call dxm%get(qName, dqm)
      dq = dqm(:,1:nCells)
    else if (dxm%has(hydro_pt_name(qName))) then
      call dxm%get(hydro_pt_name(qName), dz)
      call mpas_pool_get_array(self%trajectory, trim(qName), qTraj)
      dq = pt_slope(qTraj(:,1:nCells), geom%hydro_pt_exponent, geom%hydro_pt_qmin) * dz(:,1:nCells)
    else
      call abor1_ftn('mpasjedi_lvc_model2geovars::multiply: increment has neither '// &
                     trim(qName)//' nor '//trim(hydro_pt_name(qName)))
    end if
  end subroutine hydro_q_tl

  !> TL layer mass content [kg/m2] of a hydrometeor
  subroutine hydro_wp_tl(qName, qGeo_tl)
    character(len=*),           intent(in)    :: qName
    type(field2DReal), pointer, intent(inout) :: qGeo_tl
    real(kind=RKIND), allocatable :: dq(:,:)
    if (dxm%has(qName)) then
      ! unchanged behavior when the mixing ratio itself is analyzed
      call q_fields_TL(qName, mFields_tl, qGeo_tl, plevels, nCells, nVertLevels)
    else
      allocate(dq(nVertLevels, nCells))
      call hydro_q_tl(qName, dq)
      qGeo_tl%array(1:nVertLevels,1:nCells) = dq * &
        ( plevels(1:nVertLevels,1:nCells) - plevels(2:nVertLevels+1,1:nCells) ) / gravity
      deallocate(dq)
    end if
  end subroutine hydro_wp_tl

end subroutine multiply

! --------------------------------------------------------------------------------------------------

subroutine multiplyadjoint(self, geom, dxg, dxm)
  class(mpasjedi_lvc_model2geovars), intent(inout) :: self
  type(mpas_geom),                   intent(inout) :: geom
  class(mpas_fields),                intent(in)    :: dxg
  class(mpas_fields),                intent(inout) :: dxm

  ! pool-related pointers
  type(mpas_pool_type), pointer :: mFields_ad, gFields_ad
  type(mpas_pool_data_type), pointer :: gdata

  ! reusable array pointers
  real(kind=RKIND), dimension(:), pointer :: ptrr1_a
  real(kind=RKIND), dimension(:,:), pointer :: ptrr2_a, ptrr2_b, traj_ptrr2_a, traj_ptrr2_b
  real(kind=RKIND), dimension(:,:), allocatable :: r2, trajr2

  ! iteration-specific variables
  character(len=MAXVARLEN) :: geovar
  integer :: nCells, nVertLevels, nVertLevelsP1
  integer :: iVar

  ! air pressure on w levels
  real(kind=RKIND), allocatable :: plevels(:,:)

  ! convenient local variables
  mFields_ad => dxm % subFields
  gFields_ad => dxg % subFields
  nCells = geom%nCellsSolve
  nVertLevels = geom%nVertLevels
  nVertLevelsP1 = geom%nVertLevelsP1

  ! pre-calculate pressure on w levels
  allocate(plevels(1:nVertLevelsP1,1:nCells))
  call mpas_pool_get_array(self%trajectory, 'air_pressure', ptrr2_a)
  call mpas_pool_get_array(self%trajectory, 'air_pressure_at_surface', ptrr1_a)
  call pressure_half_to_full(ptrr2_a(:,1:nCells), geom%zgrid(:,1:nCells), ptrr1_a(1:nCells), &
                             nCells, nVertLevels, plevels)

  ! populate dxg
  do iVar = 1, dxg % nf

    geovar = trim(dxg % fldnames(iVar))

    if (geom%has_identity(geovar)) then
      ! "identity" variable changes are controlled in an external configuration file
      ! through the geom object. Refer to the mpas_geom type for more information.
      if (dxm%has(geom%identity(geovar))) then
        call dxm%copy_to_ad(geom%identity(geovar), dxg, geovar)
      else if (is_hydro_base(geom%identity(geovar)) .and. &
               dxm%has(hydro_pt_name(geom%identity(geovar)))) then
        call dxg%get(geovar, gdata)
        call hydro_q_ad(geom%identity(geovar), gdata%r2%array(:,1:nCells))
      else
        !note: this warning is specifically for pressure => air_pressure in GNSSRO
        !write(message,'(2A)') &
        !   'WARNING: mpasjedi_lvc_model2geovars::multiplyadjoint: '&
        !  &'increment missing identity field for geovar => ', trim(geovar)
        !call fckit_log%warning(message)
        continue
      end if
    else

      call dxg%get(geovar, gdata)

      select case (trim(geovar))

        case ( var_tv ) !-virtual_temperature
          ! get AD variables
          call dxm%get('air_temperature', ptrr2_a)
          call dxm%get('water_vapor_mixing_ratio_wrt_moist_air', ptrr2_b)

          ! temporary work fields
          allocate(r2(1:nVertLevels,1:nCells))     ! AD of mixing_ratio
          allocate(trajr2(1:nVertLevels,1:nCells)) ! NL mixing_ratio

          ! get linearization state
          call mpas_pool_get_array(self%trajectory, 'air_temperature', traj_ptrr2_a)
          call mpas_pool_get_array(self%trajectory, 'water_vapor_mixing_ratio_wrt_moist_air', traj_ptrr2_b)
          call q_to_w(traj_ptrr2_b(:,1:nCells), trajr2(:,1:nCells)) !NL coeff.

          ! calculations
          r2(:,1:nCells) = MPAS_JEDI_ZERO_kr !initialize local var.
          call tw_to_tv_ad(ptrr2_a(:,1:nCells), r2(:,1:nCells), &
                           traj_ptrr2_a(:,1:nCells), trajr2(:,1:nCells), &
                           gdata%r2%array(:,1:nCells) )
          call q_to_w_ad(ptrr2_b(:,1:nCells), traj_ptrr2_b(:,1:nCells), r2(:,1:nCells))

          ! cleanup
          deallocate(r2, trajr2)

        case ( var_mixr ) !-water_vapor_mixing_ratio_wrt_dry_air
          ! get AD variables
          call dxm%get('water_vapor_mixing_ratio_wrt_moist_air', ptrr2_a)

          ! temporary work fields
          allocate(r2(1:nVertLevels,1:nCells)) ! AD of var_mixr

          ! get linearization state
          call mpas_pool_get_array(self%trajectory, 'water_vapor_mixing_ratio_wrt_moist_air', traj_ptrr2_a)

          ! Ensure positive-definite mixing ratios
          !  with respect to precision of crtm::CRTM_Parameters::ZERO.
          ! TODO: this should be moved to the mpas_fields%read step
          !       only other place it should be needed is add_incr (already there)
          r2 = gdata%r2%array(1:nVertLevels,1:nCells)
          where (traj_ptrr2_a(:,1:nCells) <= MPAS_JEDI_ZERO_kr)
            r2(:,1:nCells) = MPAS_JEDI_ZERO_kr
          end where

          ! calculations
          call q_to_w_ad(ptrr2_a(:,1:nCells), traj_ptrr2_a(:,1:nCells), r2(:,1:nCells))

          ! cleanup
          deallocate(r2)

        case ( var_clw_wp ) !-mass_content_of_cloud_liquid_water_in_atmosphere_layer
          call hydro_wp_ad('cloud_liquid_water', gdata%r2)

        case ( var_cli_wp ) !-mass_content_of_cloud_ice_in_atmosphere_layer
          call hydro_wp_ad('cloud_liquid_ice', gdata%r2)

        case ( var_clr_wp ) !-mass_content_of_rain_in_atmosphere_layer
          call hydro_wp_ad('rain_water', gdata%r2)

        case ( var_cls_wp ) !-mass_content_of_snow_in_atmosphere_layer
          call hydro_wp_ad('snow_water', gdata%r2)

        case ( var_clg_wp ) !-mass_content_of_graupel_in_atmosphere_layer
          call hydro_wp_ad('graupel', gdata%r2)

        case ( var_clh_wp ) !-mass_content_of_hail_in_atmosphere_layer
          call hydro_wp_ad('hail', gdata%r2)

      end select

    end if

  end do !iVar

  deallocate(plevels)

contains

  !> AD of the power-transformed path of hydro_q_tl: dz_ad += (dq/dz)|traj * dq_ad
  subroutine hydro_q_ad(qName, dq_ad)
    character(len=*), intent(in) :: qName
    real(kind=RKIND), intent(in) :: dq_ad(:,:)   ! (nVertLevels, nCells)
    real(kind=RKIND), dimension(:,:), pointer :: dz_ad, qTraj
    call dxm%get(hydro_pt_name(qName), dz_ad)
    call mpas_pool_get_array(self%trajectory, trim(qName), qTraj)
    dz_ad(:,1:nCells) = dz_ad(:,1:nCells) + &
      pt_slope(qTraj(:,1:nCells), geom%hydro_pt_exponent, geom%hydro_pt_qmin) * dq_ad
  end subroutine hydro_q_ad

  !> AD of hydro_wp_tl
  subroutine hydro_wp_ad(qName, qGeo_ad)
    character(len=*),           intent(in) :: qName
    type(field2DReal), pointer, intent(in) :: qGeo_ad
    real(kind=RKIND), allocatable :: dq_ad(:,:)
    if (dxm%has(qName)) then
      call q_fields_AD(qName, mFields_ad, qGeo_ad, plevels, nCells, nVertLevels)
    else if (dxm%has(hydro_pt_name(qName))) then
      allocate(dq_ad(nVertLevels, nCells))
      dq_ad = qGeo_ad%array(1:nVertLevels,1:nCells) * &
        ( plevels(1:nVertLevels,1:nCells) - plevels(2:nVertLevels+1,1:nCells) ) / gravity
      call hydro_q_ad(qName, dq_ad)
      deallocate(dq_ad)
    else
      call abor1_ftn('mpasjedi_lvc_model2geovars::multiplyadjoint: increment has neither '// &
                     trim(qName)//' nor '//trim(hydro_pt_name(qName)))
    end if
  end subroutine hydro_wp_ad

end subroutine multiplyadjoint

! --------------------------------------------------------------------------------------------------

end module mpasjedi_lvc_model2geovars_mod
! --------------------------------------------------------------------------------------------------
