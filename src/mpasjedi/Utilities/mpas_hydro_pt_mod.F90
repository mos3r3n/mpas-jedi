! (C) Copyright 2026 UCAR
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

!> \brief Power-transformed (Box-Cox) hydrometeor mixing ratios
!> \details
!>  Provides the forward, inverse, and linearized transforms used to analyze
!>  hydrometeor mixing ratios q [kg/kg] in a power-transformed space
!>
!>     z = ( max(q, qmin)**p - 1 ) / p           (forward, q -> z)
!>     q = max( p*z + 1, 0 )**(1/p)              (inverse, z -> q)
!>     dq/dz = max(q, qmin)**(1-p)               (Jacobian at the trajectory)
!>
!>  For each hydrometeor field "<name>" (e.g., rain_water), the transformed
!>  field is named "<name>_pt" (e.g., rain_water_pt). The transformed fields are
!>  diagnosed from q whenever model fields are read (update_diagnostic_fields),
!>  analyzed directly (they can be control/analysis/ensemble variables), and
!>  mapped back to q nonlinearly in State::add_incr.
!>
!>  p and qmin are set via the Geometry YAML keys
!>     "hydrometeor power transform exponent"  (default 0.4)
!>     "hydrometeor power transform qmin"      (default 1.0e-8 kg/kg)
module mpas_hydro_pt_mod

use mpas_derived_types
use mpas_kind_types, only: RKIND
use mpas_pool_routines
use ufo_vars_mod, only: MAXVARLEN

implicit none

private
public :: nHydroPT, hydroBaseNames, hydro_pt_name, is_hydro_base, pool_has_field, &
          q_to_pt, pt_to_q, pt_slope, &
          hydro_pt_diagnose, hydro_pt_add_incr

!> Both naming families used in geovars.yaml are supported; use the one that
!> matches the rest of your configuration (state, ensemble, and obs geovars).
integer, parameter :: nHydroPT = 12
character(len=MAXVARLEN), parameter :: hydroBaseNames(nHydroPT) = &
   [ character(len=MAXVARLEN) :: &
     'cloud_liquid_water', 'cloud_liquid_ice', 'rain_water', 'snow_water', 'graupel', 'hail', &
     'mixing_ratio_of_cloud_liquid_water', 'mixing_ratio_of_cloud_ice', 'mixing_ratio_of_rain', &
     'mixing_ratio_of_snow', 'mixing_ratio_of_graupel', 'mixing_ratio_of_hail' ]
character(len=*), parameter :: ptSuffix = '_pt'

contains

! ------------------------------------------------------------------------------

pure function hydro_pt_name(baseName) result(ptName)
   character(len=*), intent(in) :: baseName
   character(len=MAXVARLEN) :: ptName
   ptName = trim(baseName)//ptSuffix
end function hydro_pt_name

! ------------------------------------------------------------------------------

pure logical function is_hydro_base(fieldName)
   character(len=*), intent(in) :: fieldName
   is_hydro_base = any(trim(fieldName) == hydroBaseNames)
end function is_hydro_base

! ------------------------------------------------------------------------------

!> Membership test for a pool. Do not call while iterating over the same pool.
logical function pool_has_field(pool, fieldName)
   type(mpas_pool_type), pointer, intent(in) :: pool
   character(len=*),              intent(in) :: fieldName
   type(mpas_pool_iterator_type) :: itr
   pool_has_field = .false.
   call mpas_pool_begin_iteration(pool)
   do while ( mpas_pool_get_next_member(pool, itr) )
      if ( itr % memberType == MPAS_POOL_FIELD .and. &
           trim(itr % memberName) == trim(fieldName) ) pool_has_field = .true.
   end do
end function pool_has_field

! ------------------------------------------------------------------------------

elemental function q_to_pt(q, p, qmin) result(z)
   real(kind=RKIND), intent(in) :: q, p, qmin
   real(kind=RKIND) :: z
   z = ( max(q, qmin)**p - 1.0_RKIND ) / p
end function q_to_pt

! ------------------------------------------------------------------------------

elemental function pt_to_q(z, p) result(q)
   real(kind=RKIND), intent(in) :: z, p
   real(kind=RKIND) :: q
   q = max(p*z + 1.0_RKIND, 0.0_RKIND)**(1.0_RKIND/p)
end function pt_to_q

! ------------------------------------------------------------------------------

!> dq/dz of the inverse transform, evaluated at the linearization state qTraj
elemental function pt_slope(qTraj, p, qmin) result(s)
   real(kind=RKIND), intent(in) :: qTraj, p, qmin
   real(kind=RKIND) :: s
   s = max(qTraj, qmin)**(1.0_RKIND - p)
end function pt_slope

! ------------------------------------------------------------------------------

!> Diagnose <name>_pt from <name> for every hydrometeor present in both forms
subroutine hydro_pt_diagnose(fields, ngrid, p, qmin)
   type(mpas_pool_type), pointer, intent(inout) :: fields
   integer,                       intent(in)    :: ngrid
   real(kind=RKIND),              intent(in)    :: p, qmin

   real(kind=RKIND), dimension(:,:), pointer :: q, z
   integer :: ih

   do ih = 1, nHydroPT
      if ( .not. pool_has_field(fields, hydroBaseNames(ih)) ) cycle
      if ( .not. pool_has_field(fields, hydro_pt_name(hydroBaseNames(ih))) ) cycle
      call mpas_pool_get_array(fields, trim(hydroBaseNames(ih)), q)
      call mpas_pool_get_array(fields, trim(hydro_pt_name(hydroBaseNames(ih))), z)
      z(:,1:ngrid) = q_to_pt(q(:,1:ngrid), p, qmin)
   end do

end subroutine hydro_pt_diagnose

! ------------------------------------------------------------------------------

!> Nonlinear update of q from a transformed-space increment dz.
!> Must be called AFTER the generic state+increment addition, i.e., when the
!> state's <name>_pt already holds z_b + dz.
!>
!>   q_a = max(0, q_b + Finv(F(q_b) + dz) - Finv(F(q_b)))
!>
!> Using the difference of inverses (rather than Finv(z_a) alone) leaves q_b
!> unchanged where dz = 0, including points with 0 <= q_b < qmin.
!> Afterwards <name>_pt is re-diagnosed from q_a for consistency.
subroutine hydro_pt_add_incr(stateFields, incrFields, ngrid, p, qmin)
   type(mpas_pool_type), pointer, intent(inout) :: stateFields
   type(mpas_pool_type), pointer, intent(in)    :: incrFields
   integer,                       intent(in)    :: ngrid
   real(kind=RKIND),              intent(in)    :: p, qmin

   real(kind=RKIND), dimension(:,:), pointer :: q, z, dz
   real(kind=RKIND), dimension(:,:), allocatable :: zb
   character(len=MAXVARLEN) :: qName, zName
   integer :: ih

   do ih = 1, nHydroPT
      qName = hydroBaseNames(ih)
      zName = hydro_pt_name(qName)
      if ( .not. pool_has_field(incrFields, zName) ) cycle
      if ( pool_has_field(incrFields, qName) ) then
         call abor1_ftn('hydro_pt_add_incr: increment contains both '//trim(qName)// &
                        ' and '//trim(zName)//'; analyze only one of them')
      end if
      if ( .not. pool_has_field(stateFields, qName) .or. &
           .not. pool_has_field(stateFields, zName) ) then
         call abor1_ftn('hydro_pt_add_incr: state must contain both '//trim(qName)// &
                        ' and '//trim(zName)//' (add both to the state variables)')
      end if

      call mpas_pool_get_array(stateFields, trim(qName), q)
      call mpas_pool_get_array(stateFields, trim(zName), z)
      call mpas_pool_get_array(incrFields,  trim(zName), dz)

      allocate(zb(size(q,1), ngrid))
      zb = q_to_pt(q(:,1:ngrid), p, qmin)
      q(:,1:ngrid) = max( 0.0_RKIND, &
                          q(:,1:ngrid) + pt_to_q(zb + dz(:,1:ngrid), p) - pt_to_q(zb, p) )
      z(:,1:ngrid) = q_to_pt(q(:,1:ngrid), p, qmin)
      deallocate(zb)
   end do

end subroutine hydro_pt_add_incr

! ------------------------------------------------------------------------------

end module mpas_hydro_pt_mod
