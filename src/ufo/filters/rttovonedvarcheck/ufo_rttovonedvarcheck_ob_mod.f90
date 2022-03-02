! (C) copyright 2020 Met Office
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.

!> Fortran module which contains the observation metadata for a single observation

module ufo_rttovonedvarcheck_ob_mod

use kinds
use missing_values_mod
use ufo_constants_mod, only: zero
use ufo_rttovonedvarcheck_constants_mod
use ufo_rttovonedvarcheck_setup_mod
use ufo_rttovonedvarcheck_pcemis_mod

implicit none
private

! Ob info type definition
type, public :: ufo_rttovonedvarcheck_ob

  character(len=max_string) :: forward_mod_name !< forward model name (RTTOV only one at the moment)
  integer              :: nlocs !< number of locations = 1
  integer              :: surface_type  !< surface type of observation
  integer              :: niter
  integer, allocatable :: channels_used(:) !< channels used for this observation
  integer, allocatable :: channels_all(:) !< all channels used for output
  integer, allocatable :: rejected_channels_ctp(:) !< a list of channels rejected based on the retrieved ctp
  real(kind_real)      :: latitude !< latitude of observation
  real(kind_real)      :: longitude !< longitude of observation
  real(kind_real)      :: elevation  !< elevation above sea level of observation
  real(kind_real)      :: sensor_zenith_angle  !< sensor zenith of observation
  real(kind_real)      :: sensor_azimuth_angle  !< sensor azimuth of observation
  real(kind_real)      :: solar_zenith_angle !< solar zenith of observation
  real(kind_real)      :: solar_azimuth_angle  !< solar azimuth of observation
  real(kind_real)      :: cloudtopp !< cloud top pressure (used in if cloudy retrieval used)
  real(kind_real)      :: cloudfrac !< cloud fraction (used in if cloudy retrieval used)
  real(kind_real)      :: final_cost !< final cost at solution
  real(kind_real)      :: LWP !< retrieved liquid water path. This is output for future filters
  real(kind_real), allocatable :: clw(:) !< cloud liquid water profile. Currently used in Var 
  real(kind_real), allocatable :: yobs(:) !< satellite BTs
  real(kind_real), allocatable :: final_bt_diff(:) !< final bt difference if converged
  real(kind_real), allocatable :: emiss(:) !< surface emissivity
  real(kind_real), allocatable :: background_T(:) !< background temperature used by qsplit
  real(kind_real), allocatable :: output_profile(:) !< retrieved state at converge as profile vector
  real(kind_real), allocatable :: output_BT(:) !< Brightness temperatures using retrieved state
  real(kind_real), allocatable :: background_BT(:) !< Brightness temperatures from 1st itreration
  real(kind_real), allocatable :: pcemiss(:) !< principal component emissivity
  logical              :: retrievecloud  !< flag to turn on retrieve cloud
  logical              :: mwscatt !< flag to use rttov-scatt model through the interface
  logical              :: mwscatt_totalice !< flag to use total ice (rather then ciw) for rttov-scatt simulations
  logical, allocatable :: calc_emiss(:) !< flag to decide if RTTOV calculates emissivity
  type(ufo_rttovonedvarcheck_pcemis), pointer :: pcemiss_object !< pointer to the IR pc emissivity object

contains
  procedure :: setup  => ufo_rttovonedvarcheck_InitOb
  procedure :: delete => ufo_rttovonedvarcheck_DeleteOb
  procedure :: info => ufo_rttovonedvarcheck_PrintOb

end type ufo_rttovonedvarcheck_ob

contains

!-------------------------------------------------------------------------------
!> Initialize observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_InitOb(self, & ! out
                                        nchans, &  ! in
                                        nlevels, & ! in
                                        nprofelements, & ! in
                                        nchans_all, & ! in
                                        storeclw)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_ob), intent(out) :: self !< observation metadata type
integer, intent(in) :: nchans !< number of channels used for this particular observation
integer, intent(in) :: nlevels !< number of levels in the profile
integer, intent(in) :: nprofelements !< number of profile elements used
integer :: nchans_all !< Size of all channels in ObsSpace
logical, intent(in) :: storeclw
character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_InitOb"
real(kind_real) :: missing

missing = missing_value(missing)

call self % delete()

allocate(self % yobs(nchans))
allocate(self % final_bt_diff(nchans))
allocate(self % channels_used(nchans))
allocate(self % channels_all(nchans_all))
allocate(self % emiss(nchans_all))
allocate(self % background_T(nlevels))
allocate(self % output_profile(nprofelements))
allocate(self % output_BT(nchans_all))
allocate(self % background_BT(nchans_all))
allocate(self % calc_emiss(nchans_all))
if (storeclw) then
  allocate(self % clw(nlevels))
  self % clw(:) = missing
endif
self % yobs(:) = missing
self % final_bt_diff(:) = missing
self % emiss(:) = missing
self % background_T(:) = missing
self % output_profile(:) = missing
self % output_BT(:) = missing
self % background_BT(:) = missing
self % calc_emiss(:) = .true.

end subroutine ufo_rttovonedvarcheck_InitOb

!------------------------------------------------------------------------------
!> Delete the single observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_DeleteOb(self)    ! inout

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_ob), intent(inout) :: self !< observation metadata type

character(len=*), parameter :: routinename = "ufo_rttovonedvarcheck_DeleteOb"
real(kind_real) :: missing

missing = missing_value(missing)

self % nlocs = 1
self % latitude = missing
self % longitude = missing
self % elevation = missing
self % surface_type = 0
self % niter = 0
self % sensor_zenith_angle = missing
self % sensor_azimuth_angle = missing
self % solar_zenith_angle = missing
self % solar_azimuth_angle = missing
self % cloudtopp = 500.0_kind_real
self % cloudfrac = zero
self % final_cost = missing
self % LWP = missing
self % retrievecloud = .false.
self % mwscatt = .false.
self % mwscatt_totalice = .false.

if (allocated(self % yobs))           deallocate(self % yobs)
if (allocated(self % final_bt_diff))  deallocate(self % final_bt_diff)
if (allocated(self % channels_used))  deallocate(self % channels_used)
if (allocated(self % channels_all))   deallocate(self % channels_all)
if (allocated(self % emiss))          deallocate(self % emiss)
if (allocated(self % background_T))   deallocate(self % background_T)
if (allocated(self % output_profile)) deallocate(self % output_profile)
if (allocated(self % output_BT))      deallocate(self % output_BT)
if (allocated(self % background_BT))  deallocate(self % background_BT)
if (allocated(self % calc_emiss))     deallocate(self % calc_emiss)
if (allocated(self % rejected_channels_ctp)) deallocate(self % rejected_channels_ctp)
if (allocated(self % clw))            deallocate(self % clw)
if (allocated(self % pcemiss))        deallocate(self % pcemiss)

self % pcemiss_object => null()

end subroutine ufo_rttovonedvarcheck_DeleteOb

!------------------------------------------------------------------------------
!> Print information about the single observation object
!!
!! \author Met Office
!!
!! \date 09/06/2020: Created
!!
subroutine ufo_rttovonedvarcheck_PrintOb(self)

implicit none

! subroutine arguments:
class(ufo_rttovonedvarcheck_ob), intent(inout) :: self !< observation metadata type

character(len=20) :: surface_type

if (self % surface_type == RTLand) then
  surface_type = "land"
else if (self % surface_type == RTSea) then
  surface_type = "sea"
else if (self % surface_type == RTIce) then
  surface_type = "ice"
else
  surface_type = "unknown"
end if

write(*,"(A,2F8.2)") "Lat,Long:",self % latitude, self % longitude
write(*,*) "Surface type for RTTOV: ",surface_type
write(*,"(A,F8.2)") "Surface height:",self % elevation
write(*,"(A,F8.2)") "Satellite zenith angle: ",self % sensor_zenith_angle
write(*,"(A,F8.2)") "Solar zenith angle: ",self % solar_zenith_angle
write(*,"(A,F8.2)") "Cloud Top Pressure",self % cloudtopp
write(*,"(A,F8.2)") "Cloud Fraction",self % cloudfrac
write(*,"(A)") "Background T profile: "
write(*,"(10F8.2)") self % background_T
write(*,"(A)") "Emissivity: "
write(*,"(10F8.2)") self % emiss(:)
write(*,"(A)") "Emissivity PC : "
write(*,"(10F18.8)") self % pcemiss(:)

end subroutine

!-------------------------------------------------------------------------------

end module ufo_rttovonedvarcheck_ob_mod
