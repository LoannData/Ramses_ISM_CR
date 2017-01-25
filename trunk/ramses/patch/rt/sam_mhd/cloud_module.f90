module cloud_module
  use amr_parameters


  !Initial conditions parameter for the dense core
  real(dp)::mass_c=0d0   !cloud mass in solar mass
  real(dp)::rap=1.      !axis ratio
  real(dp)::cont=1.     !density contras
  real(dp)::ff_sct=1.   !freefall time / sound crossing time
  real(dp)::ff_rt=1.    !freefall time / rotation time
  real(dp)::ff_act=1.   !freefall time / Alfven crossing time
  real(dp)::ff_vct=1.   !freefall time / Vrms crossing time
  real(dp)::thet_mag=0. !angle between magnetic field and rotation axis
  real(dp)::bl_fac=1.   !multiply calculated boxlen by this factor


end module cloud_module

subroutine read_cloud_params(nml_ok)

  use amr_parameters
  use feedback_module
  use clfind_commons
  use cloud_module

  implicit none
  logical::nml_ok
  real(dp)::cellsize
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),parameter::pcincm=3.086d18

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/cloud_params/mass_c,rap,cont,ff_sct,ff_rt,ff_act,ff_vct,thet_mag &
       & ,bl_fac

  ! Read namelist file
  rewind(1)
  read(1,NML=cloud_params,END=101)
101 continue                                   ! No harm if no namelist

  ! Get some units out there
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)



end subroutine read_cloud_params

