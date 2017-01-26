module cloud_module
  use amr_parameters
  use hydro_parameters,only:Msun
  
  real(dp)::bl_fac=1.   !multiply calculated boxlen by this factor

  !Initial conditions parameters for the dense core
  logical ::bb_test=.false. ! Activate Boss & Bodenheimer inital conditions instead of 1/R^2 density profile
  logical ::uniform_bmag=.false. ! Activate uniform magnetic field initial conditions for BE-like initial density profile
  real(dp)::mass_c=1.         !cloud mass in solar mass
  real(dp)::contrast=100.d0   !density contrast (used when bb_test=.true.)
  real(dp)::cont=1.           !density contrast (used when bb_test=.false.)
  real(dp)::rap=1.            !axis ratio
  real(dp)::ff_sct=1.         !freefall time / sound crossing time
  real(dp)::ff_rt=1.          !freefall time / rotation time
  real(dp)::ff_act=1.         !freefall time / Alfven crossing time
  real(dp)::ff_vct=1.         !freefall time / Vrms crossing time
  real(dp)::theta_mag=0.      !angle between magnetic field and rotation axis

  real(dp):: C2_vis=0.0d0 !Von Neumann & Richtmeyer artificial viscosity coefficient 3 en principe
  real(dp):: alpha_dense_core=0.5d0
  real(dp):: beta_dense_core=0.0d0
  real(dp):: crit=0.0d0
  real(dp):: delta_rho=0.0d0
  real(dp):: Mach=0.0d0

  ! PMS evolution related stuff
  logical :: rt_feedback=.false.       ! take into account RT feedback
  logical :: PMS_evol=.false.          ! Take into account PMS evolution subgrid model
  logical :: Hosokawa_track=.false.    ! Take into account PMS evolution subgrid model
  real(dp):: dt_lsink_update=50        ! frequency of the sink luminosity update with PMS evolution (in yr)
  real(dp):: epsilonlib=0.0            ! Fraction of energy absorbed by the prostostars at the accretion shock
  real(dp):: mprotostar=0.0009546*Msun ! initial mass of the protostar (1 Mjup)
  real(dp):: rstar_init=2.5            ! Initial radius of the protostar in Rsun
  integer :: modell=0
  integer :: modrestart=0              ! name of model you want to restart from, this is an input
  real(dp):: facc_star_lum=0.75d0      ! fraction of the accretion luminosity radiated by the sinks
  real(dp):: facc_star=0.5d0           ! fraction of the sink accreted mass actually accreted by the star
  integer::nmdot_PMS,nm_PMS,ndata_PMS
  integer ,allocatable,dimension(:)::nb_ligne_PMS
  real(dp),allocatable,dimension(:,:,:)::data_PMS

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
  namelist/cloud_params/bl_fac
!  namelist/cloud_params/alpha_dense_core,beta_dense_core,crit,delta_rho &
!       & ,mass_c,rap,cont,ff_sct,ff_rt,ff_act,ff_vct,theta_mag,bb_test &
!       & ,contrast,Mach,uniform_bmag

  ! Read namelist file
  rewind(1)
  read(1,NML=cloud_params,END=101)
101 continue                                   ! No harm if no namelist

  ! Get some units out there
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)



end subroutine read_cloud_params

