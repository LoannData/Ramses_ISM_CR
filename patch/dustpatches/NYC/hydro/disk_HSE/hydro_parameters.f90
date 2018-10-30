module hydro_parameters

#ifdef grackle
  use grackle_parameters
#endif
  use amr_parameters

  ! Number of independant variables
#ifndef NENER
  integer,parameter::nener=0
#else
  integer,parameter::nener=NENER
#endif
#ifndef NVAR
  integer,parameter::nvar=ndim+2+nener
#else
  integer,parameter::nvar=NVAR
#endif
  ! Size of hydro kernel
  integer,parameter::iu1=-1
  integer,parameter::iu2=+4
  integer,parameter::ju1=(1-ndim/2)-1*(ndim/2)
  integer,parameter::ju2=(1-ndim/2)+4*(ndim/2)
  integer,parameter::ku1=(1-ndim/3)-1*(ndim/3)
  integer,parameter::ku2=(1-ndim/3)+4*(ndim/3)
  integer,parameter::if1=1
  integer,parameter::if2=3
  integer,parameter::jf1=1
  integer,parameter::jf2=(1-ndim/2)+3*(ndim/2)
  integer,parameter::kf1=1
  integer,parameter::kf2=(1-ndim/3)+3*(ndim/3)

  ! Imposed boundary condition variables
  real(dp),dimension(1:MAXBOUND,1:nvar)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0
#if NENER>0
  real(dp),dimension(1:MAXBOUND,1:NENER)::prad_bound=0.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:MAXBOUND,1:NVAR-NDIM-2-NENER)::var_bound=0.0
#endif
  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::floor_d=1.d-10   ! Density floor
  real(dp)::floor_u=1.d-10   ! Velocity floor
  real(dp)::floor_p=1.d-10   ! Pressure floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
#if NENER>0
  real(dp),dimension(1:NENER)::err_grad_prad=-1.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:NVAR-NDIM-2)::err_grad_var=-1.0
#endif
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Initial conditions hydro variables
  ! 1 = Gressel; 2 = Bai; 3 = NFW galaxy halo.
  integer::disk_setup=1

  ! Disk parameters: Center, char. radius, density, power indices,
  real(dp)::disk_xc=0.d0,disk_yc=0.d0,disk_zc=0.d0
  real(dp)::disk_R0=5.0d0
  real(dp)::disk_D0=1.0d0
  real(dp)::disk_rho_pow=-1.5d0
  real(dp)::disk_T0=100.0d0
  real(dp)::disk_cs_pow=-1.0d0
  real(dp)::disk_Mass=1.0d0 ! Solar units
  real(dp)::disk_abar=1.22d0
  real(dp)::disk_vrms=0.0d0 ! < 0 means fraction of Cs, Gressel = -0.01
  real(dp)::disk_V0=2.0
  real(dp)::disk_P0=1.0
  real(dp)::disk_Textra=1.05d0
  real(dp)::disk_aspect=0.05d0  ! Aspect ratio H/R
  real(dp)::disk_ratio=0.05d0   ! Used for inforcing refinement region, typically ~ 7xdisk_aspect
  real(dp)::disk_innerR=0.5d0   ! Transition radius to inner boundary
  real(dp)::disk_outerR=1.d30   ! Transition radius to outer boundary
  real(dp)::disk_outerR2=19.8d0 ! Transition radius to second outer boundary, i.e., HTL17
  real(dp)::disk_innerR0=0.1d0  ! Transition to zero gravity in inner boundary
  real(dp)::disk_outerR0=1.d30  !  Transition to zero gravity in outer boundary
  real(dp)::disk_vramp=1.0d0
  real(dp)::disk_shear_forcing=1.d0
  real(dp)::disk_raiseP_factor=0.d0
  real(dp)::disk_delta_alpha=1.d0
  integer ::disk_nsubsample=0
  real(dp)::disk_testR=-1.d0
  real(dp)::disk_testz=-1.d0
  real(dp)::disk_m=0.5d0
  real(dp)::disk_trans_fact=1.d0
  real(dp)::disk_exp_limit=20.0
  real(dp)::disk_R200=300.0
  real(dp)::disk_R1=6.0
  real(dp)::disk_H1=0.15
  real(dp)::disk_D1=-1.0
  real(dp)::disk_conc=13.0
  real(dp)::disk_thetaMax=0.5236
  real(dp)::grav_angle=-1.d0
  real(dp)::grav_width=0.3d0

! table for interpolating rotation velocity
  integer, parameter :: n_zint=2**20
  real(dp), dimension(n_zint) :: vrot_table = 0.d0, p_table, d_table, P0_table, P1_table, P2_table, P3_table
  logical :: disk_ReadTable=.true.
  logical :: disk_smooth_vertRho=.true.
  integer :: disk_TableLevel=-1

  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.
#if NENER>0
  real(dp),dimension(1:MAXREGION,1:NENER)::prad_region=0.0
#endif
#if NVAR>NDIM+2+NENER
  real(dp),dimension(1:MAXREGION,1:NVAR-NDIM-2-NENER)::var_region=0.0
#endif
  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  real(dp)::slope_theta=1.5d0
  real(dp)::gamma=1.4d0
  real(dp),dimension(1:512)::gamma_rad=1.33333333334d0
  real(dp)::courant_factor=0.5d0
  logical::apply_relaxation=.false.
  logical::relaxation_unew_uold=.false.
  real(dp)::relaxation_time=0.d0
  real(dp)::relaxation_OuterRadius=0.d0 ! > 0 spherical, < 0 cylindrical
  real(dp)::relaxation_InnerRadius=0.d0 ! > 0 spherical, < 0 cylindrical
  real(dp)::relax_factor=2.d0
  real(dp)::relax_revert=1.d0
  integer ::relaxation_step=0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1

  ! Passive variables index
  integer::imetal=6
  integer::idelay=6
  integer::ixion=6
  integer::ichem=6
  integer::ivirial1=6
  integer::ivirial2=6
  integer::inener=6

end module hydro_parameters
