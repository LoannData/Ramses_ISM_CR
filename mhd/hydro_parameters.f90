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
#ifndef NGRP
  integer,parameter::ngrp=0   ! Number of radiative energy groups
#else
  integer,parameter::ngrp=NGRP
#endif
#if USE_M_1==0
  integer,parameter::nrad=ngrp          ! Number of pure radiative variables (= radiative energies)
  integer,parameter::nvar_bicg=nrad     ! Number of variables in BICG (= radiative variables)
#endif
#if USE_M_1==1
  integer,parameter::nrad=(1+ndim)*ngrp  ! Number of pure radiative variables (= radiative energies + radiative fluxes)
  integer,parameter::nvar_bicg=nrad+1    ! Number of variables in BICG (= temperature + radiative variables)
#endif
  integer,parameter::nvar_trad=nrad+1   ! Total number of radiative variables (= temperature + radiative energies)

#ifndef NEXTINCT
  integer,parameter::nextinct = 0       ! Add a variable to store extinction coefficient [0,1]
#else
  integer,parameter::nextinct = NEXTINCT
#endif

  ! Advect internal energy as a passive scalar, in a supplementary index
#ifndef NPSCAL
  integer,parameter::npscal=1
#else
  integer,parameter::npscal=NPSCAL
#endif
! Cosmic rays energy groups
#ifndef NCR
  integer,parameter::ncr=0
#else
  integer,parameter::ncr=NCR
#endif

  integer,parameter::nent=nener-ngrp      ! Number of non-thermal energies
#if USE_M_1==0
  integer,parameter::nfr = 0              ! Number of radiative fluxes for M1
#else
  integer,parameter::nfr =ndim*ngrp       ! Number of radiative fluxes for M1
#endif

  ! First index of variables (in fact index just before the first index)
  ! so that we can loop over 1,nener for instance
  integer,parameter::firstindex_ent=8     ! for non-thermal energies
  integer,parameter::firstindex_er=8+nent ! for radiative energies
  integer,parameter::firstindex_fr=8+nener ! for radiative fluxes (if M1)
  integer,parameter::firstindex_extinct=8+nent+nrad ! for extinction
  integer,parameter::firstindex_pscal=8+nent+nrad+nextinct ! for passive scalars
  integer::lastindex_pscal ! last index for passive scalars other than internal energy
  ! Initialize NVAR
#ifndef NVAR
  integer,parameter::nvar=8+nent+nrad+nextinct+npscal
#else
  integer,parameter::nvar=NVAR
#endif

  ! EOS parameters
  integer  :: nRho,nEnergy,nTemp
  real(dp) :: rhomin,rhomax,Emax,emin,yHe,Tmax,Tmin

  ! barotrop parameters
  integer  :: nrho_barotrop
  logical  :: analytical_barotrop = .true.
  real(dp) :: rhomin_barotrop,rhomax_barotrop,drho_barotrop

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
  real(dp),dimension(1:MAXBOUND,1:nvar+3)::boundary_var
  real(dp),dimension(1:MAXBOUND)::d_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::p_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::u_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::v_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::w_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::A_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::B_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::C_bound=0.0d0
  real(dp),dimension(1:MAXBOUND)::T_bound=0.0d0
  real(dp),dimension(1:MAXBOUND,1:ngrp)::E_bound=0.0d0
  real(dp),dimension(1:MAXBOUND,1:ngrp)::fx_bound=0.0d0
  real(dp),dimension(1:MAXBOUND,1:ngrp)::fy_bound=0.0d0
  real(dp),dimension(1:MAXBOUND,1:ngrp)::fz_bound=0.0d0
#if NENER>0
  real(dp),dimension(1:MAXBOUND,1:NENER)::prad_bound=0.0
#endif
#if NPSCAL>0
  real(dp),dimension(1:MAXBOUND,1:npscal)::var_bound=0.0
#endif

  ! Refinement parameters for hydro
  real(dp)::err_grad_d=-1.0  ! Density gradient
  real(dp)::err_grad_u=-1.0  ! Velocity gradient
  real(dp)::err_grad_p=-1.0  ! Pressure gradient
  real(dp)::err_grad_A=-1.0  ! Bx gradient
  real(dp)::err_grad_B=-1.0  ! By gradient
  real(dp)::err_grad_C=-1.0  ! Bz gradient
  real(dp)::err_grad_B2=-1.0 ! B L2 norm gradient
  real(dp)::err_grad_E=-1.0  ! Radiative energy norm gradient
  real(dp)::err_grad_F=-1.0  ! Radiative flux norm gradient
  real(dp)::floor_d=1.d-10   ! Density floor
  real(dp)::floor_u=1.d-10   ! Velocity floor
  real(dp)::floor_p=1.d-10   ! Pressure floor
  real(dp)::floor_A=1.d-10   ! Bx floor
  real(dp)::floor_B=1.d-10   ! By floor
  real(dp)::floor_C=1.d-10   ! Bz floor
  real(dp)::floor_b2=1.d-10  ! B L2 norm floor
  real(dp)::floor_E=1.d-10   ! Radiative energy floor
  real(dp)::floor_F=1.d-10   ! Radiative flux floor
  real(dp)::mass_sph=0.0D0   ! mass_sph
#if NENER>0
  real(dp),dimension(1:NENER)::err_grad_prad=-1.0
#endif
#if NPSCAL>0
#if USE_M_1==0
  real(dp),dimension(1:NVAR-8-NENER)::err_grad_var=-1.0
#endif
#if USE_M_1==1
  real(dp),dimension(1:NVAR-8-NENER-nfr)::err_grad_var=-1.0
#endif
#endif
  real(dp),dimension(1:MAXLEVEL)::jeans_refine=-1.0

  ! Initial conditions hydro variables
  real(dp),dimension(1:MAXREGION)::d_region=0.
  real(dp),dimension(1:MAXREGION)::u_region=0.
  real(dp),dimension(1:MAXREGION)::v_region=0.
  real(dp),dimension(1:MAXREGION)::w_region=0.
  real(dp),dimension(1:MAXREGION)::p_region=0.
  real(dp),dimension(1:MAXREGION)::A_region=0.
  real(dp),dimension(1:MAXREGION)::B_region=0.
  real(dp),dimension(1:MAXREGION)::C_region=0.
  real(dp),dimension(1:MAXREGION)::T_region=0.
  real(dp),dimension(1:MAXREGION,1:ngrp)::E_region=0.
  real(dp),dimension(1:MAXREGION,1:ngrp)::fx_region=0.0d0
  real(dp),dimension(1:MAXREGION,1:ngrp)::fy_region=0.0d0
  real(dp),dimension(1:MAXREGION,1:ngrp)::fz_region=0.0d0
#if NENER>0
  real(dp),dimension(1:MAXREGION,1:NENER)::prad_region=0.0
#endif
#if NPSCAL>0
  real(dp),dimension(1:MAXREGION,1:npscal)::var_region=0.0
#endif

  ! Hydro solver parameters
  integer ::niter_riemann=10
  integer ::slope_type=1
  integer ::slope_mag_type=-1
  integer ::positivity_type=0
  real(dp)::slope_theta=1.5d0
  real(dp)::gamma=1.4d0
  real(dp),dimension(1:nener)::gamma_rad=1.33333333334d0
  real(dp)::courant_factor=0.5d0
  real(dp)::difmag=0.0d0
  real(dp)::smallc=1.d-10
  real(dp)::smallr=1.d-10
  real(dp)::smallcr=1.d-10
  real(dp)::eta_mag=0.0d0
  character(LEN=10)::scheme='muscl'
  character(LEN=10)::riemann='llf'
  character(LEN=10)::riemann2d='llf'
  real(dp)::switch_solv=1.d20
  real(dp)::switch_solv_dens=1.d20
  integer ::ischeme=0
  integer ::iriemann=0
  integer ::iriemann2d=0

  ! Interpolation parameters
  integer ::interpol_var=0
  integer ::interpol_type=1
  integer ::interpol_mag_type=-1

  ! EXTINCTION RELATED PARAMETERS
  ! get_dx
  real(dp)                 :: pi_g           !pi for a global calculation (done in cooling_fine)
  real(dp),dimension(1:4)  :: mod13
  real(dp),dimension(1:4)  :: mod23
  real(dp),allocatable, dimension(:,:,:,:)   :: xalpha
  !  integer, allocatable, dimension(:,:,:)     :: dirM, dirN
  integer ,allocatable, dimension(:,:)       :: Mdirection, Ndirection
  integer ,allocatable, dimension(:,:,:,:)   :: dirM_ext, dirN_ext, dirMN_ext
  real(dp),allocatable, dimension(:,:)       :: Mdx_cross_int
  real(dp),allocatable, dimension(:,:,:,:)   :: Mdx_cross_loc
  real(dp),allocatable, dimension(:,:,:,:,:) :: Mdx_ext
  logical ,allocatable, dimension(:,:,:,:,:) :: Mdx_ext_logical

  ! Passive variables index
  integer::imetal=9
  integer::idelay=9
  integer::ixion=9
  integer::ichem=9
  integer::ivirial1=9
  integer::ivirial2=9
  integer::inener=9

  real(dp),parameter::Grav=6.67e-08_dp   !Gravitational constant
  real(dp),parameter::Msun=1.9889e33_dp  !Sun mass in g
  real(dp),parameter::Rsun=6.95508e10_dp !Sun radius in cm
  real(dp),parameter::Lsun=3.846e33_dp   !Sun luminosity in erg/s
  real(dp),parameter::year=3.15576e7_dp  !1 year in s

  ! modif nimhd
  integer:: nxx=1
  integer:: nyy=2
  integer:: nzz=3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Artificial pseudo-viscosity ? Yes=1 No=0
  integer :: nvisco = 0
! coefficient of pseudo-viscosity
  real(dp):: visco=1.d0
! security factor for time-step
  real(dp):: coefvisco=0.1d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NON IDEAL MHD parameters
! Ambipolar diffusion ? Yes=1 No=0
  integer :: nambipolar = 0
  integer :: nambipolar2 = 0

! Magnetic  diffusion ? Yes=1 No=0 
  integer :: nmagdiffu  = 0                ! magnetic diffusion with multiple time stepping
  integer :: nmagdiffu2 = 0                 ! magnetic diffusion with subcycling

! Hall effect ? Yes=1 No=0 
  integer :: nhall = 0

! Magnetic  diffusion alone then magohm=1. (or flux=0.)
! predictor including eta*J*B
! magohm=0.d0
! no predictor including eta*J*B
! THIS SEEMS TO BE BETTER for oblique choc
  integer :: magohm=1

! Mellon & Li 2009 (?) or Hennebelle & Teyssier 2007
! WARNING this value is in CGS. The connection with user units
! is made in function gammaadbis in umsucl
  real(dp):: gammaAD=1.0d0
  real(dp):: rho_threshold=1.d-10     ! safeguard for the ambipolar flux in high density contrast cases (in code units)
  integer :: use_x1d=0               ! use abundances
  integer :: use_x2d=0               ! use abundances
  integer :: use_x3d=0               ! use abundance table with rho, T, Xi
  integer :: use_res=0             ! use resistivities

! magnetic diffusion coefficient see function etaohmdiss in umsucl
  real(dp):: etaMD=1.d0

! Hall resistivity
  real(dp):: rHall=1.d0

! Making a test or not Yes=1 No=0
  real(dp):: ntestDADM=0


  real(dp):: xmion=29.*1.667d-24 ! mp*mu_ion of the most abundant ion (HCO+)
  real(dp), parameter:: H2_fraction = 0.844d0 ! H2 fraction in number of particules (equals 0.73 in mass)
! WARNING !! Think to change xmolaire if proportion are changed
  real(dp):: xmneutre=(H2_fraction*2.d0 +(1.d0-H2_fraction)*4.)*1.667d-24

! Mellon & Li 2009 (?) or Hennebelle & Teyssier 2007
! WARNING this value is in CGS. The connection with user units
! is made in function densionbis in umsucl
  real(dp):: coefionis=3.d-16 ! coefionis*sqrt(n_H)=n_i , empirical value from Shu book 2, p. 363
  real(dp):: coefad = 0.1d0
  integer :: nminitimestep = 0
  real(dp):: coefalfven = 1.d-10  ! Meme si ca n'a rien a voir avec alfven : c'est le coefficient de seuil. Par defaut, on ne seuille pas.
  real(dp):: coefdtohm = 1.d-10
  real(dp):: coefohm = 0.05d0
  real(dp):: coefhall=0.05d0
  real(dp)::default_ionisrate=1d-17

! WARNING FOLLOWING VALUES IN CGS
! IF NEW ADDED MODIFY units.f90 TO CONVERT IN USER UNITS

  real(dp):: rhoi0=1.d0
  
  logical :: use_nonideal_mhd
  real(dp)::nu_sts=0.001
! fin modif nimhd


  !Cosmic rays related variables
  real(dp)::k_perp=0.01       ! Perpendicular diffusion coefficient (k_perp*kspitzer_para)
  real(dp)::R_length = 1.0d0
  real(dp)::M0 = 1.0 
  real(dp)::gamma_ei=1.0       
  real(dp)::Dcr=1.0d29 ! Classical value, in cm^2/s (e.g., Jockipii 1999)
  real(dp)::flinj=1.0d0 ! fraction of boxlen for CR injection
  ! Interpolation parameters for anisotropic diffusion
  integer ::interpol_var_cond=0
  integer ::interpol_type_cond=0
  integer ::interpol_mag_type_cond=-1
  real(dp)::frac_pcr=0.1

end module hydro_parameters
