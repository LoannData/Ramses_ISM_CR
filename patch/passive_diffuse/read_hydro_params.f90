subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  use radiation_parameters
  use pm_commons
  use cooling_module,ONLY:kB,mH,clight
  use const
  use units_commons
  use mod_opacities
  use cloud_module
#if NIMHD==1
  use variables_X,ONLY:nvarchimie,nchimie,tchimie,&
      &nminchimie,tminchimie,dnchimie,dtchimie,&
      &xichimie,ximinchimie,dxichimie,&
      nislin,tislin,xiislin
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables  
  !--------------------------------------------------
  integer::i,j,idim,irad,nboundary_true=0,ht
  integer ,dimension(1:MAXBOUND)::bound_type
  real(dp)::ek_bound,em_bound,er_bound
  real(dp)::radiation_source
  character(len=2):: rad_trans_model='m1'

  integer::ii,jj,kk,ee,hh,gg,ie,ir,k,it
  real(dp)::dummy,compute_db,d0
  real(dp)::xx,yy,vv,ww,zz
  real(dp)::dtemp1,Temp_new2,epsilon_n,eint_old,T0,temp_new,d_loc,eint_new,pi

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini &
       & ,length_x,length_y,length_z,exp_region &
       & ,d_region,u_region,v_region,w_region,p_region &
#if NENER>NGRP
       & ,prad_region &
#endif
#if NGRP>0
       & ,E_region &
#endif
#if NVAR>8+NENER
       & ,var_region &
#endif
       & ,A_region,B_region,C_region &
       & ,alpha_dense_core,beta_dense_core,crit,delta_rho,mass_c,rap,cont &
       & ,ff_sct,ff_rt,ff_act,ff_vct,theta_mag,bb_test &
       & ,contrast,Mach,uniform_bmag,r0_box
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,slope_mag_type,switch_solv,switch_solv_dens &
#if NENER>0
       & ,gamma_rad &
#endif
       & ,pressure_fix,beta_fix,scheme,riemann,riemann2d &
       & ,positivity_type
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
       & ,iso_jeans,Tp_jeans &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,err_grad_A,err_grad_B,err_grad_C,err_grad_B2,err_grad_E &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,floor_A,floor_B,floor_C,floor_B2,floor_E &
       & ,interpol_var,interpol_type,sink_refine,interpol_mag_type
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max &
       & ,kbound_min,kbound_max &
       & ,d_bound,u_bound,v_bound,w_bound,p_bound &
#if NENER>NGRP
       & ,prad_bound &
#endif
#if NGRP>0
       & ,E_bound &
#endif
#if NVAR>8+NENER
       & ,var_bound &
#endif
!-----------ynlee june 2018
#if NPSCAL>0
       & ,seed_pscal,seed_high_T,i_seed,T_seed,evaporate_pscal,T_evaporate & 
!---------------------------
#endif
       & ,A_bound,B_bound,C_bound,T_bound ,no_inflow
  namelist/physics_params/cooling,haardt_madau,metal,isothermal,barotrop,eos &
       & ,m_star,t_star,n_star,T2_star,g_star,del_star,eps_star,jeans_ncells &
       & ,eta_sn,yield,rbubble,f_ek,ndebris,f_w,mass_gmc,kappa_IR &
       & ,J21,a_spec,z_ave,z_reion,eta_mag,delayed_cooling,T2max &
       & ,self_shielding,smbh,agn,B_ave,t_diss &
!       & ,rsink_max,msink_max,merge_stars &
       & ,units_density,units_time,units_length,neq_chem,ir_feedback,ir_eff &
       & ,larson_lifetime,flux_accretion,t_diss &
       & ,mu_gas,analytical_barotrop
  namelist/radiation_params/grey_rad_transfer,dtdiff_params,dt_control &
       & ,rosseland_params,planck_params,epsilon_diff,fld_limiter &
       & ,freqs_in_Hz,read_groups,split_groups_log,extra_end_group  &
       & ,numin,numax,Tr_floor,robin,rad_trans_model,min_optical_depth,rt_feedback &
       & ,PMS_evol,Hosokawa_track,energy_fix,facc_star,facc_star_lum,valp_min,store_matrix,external_radiation_field &
       & ,opacity_type,rad_trans_model,min_optical_depth &
       & ,rt_feedback,PMS_evol,Hosokawa_track,energy_fix &
       & ,facc_star,facc_star_lum,store_matrix &
       & ,external_radiation_field,stellar_photon
  ! modif nimhd
  namelist/nonidealmhd_params/nambipolar,gammaAD &
       & ,nmagdiffu,etaMD,nhall,rHall,ntestDADM &
       & ,coefad, nminitimestep, coefalfven,nmagdiffu2,nambipolar2,nu_sts,coefdtohm &
       & ,rho_threshold,use_x1d,use_x2d,use_x3d,use_res,default_ionisrate
  namelist/pseudovisco_params/nvisco,visco
  ! fin modif nimhd

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call clean_stop
102 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
#if USE_FLD==1 || USE_M_1==1
  if(FLD)read(1,NML=radiation_params)
#endif
  rewind(1)
  read(1,NML=boundary_params,END=103)
  simple_boundary=.true.
  goto 104
103 simple_boundary=.false.
104 if(nboundary>MAXBOUND)then
    write(*,*) 'Error: nboundary>MAXBOUND'
    call clean_stop
  end if
  rewind(1)
  read(1,NML=physics_params,END=105)
105 continue

  ! Conversion factor from user units to cgs units (to be done after read physics_params with units_density...)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(barotrop)fld=.false.

  !------------------------------------------------
  ! set ischeme
  !------------------------------------------------
  SELECT CASE (scheme)
  CASE ('muscl')
    ischeme = 0
  CASE ('induction')
    ischeme = 1

  CASE DEFAULT
    write(*,*)'unknown scheme'
    call clean_stop
  END SELECT
  !------------------------------------------------
  ! set iriemann
  !------------------------------------------------
  SELECT CASE (riemann)
  CASE ('llf')
    iriemann = 0
  CASE ('roe')
    iriemann = 1
  CASE ('hll')
    iriemann = 2
  CASE ('hlld')
    iriemann = 3
  CASE ('upwind')
    iriemann = 4
  CASE ('hydro')
    iriemann = 5

  CASE DEFAULT
    write(*,*)'unknown riemann solver'
    call clean_stop
  END SELECT
  !------------------------------------------------
  ! set iriemann
  !------------------------------------------------
  SELECT CASE (riemann2d)
  CASE ('llf')
    iriemann2d = 0
  CASE ('roe')
    iriemann2d = 1
  CASE ('upwind')
    iriemann2d = 2
  CASE ('hll')
    iriemann2d = 3
  CASE ('hlla')
    iriemann2d = 4
  CASE ('hlld')
    iriemann2d = 5
  CASE DEFAULT
    write(*,*)'unknown 2D riemann solver'
    call clean_stop
  END SELECT

#if NIMHD==1
  ! modif nimhd
  rewind(1)
  read(1,NML=nonidealmhd_params,END=106)
106 continue
  if((nambipolar.ne.0).and.(nambipolar.ne.1)) then
     write(*,*)'Wrong choice for nambipolar'
     call clean_stop
  end if
  if((nmagdiffu.ne.0).and.(nmagdiffu.ne.1)) then
     write(*,*)'Wrong choice for nmagdiffu'
     call clean_stop
  end if
  if((nhall.ne.0).and.(nhall.ne.1)) then
     write(*,*)'Wrong choice for nhall'
     call clean_stop
  end if

  if((nmagdiffu.eq.1).and.(nmagdiffu2.eq.1)) then
     write(*,*)'Wrong choice for nmagdiffu : choose one kind not both'
     call clean_stop
  end if

  if((nambipolar.eq.1).and.(nambipolar2.eq.1)) then
     write(*,*)'Wrong choice for nambipolar : choose one kind not both'
     call clean_stop
  end if

  if( (nambipolar.eq.1) .or. (nambipolar2.eq.1) .or. &
      (nmagdiffu .eq.1) .or. (nmagdiffu2 .eq.1) .or. &
      (nhall.eq.1) )then
     use_nonideal_mhd = .true.
  else
     use_nonideal_mhd = .false.
  endif

  if(myid==1) then
     write(*,*)'!!!!!!!!!!!!!!!  Non Ideal MHD   !!!!!!!!!!!!!!!!'
     write(*,*)'Non ideal MHD parameters'
     write(*,*)'Making a test ? Yes=1 No=0',ntestDADM
     if(nambipolar.eq.1) then
        write(*,*)'Ambipolar diffusion switched ON'
        write(*,*)'Ambipolar diffusion coefficient',gammaAD
        write(*,*)'Ambipolar diffusion time coefficient',coefad
        write(*,*)'Ionisation coefficient',coefionis
        if(nminitimestep.eq.1) then
           write(*,*)'Mini time step switched ON'
           write(*,*)'Mini time step coefficient',coefalfven
        else
           write(*,*)'Mini time step switched OFF'
        endif
     endif
     if(nambipolar2==1) then
        write(*,*)'Ambipolar diffusion switched ON : subcylcing'
        write(*,*)'Ambipolar diffusion coefficient',gammaAD
        write(*,*)'Ambipolar diffusion time coefficient',coefad
        write(*,*)'Ionisation coefficient',coefionis
        if(nminitimestep.eq.1) then
           write(*,*)'Mini time step switched ON'
           write(*,*)'Mini time step coefficient',coefalfven
        else
           write(*,*)'Mini time step switched OFF'
        endif
     endif

     if((nambipolar.eq.0) .and. (nambipolar2 == 0)) write(*,*)'Ambipolar diffusion switched OFF'

     if(nmagdiffu.eq.1)then
        write(*,*)'Magnetic diffusion switched ON : multiple time stepping'
        write(*,*)'Magnetic diffusion coefficient',etaMD
        write(*,*)'Magnetic diffusion  time coefficient',coefohm
     endif
     if(nmagdiffu2.eq.1)then
        write(*,*)'Magnetic diffusion switched ON : subcycling'
        write(*,*)'Magnetic diffusion coefficient',etaMD
        write(*,*)'Magnetic diffusion  time coefficient',coefohm
     endif
     if((nmagdiffu.eq.0).and.(nmagdiffu2.eq.0))write(*,*)'Magnetic diffusion switched OFF'

     if(nhall.eq.1)then
        write(*,*)'Hall effect switched ON'
        write(*,*)'Hall resistivity',rHall
        write(*,*)'Hall effect time coefficient',coefhall
     endif
     if(nhall.eq.0)write(*,*)'Hall effect switched OFF'

    !if(change_solver.eq.1)then             ! change solver is always used in this version
        write(*,*)'Solveur change when the time step becomes too small'
        write(*,*)'switch_solv', switch_solv
    !endif
  endif

  rewind(1)
  read(1,NML=pseudovisco_params,END=107)
107 continue
  if((nvisco.ne.0).and.(nvisco.ne.1)) then
     write(*,*)'Wrong choice for nvisco'
     call clean_stop
  end if

  if(myid==1) then
     write(*,*)'!!!!!!!!! Pseudo Viscosity Parameters  !!!!!!!!!!'
     if(nvisco.eq.1) then
        write(*,*)'Pseudo viscosity switched ON'
        write(*,*)'Pseudo viscosity coefficient',visco
        write(*,*)'Pseudo viscosity time coefficient',coefvisco
     endif

     if(nvisco.eq.0) then
        write(*,*)'Pseudo viscosity switched OFF'
     endif
  endif

  if((nambipolar2.eq.1).or.(nmagdiffu2.eq.1))then
     if(pressure_fix.eqv..false.) then
        write(*,*)'STS needs pressure_fix=.true. to work....'
        call clean_stop
     end if
  end if
#else
  rewind(1)
  read(1,NML=nonidealmhd_params,END=108)
108 continue
  if( (nambipolar.eq.1) .or. (nambipolar2.eq.1) .or. &
       (nmagdiffu .eq.1) .or. (nmagdiffu2 .eq.1) .or. &
       (nhall.eq.1) )then
     if (myid==1) write(*,*)'You must recompile with NIMHD=1 for non-ideal MHD...'
     call clean_stop
  endif
  ! fin modif nimhd
#endif

#if USE_FLD==1 || USE_M_1==1
  ! Initialize multigroup
  allocate(nu_min_hz(1:ngrp),nu_max_hz(1:ngrp),nu_min_ev(1:ngrp),nu_max_ev(1:ngrp))
  call create_groups
  call tabulate_art4
  call read_omegas
  if(myid==1 .and. grey_rad_transfer .and. ngrp .gt.1) then
     print*,'Warning: Grey Radiation Transfer with NRAD>1'
     call clean_stop
  endif
  scale_E0 = aR*(Tr_floor**4)
  P_cal = scale_E0 / (scale_d * scale_v**2)
  C_cal = clight / scale_v
  is_radiative_energy = .false.
#endif

#if USE_FLD==1
  ! Set i_fld_limiter
  i_fld_limiter=i_fld_limiter_nolim
  if(fld_limiter=='levermore') i_fld_limiter=i_fld_limiter_levermore
  if(fld_limiter=='minerbo')  i_fld_limiter=i_fld_limiter_minerbo
  ! Index array for radiative variables and temperature
  ! Needed in M1 because temperature is stored in uold(:,nvar)
  do irad = 1,nvar_bicg
     ind_bicg (irad) = firstindex_er+irad
     norm_bicg(irad) = P_cal
  enddo
  ind_trad(1) = nvar
  norm_trad(1) = Tr_floor
  do irad = 2,nvar_trad
     ind_trad(irad) = firstindex_er-1+irad
     norm_trad(irad) = P_cal
     is_radiative_energy(irad) = .true.
  enddo
#endif

#if USE_M_1==1
  ! Set radiative transfer model
  select case(rad_trans_model)
  case('P1','p1')
     irad_trans_model = irad_trans_model_p1
  case('M1','m1')
     irad_trans_model = irad_trans_model_m1
  case default
     if(myid==1) write(*,*) 'unknown radiative transfer model: '//rad_trans_model
     call clean_stop
  end select
  call compute_valp
  ! Index array for radiative variables and temperature
  ! Needed in M1 because temperature is stored in uold(:,nvar)
  ind_bicg(1) = nvar
  norm_bicg(1) = Tr_floor
  do irad = 2,nvar_bicg
     ind_bicg(irad) = firstindex_er-1+irad
     norm_bicg(irad) = P_cal
  enddo
  do irad = ngrp+2,nvar_bicg
     norm_bicg(irad) = norm_bicg(irad)*C_cal
  enddo
  ind_trad=ind_bicg
  norm_trad=norm_bicg
  is_radiative_energy(2:ngrp+1) = .true.
#endif

  ! Compute the size of the box early,
  ! to avoid problems in the initial build of the amr grid
  call calc_boxlen

  !--------------------------------------------------
  ! Make sure virtual boundaries are expanded to 
  ! account for staggered mesh representation
  !--------------------------------------------------
  nexpand_bound=2

  !--------------------------------------------------
  ! Check for star formation
  !--------------------------------------------------
  if(t_star>0)then
     star=.true.
     pic=.true.
  else if(eps_star>0)then
     t_star=0.1635449*(n_star/0.1)**(-0.5)/eps_star
     star=.true.
     pic=.true.
  endif

  !--------------------------------------------------
  ! Check for metal
  !--------------------------------------------------
  if(metal.and.nvar<(ndim+6))then
     if(myid==1)write(*,*)'Error: metals need nvar >= ndim+6'
     if(myid==1)write(*,*)'Modify hydro_parameters.f90 and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Check for non-thermal energies
  !--------------------------------------------------
#if NENER>NGRP
  if(nvar<(8+nent))then
     if(myid==1)write(*,*)'Error: non-thermal energy need nvar >= 8+nent'
     if(myid==1)write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif
  !--------------------------------------------------
  ! Check for radiative variables
  !--------------------------------------------------
#if NGRP>0
#if USE_FLD==1
  if(nvar<(8+nener))then
     if(myid==1)write(*,*)'Error: radiative energies need nvar >= 8+nent+ngrp'
#else
  if(nvar<(8+nener+nfr))then
     if(myid==1)write(*,*)'Error: radiative variables need nvar >= 8+nent+ngrp+nfr'
#endif
     if(myid==1)write(*,*)'Modify NENER, NGRP and recompile'
     nml_ok=.false.
  endif
#endif

  !-------------------------------------------------
  ! This section deals with hydro boundary conditions
  !-------------------------------------------------
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  if (simple_boundary)then

     ! Compute new coarse grid boundaries
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           nx=nx+1
           if(ibound_min(i)==-1)then
              icoarse_min=icoarse_min+1
              icoarse_max=icoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ny=ny+1
           if(jbound_min(i)==-1)then
              jcoarse_min=jcoarse_min+1
              jcoarse_max=jcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           nz=nz+1
           if(kbound_min(i)==-1)then
              kcoarse_min=kcoarse_min+1
              kcoarse_max=kcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do

     ! Compute boundary geometry
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           if(ibound_min(i)==-1)then
              ibound_min(i)=icoarse_min+ibound_min(i)
              ibound_max(i)=icoarse_min+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=1
              if(bound_type(i)==2)boundary_type(i)=11
              if(bound_type(i)==3)boundary_type(i)=21
           else
              ibound_min(i)=icoarse_max+ibound_min(i)
              ibound_max(i)=icoarse_max+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=2
              if(bound_type(i)==2)boundary_type(i)=12
              if(bound_type(i)==3)boundary_type(i)=22
           end if
           if(ndim>1)jbound_min(i)=jcoarse_min+jbound_min(i)
           if(ndim>1)jbound_max(i)=jcoarse_max+jbound_max(i)
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           if(jbound_min(i)==-1)then
              jbound_min(i)=jcoarse_min+jbound_min(i)
              jbound_max(i)=jcoarse_min+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=3
              if(bound_type(i)==2)boundary_type(i)=13
              if(bound_type(i)==3)boundary_type(i)=23
           else
              jbound_min(i)=jcoarse_max+jbound_min(i)
              jbound_max(i)=jcoarse_max+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=4
              if(bound_type(i)==2)boundary_type(i)=14
              if(bound_type(i)==3)boundary_type(i)=24
           end if
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           jbound_min(i)=jcoarse_min+jbound_min(i)
           jbound_max(i)=jcoarse_max+jbound_max(i)
           if(kbound_min(i)==-1)then
              kbound_min(i)=kcoarse_min+kbound_min(i)
              kbound_max(i)=kcoarse_min+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=5
              if(bound_type(i)==2)boundary_type(i)=15
              if(bound_type(i)==3)boundary_type(i)=25
           else
              kbound_min(i)=kcoarse_max+kbound_min(i)
              kbound_max(i)=kcoarse_max+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=6
              if(bound_type(i)==2)boundary_type(i)=16
              if(bound_type(i)==3)boundary_type(i)=26
           end if
        end if
     end do
     do i=1,nboundary
        ! Check for errors
        if( (ibound_min(i)<0.or.ibound_max(i)>(nx-1)) .and. (ndim>0) .and.bound_type(i)>0 )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along X direction',i
           nml_ok=.false.
        end if
        if( (jbound_min(i)<0.or.jbound_max(i)>(ny-1)) .and. (ndim>1) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Y direction',i
           nml_ok=.false.
        end if
        if( (kbound_min(i)<0.or.kbound_max(i)>(nz-1)) .and. (ndim>2) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Z direction',i
           nml_ok=.false.
        end if
     end do
  end if
  nboundary=nboundary_true
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  !--------------------------------------------------
  ! Compute boundary conservative variables
  !--------------------------------------------------
  do i=1,nboundary
     ! Do imposed BC for radiative transfer
     d0=compute_db()
     d_bound(i)=d0
     T_bound(i)=Tr_floor
     P_bound(i)=T_bound(i)*d_bound(i)*kb/(mu_gas*mH*scale_v**2)
     
     boundary_var(i,1)=MAX(d_bound(i),smallr)
     boundary_var(i,2)=d_bound(i)*u_bound(i)
     boundary_var(i,3)=d_bound(i)*v_bound(i)
     boundary_var(i,4)=d_bound(i)*w_bound(i)
     boundary_var(i,6)=A_bound(i)
     boundary_var(i,7)=B_bound(i)
     boundary_var(i,8)=C_bound(i)
     boundary_var(i,nvar+1)=A_bound(i)
     boundary_var(i,nvar+2)=B_bound(i)
     boundary_var(i,nvar+3)=C_bound(i)

     er_bound=0.0D0
#if NENER>0
     do j=1,nent
        boundary_var(i,firstindex_er+j)=prad_bound(i,j)
        er_bound=er_bound+boundary_var(i,8+j)/(gamma_rad(j)-1.0d0)
     end do
#endif
#if USE_FLD==1 || USE_M_1==1
     !     T_bound(i)=P_bound(i)*mu_gas*mH/kb/d_bound(i) *scale_v**2
     call temperature_eos(d_bound(i),P_bound(i)/(gamma-1.0d0),T_bound(i),ht)
     do j=1,ngrp
        boundary_var(i,firstindex_er+j)=radiation_source(T_bound(i),j)/(scale_d*scale_v**2)
        er_bound=er_bound+boundary_var(i,firstindex_er+j)
#if USE_M_1==1
        !Radiative fluxes
                   boundary_var(i,firstindex_er+  ngrp+j)=boundary_var(i,firstindex_er+j)*clight/scale_v*fx_bound(i,j)
        if(ndim>1) boundary_var(i,firstindex_er+2*ngrp+j)=boundary_var(i,firstindex_er+j)*clight/scale_v*fy_bound(i,j)
        if(ndim>2) boundary_var(i,firstindex_er+3*ngrp+j)=boundary_var(i,firstindex_er+j)*clight/scale_v*fz_bound(i,j)
#endif
     end do
#endif

     ek_bound=0.5d0*d_bound(i)*(u_bound(i)**2+v_bound(i)**2+w_bound(i)**2)
     em_bound=0.5d0*(A_bound(i)**2+B_bound(i)**2+C_bound(i)**2)
     boundary_var(i,5)=ek_bound+em_bound+er_bound+P_bound(i)/(gamma-1.0d0)
  end do

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     jeans_refine(i)=-1.0
  end do

  !-----------------------------------
  ! Sort out passive variable indices
  !-----------------------------------
  inener=9 ! MUST BE THIS VALUE !!! RT variable
  imetal=firstindex_pscal+1
  lastindex_pscal=nvar
  if(energy_fix)lastindex_pscal=nvar-1
  idelay=imetal
  if(metal)idelay=imetal+1
  ivirial1=idelay
  ivirial2=idelay
  if(delayed_cooling)then
     ivirial1=idelay+1
     ivirial2=idelay+1
  endif
  if(sf_virial)then
     if(sf_compressive) ivirial2=ivirial1+1
  endif
  ixion=ivirial2
  if(delayed_cooling)ixion=ivirial2+1
  ichem=ixion
  if(aton)ichem=ixion+1

  !-----------------------------------
  ! Set magnetic slope limiters
  !-----------------------------------
  if (slope_mag_type == -1) then
    slope_mag_type = slope_type
  endif
  if (interpol_mag_type == -1) then
    interpol_mag_type = interpol_type
  endif

#if NIMHD==1
  !------------------------------------------
  ! Read resistivity tables for non-ideal MHD
  !------------------------------------------
  if(use_nonideal_mhd)then
     if(use_res==1)then
        open(10,file='res_sig.dat', status='old')
        read(10,*) nchimie
        allocate(resistivite_chimie_res(8,nchimie))
        do i=1,nchimie
           read(10,*)resistivite_chimie_res(:,i)
        end do
        close(10)
        rho_threshold=max(rho_threshold,resistivite_chimie_res(1,1)*(mu_gas*mH)/scale_d) ! input in part/cc, output in code units
        resistivite_chimie_res(7,:)=resistivite_chimie_res(7,:)*clight**2/(4.d0*acos(-1.))
        resistivite_chimie_res(6,:)=1.43d-7**2/(&
        &max(resistivite_chimie_res(6,:)*((1.0d0-tanh(resistivite_chimie_res(1,:)/5.0d13))),1.e-36)&
        &*3.d-16*sqrt(resistivite_chimie_res(1,:))*(2.34d-24**1.5)*clight**2)
        !     open(1010,file='res_sig_v.dat', status='new')
        nminchimie=(resistivite_chimie_res(1,1))
        dnchimie=(log10(resistivite_chimie_res(1,nchimie))-log10(resistivite_chimie_res(1,1)))/&
                 &(nchimie-1)
!                 print*, dnchimie,17.d0/35.d0
        !     do i=1,nchimie
        !     write(1010,*)resistivite_chimie(1,i),resistivite_chimie(6,i)
        !  end do
        !     close(1010)
        !     stop
     else if(use_x2d==1)then
        open(42,file='resnh.dat', status='old')
        read(42,*) nchimie, tchimie, nvarchimie
        read(42,*)
        read(42,*)
        allocate(resistivite_chimie_x(-1:nvarchimie,nchimie,tchimie,1))
        do i=1,tchimie
           do j=1,nchimie
              read(42,*)resistivite_chimie_x(0:nvarchimie,j,i,1),dummy,dummy,dummy,dummy,resistivite_chimie_x(-1,j,i,1)
!              print *, resistivite_chimie_x(:,j,i)
           end do
           read(42,*)
        end do
        close(42)
        rho_threshold=max(rho_threshold,resistivite_chimie_x(0,1,1,1)*(mu_gas*mH)/scale_d) ! input in part/cc, output in code units
        nminchimie=(resistivite_chimie_x(0,1,1,1))
        dnchimie=(log10(resistivite_chimie_x(0,nchimie,1,1))-log10(resistivite_chimie_x(0,1,1,1)))/&
                 &(nchimie-1)
!                 print*, dnchimie,15.d0/50.d0
        tminchimie=(resistivite_chimie_x(-1,1,1,1))
        dtchimie=(log10(resistivite_chimie_x(-1,1,tchimie,1))-log10(resistivite_chimie_x(-1,1,1,1)))/&
                 &(tchimie-1)
!                 print*, dtchimie,3.d0/50.d0
!         close(333)
        call rq
        call nimhd_3dtable
     else if(use_x3d==1)then

        open(42,file='marchand2016_table.dat',form='unformatted')
        read(42) nchimie, tchimie, xichimie, nvarchimie
        allocate(resistivite_chimie_x(-2:nvarchimie+4,nchimie,tchimie,xichimie))
        read(42) resistivite_chimie_x
        close(42)

        rho_threshold=max(rho_threshold,resistivite_chimie_x(-2,1,1,1)*(mu_gas*mH)/scale_d) ! input in part/cc, output in code units
        nminchimie=(resistivite_chimie_x(-2,1,1,1))
        dnchimie=(log10(resistivite_chimie_x(-2,nchimie,1,1))-log10(resistivite_chimie_x(-2,1,1,1)))/&
                 &(nchimie-1)
!                 print*, dnchimie,15.d0/50.d0
        tminchimie=(resistivite_chimie_x(-1,1,1,1))
        dtchimie=(log10(resistivite_chimie_x(-1,1,tchimie,1))-log10(resistivite_chimie_x(-1,1,1,1)))/&
                 &(tchimie-1)
!                 print*, dtchimie,3.d0/50.d0
        ximinchimie=(resistivite_chimie_x(0,1,1,1))
        dxichimie=(log10(resistivite_chimie_x(0,1,1,xichimie))-log10(resistivite_chimie_x(0,1,1,1)))/&
                 &(xichimie-1)
        call rq_3d
        call nimhd_4dtable
     else
        print*, 'must choose an input for abundances or resistivities'
        stop
     endif
  endif
#endif

  if(barotrop)fld=.false.

  if(barotrop .and. (.not. analytical_barotrop))then
     open(101,file='barotropic_eos.dat', status='old')
     read(101,*)nrho_barotrop,rhomin_barotrop,rhomax_barotrop,drho_barotrop
     allocate(rho_barotrop(nrho_barotrop))
     allocate(temp_barotrop(nrho_barotrop))
     do i=1,nrho_barotrop
        read(101,*)rho_barotrop(i),temp_barotrop(i)
     end do
     close(101)
  end if

  if(eos)then
  
     !--------------------------------
     ! Read eos tables
     !--------------------------------
!      open(14,file='verif.dat')
     open(10,file='tab_eos.dat',status='old',form='unformatted')
     read(10) nRho,nEnergy
     read(10) rhomin,rhomax,emin,Emax,yHe
     
     allocate(Rho_eos(nRho,nEnergy),Ener_eos(nRho,nEnergy),Temp_eos(nRho,nEnergy),P_eos(nRho,nEnergy))
     allocate(  Cs_eos(nRho,nEnergy),S_eos(nRho,nEnergy),  xH_eos(nRho,nEnergy), xH2_eos(nRho,nEnergy)                  )
     allocate(xHe_eos(nRho,nEnergy),xHep_eos(nRho,nEnergy),Cv_eos(nRho,nEnergy)                                       )
     !inversion de la table eos
     nTemp=nEnergy
     allocate(eint_eos(nRho,nTemp))
     
     read(10)  rho_eos
     read(10) Ener_eos
     read(10) Temp_eos
     read(10)    P_eos
     read(10)    S_eos
     read(10)   Cs_eos
     read(10)   xH_eos
     read(10)  xH2_eos
     read(10)  xHe_eos
     read(10) xHep_eos
     close(10)
     
     rho_eos(:,:) = log10(rho_eos(:,:))
     ener_eos(:,:) = log10(ener_eos(:,:))
     
     do k=1,5
        ii=0
        jj=0
        kk=0
        hh=0
        ee=0
        gg=0
        do ir=2,nRho-1
           do ie=2,nEnergy-1
              if (P_eos(ir,ie) .eq. 0.0d0) then
                 ii = ii+1
                 xx = P_eos(ir,ie+1) * P_eos(ir,ie-1) *  P_eos(ir-1,ie) * P_eos(ir+1,ie)
                 yy = P_eos(ir+1,ie+1) * P_eos(ir+1,ie-1) *  P_eos(ir-1,ie-1) * P_eos(ir-1,ie+1)
                 if(ie > 2 .and. ie < nEnergy-1 .and. ir > 2 .and. ir < nRho-1)then
                    ww = P_eos(ir,ie+2) * P_eos(ir,ie-2) *  P_eos(ir-2,ie) * P_eos(ir+2,ie)
                 else
                    ww = 0.0_dp
                 endif
                 if(ie > 3 .and. ie < nEnergy-2 .and. ir > 3 .and. ir < nRho-2)then
                    zz = P_eos(ir+3,ie+3) * P_eos(ir-3,ie-3) *  P_eos(ir-3,ie+3) * P_eos(ir+3,ie-3)
                 else
                    zz = 0.0_dp
                 endif
                 if (xx .ne. 0.) then
                    P_eos(ir,ie) = 0.25d0*(P_eos(ir,ie+1) + P_eos(ir,ie-1) + P_eos(ir-1,ie) + P_eos(ir+1,ie))
                    jj=jj+1              
                 else if (yy .ne. 0. .and. k > 0) then
                    P_eos(ir,ie) = 0.25d0*(P_eos(ir+1,ie+1) + P_eos(ir+1,ie-1) + P_eos(ir-1,ie+1)+P_eos(ir-1,ie-1))
                    kk=kk+1
                 else if (ww .ne. 0 .and. k > 1) then
                    ee = ee +1
                    P_eos(ir,ie) = 0.25d0*(P_eos(ir,ie+2) + P_eos(ir,ie-2) + P_eos(ir-2,ie) + P_eos(ir+2,ie))
                 else if (zz .ne. 0 .and. k > 2) then
                    hh=hh+1
                    P_eos(ir,ie) = 0.25d0*(P_eos(ir+3,ie+3) + P_eos(ir+3,ie-3) + P_eos(ir-3,ie+3)+P_eos(ir-3,ie-3))
                 else 
                    gg=gg+1
                 endif
              endif
           enddo
        end do
        if (myid == 1) print*, "on bouche les trous P_eos", ii,jj,kk,ee,hh,gg, "iter", k
     end do
     
     do k=1,5
        ii=0
        jj=0
        kk=0
        hh=0
        ee=0
        gg=0
        do ir=2,nRho-1
           do ie=2,nEnergy-1
              if (Cs_eos(ir,ie) .eq. 0.0d0) then           
                 ii = ii+1
                 xx = Cs_eos(ir,ie+1) * Cs_eos(ir,ie-1) *  Cs_eos(ir-1,ie) * Cs_eos(ir+1,ie)
                 yy = Cs_eos(ir+1,ie+1) * Cs_eos(ir+1,ie-1) *  Cs_eos(ir-1,ie-1) * Cs_eos(ir-1,ie+1)
                 if(ie > 2 .and. ie < nEnergy-1 .and. ir > 2 .and. ir < nRho-1)then
                    ww = Cs_eos(ir,ie+2) * Cs_eos(ir,ie-2) *  Cs_eos(ir-2,ie) * Cs_eos(ir+2,ie)
                 else
                    ww = 0.0_dp
                 endif
                 if(ie > 3 .and. ie < nEnergy-2 .and. ir > 3 .and. ir < nRho-2)then
                    zz = Cs_eos(ir+3,ie+3) * Cs_eos(ir-3,ie-3) *  Cs_eos(ir-3,ie+3) * Cs_eos(ir+3,ie-3)
                 else
                    zz = 0.0_dp
                 endif
                 if (xx .ne. 0.) then
                    Cs_eos(ir,ie) = 0.25d0*(Cs_eos(ir,ie+1) + Cs_eos(ir,ie-1) + Cs_eos(ir-1,ie) + Cs_eos(ir+1,ie))
                    jj=jj+1              
                 else if (yy .ne. 0. .and. k > 0) then
                    Cs_eos(ir,ie) = 0.25d0*(Cs_eos(ir+1,ie+1) + Cs_eos(ir+1,ie-1) + Cs_eos(ir-1,ie+1)+Cs_eos(ir-1,ie-1))
                    kk=kk+1
                 else if (ww .ne. 0 .and. k > 1) then
                    ee = ee +1
                    Cs_eos(ir,ie) = 0.25d0*(Cs_eos(ir,ie+2) + Cs_eos(ir,ie-2) + Cs_eos(ir-2,ie) + Cs_eos(ir+2,ie))
                 else if (zz .ne. 0 .and. k > 2) then
                    hh=hh+1
                    Cs_eos(ir,ie) = 0.25d0*(Cs_eos(ir+3,ie+3) + Cs_eos(ir+3,ie-3) + Cs_eos(ir-3,ie+3)+Cs_eos(ir-3,ie-3))
                 else 
                    gg=gg+1
                 endif
              endif
           enddo
        end do
        if (myid == 1) print*, "on bouche les trous Cs_eos", ii,jj,kk,ee,hh,gg, "iter", k
     end do
     
     do k=1,5
        ii=0
        jj=0
        kk=0
        hh=0
        ee=0
        gg=0
        do ir=2,nRho-1
           do ie=2,nEnergy-1
              if (Temp_eos(ir,ie) .eq. 0.0d0) then           
                 ii = ii+1
                 xx = Temp_eos(ir,ie+1) * Temp_eos(ir,ie-1) *  Temp_eos(ir-1,ie) * Temp_eos(ir+1,ie)
                 yy = Temp_eos(ir+1,ie+1) * Temp_eos(ir+1,ie-1) *  Temp_eos(ir-1,ie-1) * Temp_eos(ir-1,ie+1)
                 if(ie > 2 .and. ie < nEnergy-1 .and. ir > 2 .and. ir < nRho-1)then
                    ww = Temp_eos(ir,ie+2) * Temp_eos(ir,ie-2) *  Temp_eos(ir-2,ie) * Temp_eos(ir+2,ie)
                 else
                    ww = 0.0_dp
                 endif
                 if(ie > 3 .and. ie < nEnergy-2 .and. ir > 3 .and. ir < nRho-2)then
                    zz = Temp_eos(ir+3,ie+3) * Temp_eos(ir-3,ie-3) *  Temp_eos(ir-3,ie+3) * Temp_eos(ir+3,ie-3)
                 else
                    zz = 0.0_dp
                 endif
                 if (xx .ne. 0.) then
                    Temp_eos(ir,ie) = 0.25d0*(Temp_eos(ir,ie+1)+Temp_eos(ir,ie-1)+Temp_eos(ir-1,ie)+Temp_eos(ir+1,ie))
                    jj=jj+1              
                 else if (yy .ne. 0. .and. k > 0) then
                    Temp_eos(ir,ie) = 0.25d0*(Temp_eos(ir+1,ie+1)+Temp_eos(ir+1,ie-1)+Temp_eos(ir-1,ie+1)+Temp_eos(ir-1,ie-1))
                    kk=kk+1
                 else if (ww .ne. 0 .and. k > 1) then
                    ee = ee +1
                    Temp_eos(ir,ie) = 0.25d0*(Temp_eos(ir,ie+2)+Temp_eos(ir,ie-2)+Temp_eos(ir-2,ie)+Temp_eos(ir+2,ie))
                 else if (zz .ne. 0 .and. k > 2) then
                    hh=hh+1
                    Temp_eos(ir,ie) = 0.25d0*(Temp_eos(ir+3,ie+3)+Temp_eos(ir+3,ie-3)+Temp_eos(ir-3,ie+3)+Temp_eos(ir-3,ie-3))
                 else 
                    gg=gg+1
                 endif
              endif
           enddo
        end do
        
        if (myid == 1) print*, "on bouche les trous Temp_eos", ii,jj,kk,ee,hh,gg, "iter", k
     end do
     
     Tmin=3.0d0
     Tmax=1.0d5
     dtemp1 =(log10(Tmax) - log10(Tmin))/ntemp
     eint_eos(:,:)=0.0d0
     do ir=2,nRho-1
        do it=1,ntemp
           d_loc = (10.**rho_eos(ir,1))
           T0 = 10.**(log10(Tmin) + (it-1.0d0)*dtemp1)
           
           eint_old = d_loc*kb*T0/(mu_gas*mh*(gamma-1.0d0))
           if (it >1) then
              eint_old = max(d_loc*kb*T0/(mu_gas*mh*(gamma-1.0d0)),eint_eos(ir,it-1))
           end if
           
           epsilon_n = 1.0d0
           
           do ii=1,1000
              call temperature_eos(d_loc/scale_d,eint_old/(scale_d*scale_v**2),temp_new,ht)
              if (ht==1) then
                 eint_old=0.d0
                 exit
              end if
              call temperature_eos(d_loc/scale_d,eint_old*1.001_dp/(scale_d*scale_v**2),temp_new2,ht)
              if (ht==1) then
                 eint_old=0.d0
                 exit
              end if
              
              if(abs(Temp_new2-Temp_new) .ne.0)then
                 eint_new = eint_old - (Temp_new-T0)/((Temp_new2-Temp_new)/(0.001*eint_old))
              else
                 eint_new = eint_old
              endif
              epsilon_n = abs(eint_new - eint_old)/eint_old
              eint_old = eint_new
              if  (abs(epsilon_n) .lt. 1.d-4) then
                 exit
              else if (ii==1000) then
                 print*, "newton for e(rho,T) did not converge at ", log10(d_loc),log10(T0)
              end if
           end do
           eint_eos(ir,it) = eint_old 
        enddo
     enddo
     
     Cv_eos(:,:)=0.0d0
     
     do  ir=2,nRho-1
        do  ie=2,nEnergy-1
           d_loc = (10.**rho_eos(ir,1))
           T0 = 10.**(ener_eos(ir,ie))
           
           call temperature_eos(d_loc/scale_d,(T0-0.001_dp*T0)/(scale_d*scale_v**2),temp_new,ht)
           call temperature_eos(d_loc/scale_d,(T0+0.001_dp*T0)/(scale_d*scale_v**2),temp_new2,ht)
           
           if((temp_new2-temp_new) .ne. 0.0_dp)then
              Cv_eos(ir,ie)=(0.002_dp*T0)/(temp_new2-temp_new)
           else
              Cv_eos(ir,ie) = 1.0_dp
           endif
        end do
     end do
     Cv_eos(:,nEnergy)=Cv_eos(:,nEnergy-1)
     
     if (myid == 1) print*, "ok pour le bouchage"
  end if

  ! Multigroup opacities initialization
  if(fld)then
     if((opacity_type == 'grey') .and. (ngrp > 1) .and. (.not.stellar_photon))then
        if(myid == 1)then
           write(*,*) 'WARNING! Trying to use grey opacity table with ngrp =',ngrp
           write(*,*) 'Switching to multigroup opacities'
        endif
        opacity_type = 'multigroup'
     endif
     call init_opacities
  end if
 
  if(PMS_evol .and. rt_feedback .and. Hosokawa_track)then
     open(101,file='Hosokawa_track.dat', status='old')
     read(101,*)nmdot_PMS,nm_PMS,ndata_PMS
     allocate(nb_ligne_PMS(nmdot_PMS))
     allocate(data_PMS(nmdot_PMS,nm_PMS,ndata_PMS))
     read(101,*)nb_ligne_PMS(1),nb_ligne_PMS(2),nb_ligne_PMS(3),nb_ligne_PMS(4),nb_ligne_PMS(5)
     do i=1,nmdot_PMS
        do j=1,nm_PMS
           read(101,*)data_PMS(i,j,1),data_PMS(i,j,2),data_PMS(i,j,3),data_PMS(i,j,4) ! mdot,mass,luminosity,radius
        end do
     end do
     close(101)
  end if
#if NPSCAL>0
  if(seed_pscal)then
  !-----------------------------------
  ! Read seeds of the passive scalars
  !-----------------------------------
     open(101,file='seeds_passive.dat', status='old')
     read(101,*)nseeds
     allocate(seed_type(nseeds))
     read(101,*)seed_type(:)
#if NIMHD==1
     nseeds = min(nseeds,npscal-4)
#else
     nseeds = min(nseeds,npscal-1)
#endif
     allocate(x_seed(nseeds,3))
     allocate(ax_seed(nseeds,3))
     allocate(r_seed(nseeds))
     do i=1,nseeds
       if(seed_type(i)==0)then
         read(101,*)x_seed(i,:)
       elseif (seed_type(i)==1)then
         read(101,*)x_seed(i,:),ax_seed(i,:),r_seed(i)
       end if
     end do
     close(101)
  end if

  if(seed_high_T)then
  !-----------------------------------
  ! Rearrange passive temperture recorder index
  !-----------------------------------
#if NIMHD==1
    do i=npscal-4,1,-1
#else
    do i=npscal-1,1,-1
#endif
      if(i_seed(i)>0)then
        T_seed(i_seed(i)) = T_seed(i)
        T_seed(i) = 0.
      end if
    end do
  end if
#endif
end subroutine read_hydro_params
!################################################################
!################################################################
!################################################################ 
!################################################################
!   Modification of original codes written by David H. Bailey    
!   This subroutine computes ddb(i) = dda(i)+ddb(i)
subroutine DDPDD (dda, ddb, len, itype)
use amr_commons
  implicit none
  real(dp):: e, t1, t2
  integer i, len, itype
!  complex*16:: dda(len), ddb(len)
  complex*16:: dda, ddb
!  print*,dda,ddb
  do i = 1, len
     !   Compute dda + ddb using Knuth's trick.
     t1 = real(dda) + real(ddb)
     e = t1 - real(dda)
!!$     t1 = real(dda(i)) + real(ddb(i))
!!$     e = t1 - real(dda(i))
     t2 = ((real(ddb) - e) + (real(dda) - (t1 - e)))&
     &     +imag(dda) + imag(ddb)
!!$     t2 = ((real(ddb(i)) - e) + (real(dda(i)) - (t1 - e)))&
!!$     &     +imag(dda(i)) + imag(ddb(i))
!!$ !    print*,t1,t2
     !   The result is t1 + t2, after normalization.
!!$     ddb(i) = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1) )
     ddb = cmplx ( t1 + t2, t2 - ((t1 + t2) - t1),dp)
  enddo
  
  return
end subroutine DDPDD
