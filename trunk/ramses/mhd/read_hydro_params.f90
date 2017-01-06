subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  use radiation_parameters
  use pm_commons
  use cooling_module,ONLY:kB,mH,clight
  use const
  use hydro_parameters
  use units_commons
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
       & ,A_region,B_region,C_region
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,slope_mag_type,switch_solv &
#if NENER>0
       & ,gamma_rad &
#endif
       & ,pressure_fix,beta_fix,scheme,riemann,riemann2d
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,err_grad_A,err_grad_B,err_grad_C,err_grad_B2,err_grad_E &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine &
       & ,floor_A,floor_B,floor_C,floor_B2,floor_E &
       & ,interpol_var,interpol_type,sink_refine,interpol_mag_type
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max &
       & ,kbound_min,kbound_max &
#if NENER>0
       & ,prad_bound &
#endif
#if NVAR>8+NENER
       & ,var_bound &
#endif
       & ,d_bound,u_bound,v_bound,w_bound,p_bound,no_inflow &
#if NENER>NGRP
       & ,prad_bound &
#endif
#if NGRP>0
       & ,E_bound &
#endif
       & ,A_bound,B_bound,C_bound
  namelist/physics_params/cooling,haardt_madau,metal,isothermal,barotrop,eos &
       & ,m_star,t_star,n_star,T2_star,g_star,del_star,eps_star,jeans_ncells &
       & ,eta_sn,yield,rbubble,f_ek,ndebris,f_w,mass_gmc,kappa_IR &
       & ,J21,a_spec,z_ave,z_reion,eta_mag,delayed_cooling,T2max &
       & ,self_shielding,smbh,agn,B_ave,t_diss &
!       & ,rsink_max,msink_max,merge_stars &
       & ,units_density,units_time,units_length,neq_chem,ir_feedback,ir_eff &
       & ,larson_lifetime,flux_accretion,t_diss &
       & ,mu_gas
  namelist/radiation_params/grey_rad_transfer,dtdiff_params,dt_control &
       & ,rosseland_params,planck_params,epsilon_diff,fld_limiter &
       & ,freqs_in_Hz,read_groups,split_groups_log,extra_end_group  &
       & ,numin,numax,Tr_floor,robin,rad_trans_model,min_optical_depth,rt_feedback &
       & ,PMS_evol,Hosokawa_track,energy_fix,facc_star,facc_star_lum,valp_min,store_matrix,external_radiation_field
  ! modif nimhd
  namelist/nonidealmhd_params/nambipolar,gammaAD &
       & ,nmagdiffu,etaMD,nhall,rHall,ntestDADM &
       & ,coefad, nminitimestep, coefalfven,nmagdiffu2,nambipolar2,nu_sts,coefdtohm
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
  ixion=idelay
  if(delayed_cooling)ixion=idelay+1
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
