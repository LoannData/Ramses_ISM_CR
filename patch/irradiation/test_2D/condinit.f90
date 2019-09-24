!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_commons,only:myid
  use amr_parameters
  use hydro_parameters
  use radiation_parameters
  use units_commons
  use rt_parameters, ONLY :rt_src_x_center,rt_src_y_center
#if NDIM>2
  use rt_parameters, ONLY :rt_src_z_center
#endif
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:3): d.u,d.v, U(i,4): E, U(i,5:7): Bleft, 
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:3):u,v, Q(i,4): P, Q(i,5:7): Bleft, 
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 7, remaining variables (8:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp):: x0,y0,rc,rs,xx,yy,rd,radiation_source
  integer :: i
  real(dp)::Temp,rho0,h,rho,ee,pi
  real(dp)::rosseland_ana,planck_ana
  real(dp)::au ! innner disc radius
  real(dp)::hr,zd
#if NDIM>2
  real(dp)::z0,zz
#endif
  
  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........
  boxlen = 1.0d0
  x0=rt_src_x_center(1)
  y0=rt_src_y_center(1)
#if NDIM>2
  z0=rt_src_z_center(1)
#endif
  pi=acos(-1.0d0)
  rd  = Rin
  zd=boxlen/4.
  ! for flarer disk (next step)
!  rd = 0.1d0
!  zd = 0.01d0

  ! AU in code units
  au = 1.496e13/scale_l

  rho0=rho_disk0/scale_d !2.874d-18/scale_d!8.321d-18/scale_d
  
  Temp=Tr_floor

  q=0.0d0
  
  dO i=1,nn
     xx=x(i,1)-x0
     yy=x(i,2)-y0
     
     q(i,2) = 0.0d0
     q(i,3) = 0.0d0
     q(i,4) = 0.0d0
  
#if NDIM>1
     !it is not a cylindrical radius
     rs  = sqrt(xx**2+yy**2)
!!$     rho = rho0*((rs/rd)**(-2.0d0))
     !rs  = sqrt(xx**2) !flarer disk, closer to pinte's test
     hr=zd*abs(xx)
     !hr = zd * (abs(xx)/rd)**1.25
     rho = rho0*((rs/rd)**(-2.0d0))*exp(-pi/4.*(yy/hr)**2.)
     !rho = rho0*((rs/rd)**(-2.625d0))*exp(-0.5d0*(yy/hr)**2.)                                                                                           
#endif

#if NDIM>2
     !Same disk as in the 2D tests (ie not much flared)
     zz  = x(i,3)-z0
     hr  = zd*sqrt(xx**2+yy**2)
     !it is not cylindrical radius
     rs  = sqrt(xx**2+yy**2+zz**2)
     rho = rho0*((rs/rd)**(-2.0d0))*exp(-pi/4.*(zz/hr)**2.)
#endif
     if (rs .lt. Rin) then
        rho = smallr
     endif

     !rho=max(rho,smallr)
     rho=max(rho,1.0d0)

     !smoothing for truncation ?
     !   if (rs .lt. 1.0d0) then
     !      rho = 6.25*10**10*rho0*((rs/rd)**(3.0d0))
     !   endif

     temp=tr_floor

     q(i,1)    = rho
     call enerint_eos(q(i,1),Temp,ee)
     q(i,5   ) = ee
     q(i,nvar) = ee

     do ivar=1,ngrp
        if (ivar==1 .and. stellar_photon) then
           q(i,firstindex_er+ivar) = eray_min/(scale_d*scale_v**2)
        else
           q(i,firstindex_er+ivar) = max(radiation_source(Temp,ivar)/(scale_d*scale_v**2),eray_min/(scale_d*scale_v**2))
        endif
     enddo

     ! print rosseland ana for position in the disk midplane and above the star, valid for p/k tests => to adapt in 2D
     !A
!     if (rs .lt. 0.124 .and. rs .gt. 0.122 .and. abs(zz) .lt. 0.0062) then
!        print *, "A : rs zz kr1(T) kr2(T)= ", rs, zz, rosseland_ana(rho*scale_d,Temp,Temp,1), rosseland_ana(rho*scale_d,Temp,Temp,2)
     !B
!     else if (zz .lt. 0.124d0 .and. zz .gt. 0.122d0 .and. rs .lt. sqrt(2.0d0)*0.025d0) then
!        print *, "B : rs zz kr1(T) kr2(T)= ", rs, zz, rosseland_ana(rho*scale_d,Temp,Temp,1), rosseland_ana(rho*scale_d,Temp,Temp,2) 
     !Abis                                                                     
!     else if (rs .lt. 1.25d0 .and. rs .gt. 1.24d0 .and. abs(zz) .lt. 0.06) then
!        print *, "A': rs zz kr1(T) kr2(T)= ", rs, zz, rosseland_ana(rho*scale_d,Temp,Temp,1), rosseland_ana(rho*scale_d,Temp,Temp,2)
     !Bbis         
!     else if (zz .lt. 1.27d0 .and. zz .gt. 1.26d0 .and. rs .lt. sqrt(2.0d0)*0.1d0) then                                      
!        print *, "B': rs zz kr1(T) kr2(T)= ", rs, zz, rosseland_ana(rho*scale_d,Temp,Temp,1), rosseland_ana(rho*scale_d,Temp,Temp,2) 
!     endif


  enddo

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=0.0d0!q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)
  ! magnetic energy -> total fluid energy
!  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
!  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
!  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
!  u(1:nn,6:8)=q(1:nn,6:8)
!  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
#if NENER>0
  ! non-thermal pressure -> non-thermal energy
  ! non-thermal energy   -> total fluid energy
  do ivar=1,nener-ngrp
     u(1:nn,8+ivar)=q(1:nn,8+ivar)/(gamma_rad(ivar)-1.0d0)
     u(1:nn,5)=u(1:nn,5)+u(1:nn,8+ivar)
  enddo
 ! Radiative transfer
#if NGRP>0
  ! radiative energy   -> total fluid energy
  do ivar=1,ngrp
     u(1:nn,firstindex_er+ivar)= q(1:nn,firstindex_er+ivar)
     u(1:nn,5)=u(1:nn,5)+ u(1:nn,firstindex_er+ivar)
  enddo
#if USE_M_1==1
  ! radiative flux
  do ivar=1,ndim*ngrp
     do i=1,ncache
        u(1:nn,fisrtindex_fr+ivar)=q(1:nn,firstindex+ivar)
     end do
     write(ilun)xdp
  end do
#endif
#endif
#endif
#if NEXTINCT>0
  ! Extinction
  if(extinction)u(1:nn,firstindex_extinct+nextinct)=1.0D0
#endif
#if NPSCAL>0
  ! passive scalars
  do ivar=1,npscal
     u(1:nn,firstindex_pscal+ivar)=q(1:nn,1)*q(1:nn,firstindex_pscal+ivar)
  end do
  ! Internal energy
  u(1:nn,nvar)=q(1:nn,5)
#endif

end subroutine condinit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters  
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

!!$     xx=x(i,1)
!!$#if NDIM > 1
!!$     yy=x(i,2)
!!$#endif
!!$#if NDIM > 2
!!$     zz=x(i,3)
!!$#endif
!!$     ! ABC
!!$     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
!!$     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
!!$     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega
     
!!$     v(i,1)=vx
!!$#if NDIM > 1
!!$     v(i,2)=vy
!!$#endif
!!$#if NDIM > 2
!!$     v(i,3)=vz
!!$#endif

  end do


end subroutine velana
!================================================================
!================================================================
!================================================================
!================================================================
subroutine calc_boxlen

  implicit none

  return

end subroutine calc_boxlen



!************************************************************************
SUBROUTINE add_radiation_sources(ilevel,dt)

! Inject radiation from RT source regions (from the RT namelist). Since
! the light sources are continuously emitting radiation, this is called
! continuously during code execution, rather than just during
! initialization.
!
! ilevel => amr level at which to inject the radiation
! dt     => timestep for injection (since injected values are per time)
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters, ONLY :rt_nsource,rt_source_type,rt_src_x_center,rt_src_y_center &
       &,rt_src_length_x,rt_src_length_y,rt_src_start,rt_src_end,rt_exp_source,rt_src_group,rt_n_source 
  use radiation_parameters
  use hydro_parameters, ONLY : firstindex_er
!  use rt_hydro_commons
  use hydro_commons, ONLY: nvar,uold
  implicit none
  integer::ilevel
  real(dp)::dt
  integer::i,igrid,ncache,iskip,ngrid
  integer::ind,idim,ivar,ix,iy,iz,nx_loc
  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:ngrp),save::uu
!------------------------------------------------------------------------
  call add_UV_background(ilevel)
  if(numbtot(1,ilevel)==0)return    ! no grids at this level
  if(rt_nsource .le. 0) return      ! no rt sources
  if(verbose)write(*,111)ilevel

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  ncache=active(ilevel)%ngrid
  ! dx (and dx_loc=dx) are just equal to 1/nx (where 1 is the boxlength)
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     ! Loop over cells
     do ind=1,twotondim
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))!*scale
           end do
        end do
        ! Read the RT variables
        do ivar=1,ngrp
           do i=1,ngrid
              uu(i,ivar)=uold(ind_cell(i),firstindex_er+ivar)
              uold(ind_cell(i),5) = uold(ind_cell(i),5) - uold(ind_cell(i),firstindex_er+ivar)
           end do
        end do
        ! find injected values per cell
        call radiation_sources_vsweep(xx,uu,dx,dt,ngrid)
        ! Write the RT variables
        do ivar=1,ngrp
           do i=1,ngrid
              uold(ind_cell(i),firstindex_er+ivar)=uu(i,ivar)
              uold(ind_cell(i),5) = uold(ind_cell(i),5) + uold(ind_cell(i),firstindex_er+ivar)
           end do
        end do
     end do
     ! End loop over cells
  end do
  ! End loop over grids

111 format('   Entering add_rt_sources for level ',I2)

END SUBROUTINE add_radiation_sources



!************************************************************************
SUBROUTINE radiation_sources_vsweep(x,uu,dx,dt,nn)

! Do a vector sweep, injecting RT source regions into cells, that is if
! they are in any of these regions.
!
! x      =>  ncells*ndim: positions of grid cells
! uu    <=  ncells*nrtvars: injected rt variables in each cell
! dx     =>  real cell width in code units
! dt     =>  real timestep length in code units
! nn     =>  int number of cells
!------------------------------------------------------------------------
  use amr_commons
  use rt_parameters, ONLY :rt_nsource,rt_source_type,rt_src_x_center,rt_src_y_center,rt_src_z_center &
       &,rt_src_length_x,rt_src_length_y,rt_src_length_z,rt_src_start,rt_src_end,rt_exp_source,rt_src_group,rt_n_source
  use radiation_parameters
  use hydro_commons, ONLY: nvar,uold
  use cooling_module,only: clight
  implicit none
  integer::ilevel
  integer ::nn
  real(dp)::dx,dt!,dx_cgs,dt_cgs
  real(dp),dimension(1:nvector,1:ngrp)::uu
  real(dp),dimension(1:nvector,1:ndim)::x
  integer::i,k,group_ind
  real(dp)::vol,r,xn,yn,zn,en,pi
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(dp)::scale_np,scale_fp
  real(dp)::radiation_source,lum_group,lum_star,rstar_adim
  integer:: igrp
  integer,dimension(1:nvector)::ind_cell

  pi=acos(-1.0d0)

!------------------------------------------------------------------------
 ! Initialize everything to zero
  !  uu=0.0d0
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call rt_units(scale_np, scale_fp)
  !dx_cgs=dx*scale_l
  !dt_cgs=dt*scale_t
  rstar_adim = rstar*6.96d10 / scale_l
#if NDIM>1
  lum_star=(5.67d-5*(Tstar**4))*4.0d0*3.1415d0*(rstar_adim)**2/(scale_d*(scale_v)**3)/(2.d0*Rin)
#endif
#if NDIM>2
  lum_star=(5.67d-5*(Tstar**4))*4.0d0*3.1415d0*(rstar_adim)**2/(scale_d*(scale_v)**3)
#endif
!  write(*,*) "lumstar=",lum_star, Tstar,rstar,rstar_adim
  ! Loop over RT regions
  do k=1,rt_nsource

     !if ((t-rt_src_start(k)) .lt. 0.) cycle
     !if(((t-rt_src_end(k)) .gt. 0.) .and. (rt_src_end(k) .gt. 0.)) cycle
     ! For "square" regions only:
     if(rt_source_type(k) .eq. 'square')then
       ! Exponent of choosen norm
        en=rt_exp_source(k)
        do i=1,nn
           ! Compute position in normalized coordinates
           xn=0.0d0; yn=0.0d0; zn=0.0d0
           xn=2.0d0*abs(x(i,1)-rt_src_x_center(k))/rt_src_length_x(k)
#if NDIM>1
           yn=2.0d0*abs(x(i,2)-rt_src_y_center(k))/rt_src_length_y(k)
#endif
#if NDIM>2
           zn=2.0d0*abs(x(i,3)-rt_src_z_center(k))/rt_src_length_z(k)
#endif
           ! Compute cell "radius" relative to region center
           if(rt_exp_source(k)<10)then
              r=(xn**en+yn**en+zn**en)**(1.0/en)
           else
              r=max(xn,yn,zn)
           end if
           ! If cell lies within region, inject value
           if(r<1.0)then
!              uu(i,group_ind) = rt_n_source(k)/rt_c_cgs/scale_Np
              ! The input flux is the fraction Fp/(c*Np) (Max 1 magnitude)
           end if
        end do
     end if

     ! For "point" regions only:
     if(rt_source_type(k) .eq. 'point')then
        ! Volume elements
        vol=dx**ndim!dx_cgs**ndim !to convert into code units
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=1.0; yn=1.0; zn=1.0
           xn=max(1.0-abs(x(i,1)-rt_src_x_center(k))/dx, 0.0_dp)
#if NDIM>1
           yn=max(1.0-abs(x(i,2)-rt_src_y_center(k))/dx, 0.0_dp)
#endif
#if NDIM>2
           zn=max(1.0-abs(x(i,3)-rt_src_z_center(k))/dx, 0.0_dp)
#endif
           r=xn*yn*zn
           if(r .gt. 0.) then
              ! If cell lies within CIC cloud, inject value.
              !weight = 1.0d0 !careful with energy injection normalization
              if(stellar_photon)then
                 igrp=1    ! Put all stellar radiative flux in the first group                      
                 uu(i,igrp) = uu(i,igrp) + lum_star*r/vol*dt
                 !uold(ind_cell(i),5     )=uold(ind_cell(i),5     ) + lum_star*r/vol*dtnew(ilevel)!/((4.0d0*pi*rmax**3)/3.0d0) !if weight=1.   
                 !uold(ind_cell(i),8+igrp)=uold(ind_cell(i),8+igrp) + lum_star*r/vol*dtnew(ilevel)!/((4.0d0*pi*rmax**3)/3.0d0) !if weight=1.
              else
                 do igrp=1,ngrp
#if NDIM>1
                    lum_group = radiation_source(Tstar,igrp)/(scale_d*scale_v**2)*(pi*rstar_adim**2*clight/scale_v)/(2.d0*Rin)
#endif
#if NDIM>2
                    lum_group = radiation_source(Tstar,igrp)/(scale_d*scale_v**2)*(pi*rstar_adim**2*clight/scale_v)
#endif
                    !/((4.0d0*pi*rmax**3)/3.0d0) !if weight=1, in 3D
                    uu(i,igrp) = uu(i,igrp) + lum_group*dt*r/vol ! for normalization
                    !uold(ind_cell(i),5)=uold(ind_cell(i),5) + lum_group*r/vol*dtnew(ilevel) !way of adding Erad with radiation patch
                    !uold(ind_cell(i),8+igrp)=uold(ind_cell(i),8+igrp) + lum_group*r/vol*dtnew(ilevel)
                 end do
              endif
           endif
        end do
     end if

     ! For shell regions only:
     if(rt_source_type(k) .eq. 'shell')then
        ! An emitting spherical shell with center coordinates given,
        ! along with inner and outer radius (rt_src_length_x,z).
        ! Compute CIC weights relative to region center
        do i=1,nn
           xn=0.0; yn=0.0; zn=0.0
           xn=max(abs(x(i,1)-rt_src_x_center(k)), 0.0_dp)
#if NDIM>1
           yn=max(abs(x(i,2)-rt_src_y_center(k)), 0.0_dp)
#endif
#if NDIM>2
           zn=max(abs(x(i,3)-rt_src_z_center(k)), 0.0_dp)
#endif
           r=sqrt(xn**2+yn**2+zn**2)
           if(r .gt. rt_src_length_x(k) .and. &
                r .lt. rt_src_length_y(k)) then
              ! If cell lies within CIC cloud, inject value
              ! photon input is in # per sec...need to convert to uu
!              uu(i,group_ind)=rt_n_source(k) / scale_np
           endif
        end do
     end if
  end do

  return
END SUBROUTINE radiation_sources_vsweep

