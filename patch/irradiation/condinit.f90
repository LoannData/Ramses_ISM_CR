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
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft, 
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft, 
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables

  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,rd,zd,radiation_source
  integer :: i
  real(dp)::Temp,rho0,h,rho,ee,pi
  real(dp)::rosseland_ana,planck_ana
  real(dp)::rin,au ! innner disc radius

  ! Call built-in initial condition generator
  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........

  x0=0.5*boxlen
  y0=0.5*boxlen
  z0=0.5*boxlen
  pi=acos(-1.0d0)

  ! AU in code units
  au = 1.496e13/scale_l

  if(irradiation_test == 'Kuiper') then
     ! Kuiper test
     !boxlen = 2000 for kuiper (namelist)
     rd  = 500.* au ! boxlen/4=500AU
     zd  = 125.* au ! boxlen/16=125AU
     rin = 1. * au

!!$     Teff_starsink1=5800.
!!$     r_starsink1=1.
  else if(irradiation_test == 'Pinte') then
     ! Pinte test
     !boxlen = 800 (namelist)
     rd  = 100.* au !boxlen/8=100AU
     zd  = 10. * au !boxlen/80=10AU
     rin = 0.1 * au

!!$     Teff_starsink1=4000.
!!$     r_starsink1=2.
!!$  else
!!$     if(myid == 1) write(*,*) 'Bad value for irradiation_test:',irradiation_test
!!$     call clean_stop
  endif

!!$  if(myid == 1) then
!!$     write(*,*) 'Warning: Teff_starsink1 and r_starsink1 have been overwritten:',Teff_starsink1,r_starsink1
!!$  endif

  rho0=rho_disk0/scale_d !2.874d-18/scale_d!8.321d-18/scale_d
  
  Temp=Tr_floor

  q=0.0d0
  
  dO i=1,nn
     xx=x(i,1)-x0
     yy=x(i,2)-y0
     zz=x(i,3)-z0
     
     q(i,2) = 0.0d0
     q(i,3) = 0.0d0
     q(i,4) = 0.0d0
  
     !cylindrical radius !
     rs=sqrt(xx**2+yy**2)!sqrt(xx**2+yy**2+zz**2)

     h=zd*(rs/rd)**1.125

     if(irradiation_test == 'Kuiper') then
        rho=rho0*((rs/rd)**(-1.0d0))*exp(-pi/4.0d0*(abs(zz)/h)**2)
     else if(irradiation_test == 'Pinte') then
        rho=rho0*((rs/rd)**(-2.625))*exp(-0.5d0*(abs(zz)/h)**2)
     endif

     ! truncated disc below Rin
     if (rs .lt. Rin) then
        rho = smallr
     endif

     !rho=max(rho,smallr)
     rho=max(rho,1.0d0)

     !smoothing for truncation ?
     !   if (rs .lt. 1.0d0) then
     !      rho = 6.25*10**10*rho0*((rs/rd)**(3.0d0))*exp(-pi/4.0d0*(abs(zz)/h)**2)
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

     ! print rosseland ana for position in the disk midplane and above the star                                                                                                                                                                 
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
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
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
