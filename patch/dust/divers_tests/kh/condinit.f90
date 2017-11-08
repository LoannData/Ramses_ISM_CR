!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use cooling_module      , only : kb,mh
  use radiation_parameters,only:mu_gas
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
  integer::ivar, idust
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp),dimension(1:nvector) :: sum_dust
  logical                       :: region_cond

  !Fo KH instability
  integer::i,j,id,iu,iv,iw,ip
  real(dp)::x0,lambday,ky,lambdaz,kz,rho1,rho2,p0,v0,v1,v2,epsilon_1,epsilon_2
  real(dp):: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

 region_cond = .true.

  IF (region_cond .eqv. .true.) then 
     ! Call built-in initial condition generator
     call region_condinit(x,q,dx,nn)
     ! Add here, if you wish, some user-defined initial conditions
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
     u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
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
     u(1:nn,nvar)=q(1:nn,5)/(gamma-1.0d0)
#endif
#if NDUST>0
     ! dust
     do ivar=1,ndust
        u(1:nn,firstindex_ndust+ivar)=q(1:nn,1)*q(1:nn,firstindex_ndust+ivar)
     end do
#endif
  else
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

     id=1; iu=2; iv=3; iw=4; ip=5
     x0=x_center(1)
     lambday=0.25
     ky=2.*acos(-1.0d0)/lambday
     lambdaz=0.25
     kz=2.*acos(-1.0d0)/lambdaz
     rho1=d_region(1)
     rho2=d_region(2)
     epsilon_1=dust_region(1,1)
     epsilon_2=dust_region(2,1)
     v1=v_region(1)
     v2=v_region(2)
     v0= 5.0
     p0= 1000

  do i=1,nn
     if(x(i,1) < x0)then
           q(i,id)=rho1*(1.0_dp+epsilon_1)
           q(i,iu)=0.0
           if(ndim>1)q(i,iu)=v0*cos(ky*(x(i,2)-lambday/2.))*exp(+ky*(x(i,1)-x0))
           if(ndim>1)q(i,iv)=v1
           if(ndim>2)q(i,iw)=0.0D0
           q(i,ip)=p0+rho1*(1.0_dp+epsilon_1)*gravity_params(1)*x(i,1)
           q(i,firstindex_ndust+1)=dust_region(1,1)

        else
          q(i,id)=rho2*(1.0_dp+epsilon_2)
          q(i,iu)=0.0
           if(ndim>1)q(i,iu)=v0*cos(ky*(x(i,2)-lambday/2.))*exp(-ky*(x(i,1)-x0))
           if(ndim>1)q(i,iv)=v2
           if (ndim>2)q(i,iw)=0.0D0
           q(i,ip)=p0+rho1*(1.0_dp+epsilon_1)*gravity_params(1)*x0+rho2*(1.0_dp+epsilon_2)*gravity_params(1)*(x(i,1)-x0)
           q(i,firstindex_ndust+1)=dust_region(2,1)

           
           
        endif
     end do

#if NDUST>0
     ! dust
        u(1:nn,firstindex_ndust+1)=q(1:nn,1)*q(1:nn,firstindex_ndust+1)
#endif
     ! Convert primitive to conservative variables
     ! density -> density
     u(1:nn,1)=q(1:nn,1)
     ! velocity -> momentum
     u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
     u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
     u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
     ! kinetic energy
     u(1:nn,5)=0.0d0
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
     ! pressure -> total fluid energy
     u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
     ! passive scalars
#if NPSCAL>0
     ! passive scalars
     do ivar=1,npscal
        u(1:nn,firstindex_pscal+ivar)=q(1:nn,1)*q(1:nn,firstindex_pscal+ivar)
     end do
#endif
 endif
  
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
