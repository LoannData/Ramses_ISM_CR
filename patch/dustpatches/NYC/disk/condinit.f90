!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use amr_commons
  use hydro_commons
  use poisson_parameters
  use units_commons
  use const
  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer :: i,j,id,iu,iv,iw,ip
  real(dp):: x0,y0,z0,rc,rs,xx,yy,zz,pi,r0,d00,B0,p0,omega0,omega
  integer :: ivar, np
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::M00,Cs,emass,Cs0
  real(dp),save::first

  id=1; iu=2; iv=3; iw=4; ip=5
  x0=0.5*boxlen!gravity_params(3)
  y0=0.5*boxlen!gravity_params(4)
  z0=0.5*boxlen!gravity_params(5)
  pi=acos(-1.0d0)



  M00=gravity_params(1)
  emass=gravity_params(2)

  if(myid .eq. 1 .and. first .eq. 0) write(*,*) 'M00,emass,x0,y0,z0',M00,emass,x0,y0,z0

  ! cloud radius equal to unity
  r0 = 5.5!boxlen / 8.
  d00 = 1e-11/scale_d
  ! cloud rotation
  omega0 = sqrt(M00)
  ! cloud sound speed
  Cs0 = sqrt(temper_iso)
  ! vertical magnetic field
!  B0 = sqrt(4.*pi/5.)/0.53*crit
   B0 = 0.


  mass_sph = d00 / 1000.  * (boxlen*(0.5**levelmin))**3

   if(myid .eq. 1 .and. first .eq. 0) write(*,*) 'd0,Cs,boxlen',d00,Cs0,boxlen
   first=1.

  DO i=1,nn
     xx=x(i,1)-x0
     yy=x(i,2)-y0
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2)


!     Cs=max(Cs0*( (rc**5+(10.*emass)**5)**(1./5.)/r0)**(-temper_expo/2.),Cs0/2.)
   
!     Cs=max(Cs0*(rc/r0)**(-temper_expo/2.),Cs0/2.)
!     if( rc .le. r0/10.) Cs = Cs0*(1/10.)**(-temper_expo/2.)
     Cs=max(Cs0*(rc/r0)**(-temper_expo/2.),Cs0/2.)

     if( rc .le. sqrt( (r0/coef_rad)**2+r0/coef_rad*abs(zz))  ) Cs = Cs0*( sqrt( (r0/coef_rad)**2+r0/coef_rad*abs(zz))/r0 )**(-temper_expo/2.)

!     if( rc .le. sqrt( (r0/10.)**2)  ) Cs = Cs0*( 1./10. )**(-temper_expo/2.)



     !Bx component
     q(i,6     ) = 0.
     q(i,nvar+1) = 0.

     !By component
     q(i,7     ) = 0.
     q(i,nvar+2) = 0.

     !Bz component
     q(i,8     ) = B0
     q(i,nvar+3) = B0

  If( temper_expo .eq. 0) then
     IF(rc .le. r0 .and. abs(zz) .le. r0/2.) THEN 
        !the exponent (3-temper_expo/2.) is there to insure that the mass profile is invariant with temper_expo
       q(i,id) = max(d00 * r0**2 / sqrt(rc**2 + emass**2)**(3.-temper_expo/2.) * exp( - M00 / Cs**2 / 2.  * zz**2 / (rc**2+emass**2)**(1.5) ) ,1.d-22/scale_d)
       q(i,iu) = omega0 * yy / (rc**2+emass**2)**(0.75)
       q(i,iv) = -omega0 * xx / (rc**2+emass**2)**(0.75)
       q(i,iw) = 0.0
       q(i,ip) = q(i,id) * Cs**2
     ELSE  
        q(i,id) = max(d00 * r0**2 / sqrt(rc**2 + emass**2)**(3.-temper_expo/2.) * exp( - M00 / Cs**2 / 2.  * zz**2 / (rc**2+emass**2)**(1.5) ) /1000. ,1.d-22/scale_d)
        q(i,iu) = omega0 * yy / (rc**2+emass**2)**(0.75)
        q(i,iv) = -omega0 * xx / (rc**2+emass**2)**(0.75)
        q(i,iw) = 0.0
        q(i,ip) = q(i,id) * Cs**2 
     ENDIF
  ELSE 
     IF(rc .le. r0 .and. abs(zz) .le. r0/2.) THEN 
        !the exponent (3-temper_expo/2.) is there to insure that the mass profile is invariant with temper_expo
       q(i,id) = max(d00 * r0**2 / sqrt(rc**2 + emass**2)**(3.-temper_expo/2.) * exp(  M00 / Cs**2 * ( 1. / ( zz**2+rc**2+emass**2)**(0.5) -  1. / ( rc**2+emass**2)**(0.5) ) ) ,1.d-22/scale_d)


       omega =  M00 / (rc**2+emass**2)**(1.5) -2.*Cs**2 * (3.-temper_expo/2.)/(rc**2+emass**2) 
!       omega = omega + temper_expo / (rc**2+emass**2) * (-1. / (zz**2+rc**2+emass**2)**(0.5) +  1. / (rc**2+emass**2)**(0.5))
!       omega = omega + temper_expo * rc**3 / (rc**5+(10.*emass)**5) * (-1. / (zz**2+rc**2+emass**2)**(0.5) +  1. / (rc**2+emass**2)**(0.5))

       if(rc .gt. r0/coef_rad) omega = omega + temper_expo / rc**2 * (-1. / (zz**2+rc**2+emass**2)**(0.5) +  1. / (rc**2+emass**2)**(0.5))
 
       omega = sqrt(omega)

       q(i,iu) = omega * yy 
       q(i,iv) = -omega * xx
       q(i,iw) = 0.0
       q(i,ip) = q(i,id) * Cs**2
     ELSE  
       q(i,id) = max(d00 * r0**2 / sqrt(rc**2 + emass**2)**(3.-temper_expo/2.) * exp(  M00 / Cs**2 * ( 1. / ( zz**2+rc**2+emass**2)**(0.5) -  1. / ( rc**2+emass**2)**(0.5) ) )/1000. ,1.d-22/scale_d)


!d00 * r0**2 / sqrt(rc**2 + emass**2)**(3.-temper_expo/2.) * exp(  M00 / Cs**2 * ( 1. / ( zz**2+rc**2+emass**2)**(0.5) -  1. / ( rc**2+emass**2)**(0.5) )


!       omega =  M00 / (rc**2+emass**2)**(1.5) -2.*Cs**2 * (3.-temper_expo/2.)/(rc**2+emass**2) 
!!       omega = omega + temper_expo / (rc**2+emass**2) * (-1. / (zz**2+rc**2+emass**2)**(0.5) +  1. / (rc**2+emass**2)**(0.5))
!       omega = omega + temper_expo * rc**3 / (rc**5+(10.*emass)**5) * (-1. / (zz**2+rc**2+emass**2)**(0.5) +  1. / (rc**2+emass**2)**(0.5))

       omega =  M00 / (rc**2+emass**2)**(1.5) -Cs**2 * (3.-temper_expo/2.)/(rc**2+emass**2)
!       omega = omega - temper_expo / (rc**2) * (-1. / (zz**2+rc**2+emass**2)**(0.5) +  1. / (rc**2+emass**2)**(0.5) + Cs**2)
!       if (rc .gt. r0/10. .and. Cs .gt. Cs0/2.) omega = omega - temper_expo / (rc**2) * (-1. / (zz**2+rc**2+emass**2)**(0.5) +  1. / (rc**2+emass**2)**(0.5) + Cs**2)
       if( rc .ge. sqrt( (r0/coef_rad)**2+r0/coef_rad*abs(zz))  ) omega = omega - temper_expo / (rc**2) * (-1. / (zz**2+rc**2+emass**2)**(0.5) +  1. / (rc**2+emass**2)**(0.5) + Cs**2)
 

       if(omega .lt. 0 .and. abs(zz) .lt. 0.5) then 
         write(*,*) 'rr,zz',rc,zz 
         write(*,*) 'first term',M00 / (rc**2+emass**2)**(1.5) -Cs**2 * (3.-temper_expo/2.)/(rc**2+emass**2)  
         write(*,*) 'omega ',omega
       endif

       if( isnan(omega) .or. isnan(q(i,id))) then 
        write(*,*) 'rr,zz',rc,zz          
       endif


       if(q(i,1) .lt. smallr * 10.) omega = 0.

       omega = max(omega,0.)
       omega = sqrt(omega)

       q(i,iu) = omega * yy 
       q(i,iv) = -omega * xx
       q(i,iw) = 0.0
       q(i,ip) = q(i,id) * Cs**2
     ENDIF
  ENDIF


  ENDDO

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
  !kinetic + magnetic energy
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,8)+q(1:nn,nvar+3))**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic field 
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  ! passive scalars
  do ivar=9,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do

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

     xx=x(i,1)
     yy=x(i,2)
     zz=x(i,3)

     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

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
     
     v(i,1)=vx
     v(i,2)=vy
     v(i,3)=vz

  end do


end subroutine velana
