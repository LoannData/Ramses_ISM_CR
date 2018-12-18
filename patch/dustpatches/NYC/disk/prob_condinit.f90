OB)

       omega = sqrt(omega2)

       q(i,iu) = omega * yy 

       q(i,iv) = -omega * xx

       q(i,iw) = 0.0

       q(i,ip) = q(i,id) * Cs**2


     ELSE  
        IF(rc .le. 3.5*r0) THEN 
        q(i,id) = max(d0 * r0**2 / sqrt(rc**2 + emass**2)**(3.-temper_expo/2.) * exp( - M0 / Cs**2 / 2.  * zz**2 / (rc**2+emass**2)**(1.5) ) /1000. ,1.d-5)
        q(i,iu) = omega0 * yy / (rc**2+emass**2)**(0.75)
        q(i,iv) = -omega0 * xx / (rc**2+emass**2)**(0.75)
        q(i,iw) = 0.0
        q(i,ip) = q(i,id) * Cs**2 
        ELSE
        q(i,id) = max(d0 * r0**2 / sqrt(rc**2 + emass**2)**(3.-temper_expo/2.) * exp( - M0 / Cs**2 / 2.  * zz**2 / (rc**2+emass**2)**(1.5) ) /1000. ,1.d-5)
        q(i,iu) = omega0 * yy / (rc**2+emass**2)**(0.75)
        q(i,iv) =  -omega0 * xx / (rc**2+emass**2)**(0.75)
        q(i,iw) = 0.0
        q(i,ip) = q(i,id) * Cs**2 
        ENDIF
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
