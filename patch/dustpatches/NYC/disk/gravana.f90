!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,RRR,rr,rx,ry,rz,rcyl,rin,radiusin,H1,sfive,hsmooth

  ! Constant vector
  if(gravity_type==1)then
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then
     gmass=Mstar_cen! GM
     emass=2.0*dx

     emass=gravity_params(2) ! Softening length
     xmass=gravity_params(3) ! Point mass coordinates
     ymass=gravity_params(4)
     zmass=gravity_params(5)
     
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        
        rx=x(i,1)-boxlen/2.0d0
#if NDIM>1
        ry=x(i,2)-boxlen/2.0d0
#endif
#if NDIM>2
        rz=x(i,3)-boxlen/2.0d0
#endif
        rr=sqrt(rx**2+ry**2+rz**2+rsmooth**2.0)

        RRR = sqrt(rx**2.0+ry**2.0+rz**2.0)
        

        f(i,1)=-gmass*rx/rrr/rr**2.0
#if NDIM>1
        f(i,2)=-gmass*ry/rrr/rr**2.0

#endif
#if NDIM>2
        f(i,3)=-gmass*rz/rrr/rr**2.0

#endif
     end do
  end if

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=4.D0*ACOS(-1.0D0)

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2d0*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2d0*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana

