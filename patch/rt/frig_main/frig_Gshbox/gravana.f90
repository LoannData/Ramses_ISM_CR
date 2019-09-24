!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,ff,dx,ncell)
  use amr_parameters
!  use poisson_parameters  
!  use amr_commons
  use poisson_commons
  use hydro_commons
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::ff ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i
  real(dp)::gmass,emass,xmass,ymass,zmass,rr,rx,ry,rz

  real(dp):: fourpi,G,sigma,a1,a2,z0,Myear_s,mp,ff_max


  ! Constant vector
  if(gravity_type==1)then 
     do idim=1,ndim
        do i=1,ncell
           ff(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2)then 
     gmass=gravity_params(1) ! GM
     emass=gravity_params(2) ! Softening length
     xmass=gravity_params(3) ! Point mass coordinates
     ymass=gravity_params(4)
     zmass=gravity_params(5)
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-xmass
#if NDIM>1
        ry=x(i,2)-ymass
#endif
#if NDIM>2
        rz=x(i,3)-zmass
#endif
        rr=sqrt(rx**2+ry**2+rz**2+emass**2)
        ff(i,1)=-gmass*rx/rr**3
#if NDIM>1
        ff(i,2)=-gmass*ry/rr**3
#endif
#if NDIM>2
        ff(i,3)=-gmass*rz/rr**3
#endif
     end do
  end if




  a1=1.42d-3
  a2=5.49d-4
  z0=0.18*1.d3

  mp=1.660531d-24
  Myear_s=1.d6*365.*3600.*24.

  a1=1.d3*a1/(Myear_s**2)/mp/6.67d-8
  a2=a2/(Myear_s**2)/mp/6.67d-8

  fourpi=4.D0*ACOS(-1.0D0)
  G=1.0d0 !Le code resout Poisson avec G=1
!  sigma=mass_tot/(boxlen**2) 

!  sigma = multipole(1)/(boxlen**2) 
!  ff_max=(a1*0.5*boxlen)/(((0.5*boxlen)**2+z0**2)**0.5) + a2*(0.5*boxlen) 
!  sigma = sigma - 2.*ff_max / fourpi

!  do i=1,ncell
!     x(i,3)=x(i,3)-0.5*boxlen
!!     if (x(i,3)>0d0) f(i,3)=f(i,3)-fourpi*G*0.5*(sigma)
!!     if (x(i,3)<0d0) f(i,3)=f(i,3)+fourpi*G*0.5*(sigma)
!     ff(i,3)=-(a1*x(i,3))/(((x(i,3))**2+z0**2)**0.5) + a2*(x(i,3)) 
!     if (x(i,3)>0d0) ff(i,3)=ff(i,3)-fourpi*G*0.5*(sigma)
!     if (x(i,3)<0d0) ff(i,3)=ff(i,3)+fourpi*G*0.5*(sigma)
!  end do



  sigma = multipole(1)/(boxlen**2) 
  ff_max=(a1*0.5*boxlen)/(((0.5*boxlen)**2+z0**2)**0.5) + a2*(0.5*boxlen) 
!  sigma = sigma - 2.*ff_max / fourpi

  do i=1,ncell
     x(i,3)=x(i,3)-0.5*boxlen
!     if (x(i,3)>0d0) f(i,3)=f(i,3)-fourpi*G*0.5*(sigma)
!     if (x(i,3)<0d0) f(i,3)=f(i,3)+fourpi*G*0.5*(sigma)
     ff(i,3)=-(a1*x(i,3))/(((x(i,3))**2+z0**2)**0.5) + a2*(x(i,3)) 

     ff(i,3) = ff(i,3) * sigma / (2.*ff_max / fourpi)

!     if (x(i,3)>0d0) ff(i,3)=ff(i,3)-fourpi*G*0.5*(sigma)
!     if (x(i,3)<0d0) ff(i,3)=ff(i,3)+fourpi*G*0.5*(sigma)
  end do



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



subroutine phi_ana_3(xx,pp,ngrid,z1,mass_all)
  use amr_commons
  use poisson_commons
  use hydro_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::pp
  real(dp),dimension(1:nvector,1:ndim)::xx,ff
  real(dp)::z1
  real(kind=8)::mass_all
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------
  integer :: i
  real(dp):: fourpi,G,sigma,a1,a2,z0,Myear_s,mp,ff_max

  a1=1.42d-3
  a2=5.49d-4
  z0=0.18*1.d3

  mp=1.660531d-24
  Myear_s=1.d6*365.*3600.*24.

  a1=1.d3*a1/(Myear_s**2)/mp/6.67d-8
  a2=a2/(Myear_s**2)/mp/6.67d-8

  fourpi=4.D0*ACOS(-1.0D0)
  G=1.0d0 !Le code resout Poisson avec G=1

  !sigma=multipole(1)/(boxlen**2) 

  !ff_max=(a1*0.5*boxlen)/(((0.5*boxlen)**2+z0**2)**0.5) + a2*(0.5*boxlen) 

  do i=1,ngrid
     xx(i,3)=xx(i,3)-0.5*boxlen
     pp(i)=a1*((xx(i,3)**2+z0**2)**0.5-z0)+a2*0.5*xx(i,3)**2
     
     pp(i) = pp(i) + ((xx(i,3)**2+z1**2)**0.5-z1)*fourpi*mass_all/(boxlen**3)*((boxlen/2)**2+z1**2)**0.5
     !pp(i) = pp(i) * sigma / (2.*ff_max / fourpi)
  end do

end subroutine phi_ana_3

