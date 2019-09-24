! ------------------------------------------------------------------------------
!  DUSTDIFF_SPLIT  This routine solves the dust flux after the operator splitting  
!              
!
!  inputs/outputs
!  uin         => (const)  input state
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,    
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
!
!  uin = (\rho eps_1, .... \rho eps_ndust, tstop1, tstopndust,\rho,P)

!
!  This routine was adapted from Yohan Dubois & Benoit Commerçon's
!  heat conduction routine by Ugo Lebreuilly
! --------------------------------------------------------------------------
subroutine dustdiff_split(uin,flux,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none
  
  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::uin

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::Xflux,Yflux,Zflux

  ! Local scalar variables
  integer::i,j,k,l,ivar, idust
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  flux=0.0d0
  Xflux=0.0d0
  Yflux=0.0d0
  Zflux=0.0d0
  
 ! Compute the dust flux in X direction
  call dustXflx(uin,Xflux,dx,dt,ngrid)

  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        do idust = 1, ndust
           flux(l,i,j,k,idust,1)=Xflux(l,i,j,k,idust)
        end do
    enddo
  enddo
  enddo
  enddo
#if NDIM>1
  ! Compute dust flux in Y direction
  call dustYflx(uin,Yflux,dy,dt,ngrid)
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        do idust = 1, ndust         
           flux(l,i,j,k,idust,2)=Yflux(l,i,j,k,idust)
        end do
     enddo
  enddo
  enddo
  enddo
#endif
  
#if NDIM>2
  ! Compute the dust flux in Z direction
  call dustZflx(uin,Zflux,dz,dt,ngrid)
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
        do idust = 1, ndust         
          flux(l,i,j,k,idust,3)=Zflux(l,i,j,k,idust)
        end do
     enddo
  enddo
  enddo
  enddo
#endif




end subroutine dustdiff_split

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

subroutine dustXflx(uin,myflux,dx,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::uin

  real(dp)::dx_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fx
  real(dp) :: speed,dspeed,speedtemp,speedtempr,speedtempl, sigma,dPdxl,dPdxr,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,sgn
  real(dp):: t_dyn,entho,pi
  integer::i,j,k,l,isl,idust, idens, ipress
  integer::jlo,jhi,klo,khi
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  entho= one /(gamma -one)
  pi =3.14159265358979323846_dp

  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        do idust=1,ndust
           !First order terms
           speed =  speedx
           fx(idust)= max(speed*uin(l,i-1,j,k,idust),0.0d0)+min(speed*uin(l,i,j,k,idust),0.0d0)
           !Second order terms
           if(speed.ge.0.0d0) isl = i-1
           if(speed.lt.0.0d0) isl = i
           if(slope_dust.eq.1)sigma=0.0d0
           if(slope_dust.eq.2)call minmod_dust((uin(l,isl,j,k,idust)-uin(l,isl-1,j,k,idust))/dx,(uin(l,isl+1,j,k,idust)-uin(l,isl,j,k,idust))/dx,sigma)
           if(slope_dust.eq.3)call vanleer(0.5d0*(uin(l,isl+1,j,k,idust)-uin(l,isl-1,j,k,idust))/dx,2.0d0*(uin(l,isl,j,k,idust)-uin(l,isl-1,j,k,idust))/dx,&
                 &2.0d0*(uin(l,isl+1,j,k,idust)-uin(l,isl,j,k,idust))/dx,sigma)             
           fx(idust) = fx(idust) + 0.5d0*abs(speed)*(dx-abs(speed)*dt)*sigma
        end do
        do idust= 1, ndust
           myflux(l,i,j,k,idust)=fx(idust)*dt/dx
        end do
    enddo
  enddo
  enddo
  enddo
  
end subroutine dustXflx

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

subroutine dustYflx(uin,myflux,dy,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dy,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::uin

  real(dp)::dPdyl,dPdyr,dy_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fy
  !Slopes and advection velocity
  real(dp) :: speed,dspeed,speedtempr,speedtempl, sigma,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp):: t_dyn,entho,pi
  ! Local scalar variables
  integer::i,j,k,l,ivar, idust, idens, ipress,isl
  integer::ilo,ihi,klo,khi
  entho= one /(gamma -one)
  pi =3.14159265358979323846_dp

  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        do idust=1,ndust
           speed =  speedy
           fy(idust)= max(speed*uin(l,i,j-1,k,idust),0.0d0)+ min(speed*uin(l,i,j,k,idust),0.0d0)
           !Second order terms
           if(speed.ge.0.0d0) isl = j-1
           if(speed.lt.0.0d0) isl = j
           if(slope_dust.eq.1)sigma=0.0d0
           if(slope_dust.eq.2)call minmod_dust((uin(l,i,isl,k,idust)-uin(l,i,isl-1,k,idust))/dy,(uin(l,i,isl+1,k,idust)-uin(l,i,isl,k,idust))/dy,sigma)
           if(slope_dust.eq.3)call vanleer(0.5d0*(uin(l,i,isl+1,k,idust)-uin(l,i,isl-1,k,idust))/dy,2.0d0*(uin(l,i,isl,k,idust)-uin(l,i,isl-1,k,idust))/dy,&
                 &2.0d0*(uin(l,i,isl+1,k,idust)-uin(l,i,isl,k,idust))/dy,sigma)           
           fy(idust) = fy(idust) + 0.5d0*abs(speed)*(dy-abs(speed)*dt)*sigma
        end do
        do idust= 1, ndust
           myflux(l,i,j,k,idust)=fy(idust)*dt/dy
        end do  
    enddo
  enddo
  enddo
  enddo

end subroutine dustYflx

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

subroutine dustZflx(uin,myflux,dz,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dz,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::uin
  real(dp)::dPdzl,dPdzr,dz_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fz

  !Slopes and advection velocity
  real(dp) :: speed,speedtempr,speedtempl,dspeed, sigma ,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp):: t_dyn ,entho, pi
  ! Local scalar variables
  integer::i,j,k,l,ivar, idust, idens, ipress, isl
  integer::ilo,ihi,jlo,jhi,klo,khi
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  entho= one /(gamma -one)
  pi =3.14159265358979323846_dp

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid

        do idust=1,ndust
           speed =  speedz
           fz(idust)= max(speed*uin(l,i,j,k-1,idust),0.0d0)+ min(speed*uin(l,i,j,k,idust),0.0d0)
           !Second order terms
           if(speed.ge.0.0d0) isl = k-1
           if(speed.lt.0.0d0) isl = k
           if(slope_dust.eq.1)sigma=0.0d0
           if(slope_dust.eq.2)call minmod_dust((uin(l,i,j,isl,idust)-uin(l,i,j,isl-1,idust))/dz,(uin(l,i,j,isl+1,idust)-uin(l,i,j,isl,idust))/dz,sigma)
           if(slope_dust.eq.3)call vanleer(0.5d0*(uin(l,i,j,isl+1,idust)-uin(l,i,j,isl-1,idust))/dz,2.0d0*(uin(l,i,j,isl,idust)-uin(l,i,j,isl-1,idust))/dz,&
                 &2.0d0*(uin(l,i,j,isl+1,idust)-uin(l,i,j,isl,idust))/dz,sigma)
           fz(idust) = fz(idust) + 0.5d0*abs(speed)*(dz-abs(speed)*dt)*sigma
        end do
        do idust= 1, ndust
           myflux(l,i,j,k,idust)=fz(idust)*dt/dz
        end do  
    enddo
  enddo
  enddo
  enddo

end subroutine dustZflx


subroutine minmod_dust(a,b,sigma)
  use amr_parameters

  implicit none
  real(dp)::a,b
  real(dp):: sigma
  if (abs(a).gt.abs(b)) sigma=b
  if (abs(b).ge.abs(a)) sigma=a
  if (a*b.le.0.0d0) sigma=0.0d0
end subroutine minmod_dust

subroutine vanleer(a,b,c,sigma)
  use amr_parameters

  implicit none
  real(dp)::a,b,c
  real(dp)::sigma
  real(dp)::sigma2
  call minmod_dust(a,b,sigma2)
  call minmod_dust(sigma2,c,sigma)

end subroutine vanleer
