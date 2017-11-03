! ---------------------------------------------------------------
!  DUSTDIFF_SPLIT  This routine solves the dust flux following the  
!              anistropic diffusion.
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
! ----------------------------------------------------------------
subroutine dustdiff_split(uin,flux,dx,dy,dz,dt,ngrid,fdx)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none
  
  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust+2)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::fdx

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
 ! Compute the dust flux in X direction
  call dustXflx(uin,Xflux,dx,dt,ngrid,fdx)
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
  call dustYflx(uin,Yflux,dy,dt,ngrid,fdx)
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
  call dustZflx(uin,Zflux,dz,dt,ngrid,fdx)
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

subroutine dustXflx(uin,myflux,dx,dt,ngrid,ffdx)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust+2)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdx

  real(dp)::dPdx1,dx_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fx
  integer::i,j,k,l,ivar, idust
  integer::jlo,jhi,klo,khi

  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        dx_loc=max(ffdx(l,i-1,j,k),ffdx(l,i,j,k))
        dPdx1=(uin(l,i,j,k,2*ndust+2)-uin(l,i-1,j,k,2*ndust+2))/dx
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp        
        do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i-1,j,k,idust)*uin(l,i-1,j,k,ndust+idust)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i-1,j,k,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        do idust= 1, ndust
           fx(idust) = 0.5d0*(uin(l,i-1,j,k,idust)*Tksleft(idust)+uin(l,i,j,k,idust)*Tksright(idust))*dPdx1/dx_loc
           if(dust_lin) fx(idust) = -D_lin_dust*(uin(l,i,j,k,idust)-uin(l,i-1,j,k,idust))/dx/dx_loc

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

subroutine dustYflx(uin,myflux,dy,dt,ngrid,ffdy)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dy,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust+2)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdy

  real(dp)::dPdy1,dy_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fy
  
  ! Local scalar variables
  integer::i,j,k,l,ivar, idust
  integer::ilo,ihi,klo,khi

  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)


  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        dy_loc=max(ffdy(l,i,j-1,k),ffdy(l,i,j,k))
        dPdy1=(uin(l,i,j,k,2*ndust+2)-uin(l,i,j-1,k,2*ndust+2))/dy
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp
        do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i,j-1,k,idust)*uin(l,i,j-1,k,ndust+idust)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i,j-1,k,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        do idust= 1, ndust
           fy(idust) =0.5d0*(uin(l,i,j-1,k,idust)*Tksleft(idust)+uin(l,i,j,k,idust)*Tksright(idust))*dPdy1/dy_loc
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

subroutine dustZflx(uin,myflux,dz,dt,ngrid,ffdz)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dz,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust+2)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ffdz

  real(dp)::dPdz1,dz_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fz

  ! Local scalar variables
  integer::i,j,k,l,ivar, idust
  integer::ilo,ihi,jlo,jhi,klo,khi


  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
 
     do l = 1, ngrid
        dPdz1=(uin(l,i,j,k,2*ndust+2)-uin(l,i,j,k-1,2*ndust+2))/dz
        dz_loc=max(ffdz(l,i,j,k),ffdz(l,i,j,k-1))
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp
        do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i,j,k-1,idust)*uin(l,i,j,k-1,ndust+idust)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i,j,k-1,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        do idust= 1, ndust
           fz(idust) =0.5d0*(uin(l,i,j,k-1,idust)*Tksleft(idust)+uin(l,i,j,k,idust)*Tksright(idust))*dPdz1/dz_loc
        end do        
        do idust= 1, ndust
           myflux(l,i,j,k,idust)=fz(idust)*dt/dz
        end do  
    enddo
  enddo
  enddo
  enddo

end subroutine dustZflx
