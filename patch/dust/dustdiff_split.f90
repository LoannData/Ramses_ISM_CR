! ---------------------------------------------------------------
!  DUSTDIFF_SPLIT  This routine solves the dust flux  
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

  real(dp)::dx_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fx
  real(dp) :: speed, sigma
  integer::i,j,k,l,isl,idust, idens, ipress
  integer::jlo,jhi,klo,khi

  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  idens=2*ndust+1
  ipress=2*ndust+2
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        dx_loc=max(ffdx(l,i-1,j,k),ffdx(l,i,j,k))
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp
        do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i-1,j,k,idust)*uin(l,i-1,j,k,ndust+idust)/uin(l,i-1,j,k,idens)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)/uin(l,i,j,k,idens)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i-1,j,k,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        
        do idust=1,ndust
           !First order terms
           speed  = 0.5d0*(Tksright(idust)/uin(l,i,j,k,idens)+Tksleft(idust)/uin(l,i-1,j,k,idens))*(uin(l,i,j,k,ipress)-uin(l,i-1,j,k,ipress))/dx
           if(speed.ge.0.0d0) fx(idust)= speed*uin(l,i-1,j,k,idust) 
           if(speed<0.0d0) fx(idust)= speed*uin(l,i,j,k,idust)
           !Second order terms
           if(speed.ge.0.0d0) isl = i-1
           if(speed<0.0d0) isl = i       
           call minmod_dust((uin(l,isl,j,k,idust)-uin(l,isl-1,j,k,idust))/dx,(uin(l,isl+1,j,k,idust)-uin(l,isl,j,k,idust))/dx,sigma)
           fx(idust) = fx(idust) + 0.5d0*abs(speed)*(dx-abs(speed)*dt)*sigma
        end do
        do idust= 1, ndust
           myflux(l,i,j,k,idust)=fx(idust)*dt/dx/dx_loc
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
  !Slopes and advection velocity
  real(dp) :: speed, sigma
  ! Local scalar variables
  integer::i,j,k,l,ivar, idust, idens, ipress,isl
  integer::ilo,ihi,klo,khi
  
  idens=2*ndust+1
  ipress=2*ndust+2
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)

  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        dy_loc=max(ffdy(l,i,j-1,k),ffdy(l,i,j,k))
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp
        dPdy1=(uin(l,i,j,k,ipress)-uin(l,i,j-1,k,ipress))/dy
         do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i,j-1,k,idust)*uin(l,i,j-1,k,ndust+idust)/uin(l,i,j-1,k,idens)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)/uin(l,i,j,k,idens)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i,j-1,k,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        do idust=1,ndust
           !First order terms
           speed  = 0.5d0*(Tksright(idust)/uin(l,i,j,k,idens)+Tksleft(idust)/uin(l,i,j-1,k,idens))*dPdy1
           if(speed.ge.0.0d0) fy(idust)= speed*uin(l,i,j-1,k,idust) 
           if(speed<0.0d0) fy(idust)= speed*uin(l,i,j,k,idust)
           !Second order terms
           if(speed.ge.0.0d0) isl = j-1
           if(speed<0.0d0) isl = j

           call minmod_dust((uin(l,i,isl,k,idust)-uin(l,i,isl-1,k,idust))/dy,(uin(l,i,isl+1,k,idust)-uin(l,i,isl,k,idust))/dy,sigma)
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

  !Slopes and advection velocity
  real(dp) :: speed, sigma 
  ! Local scalar variables
  integer::i,j,k,l,ivar, idust, idens, ipress, isl
  integer::ilo,ihi,jlo,jhi,klo,khi

  idens=2*ndust+1
  ipress=2*ndust+2
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
 
     do l = 1, ngrid
        dz_loc=max(ffdz(l,i,j,k),ffdz(l,i,j,k-1))
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp
        dPdz1=(uin(l,i,j,k,ipress)-uin(l,i,j,k-1,ipress))/dz
         do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i,j,k-1,idust)*uin(l,i,j,k-1,ndust+idust)/uin(l,i,j,k-1,idens)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)/uin(l,i,j,k,idens)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i,j,k-1,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        do idust=1,ndust
           !First order terms
           speed  = 0.5d0*(Tksright(idust)/uin(l,i,j,k,idens)+Tksleft(idust)/uin(l,i,j,k-1,idens))*dPdz1
           if(speed.ge.0.0d0) fz(idust)= speed*uin(l,i,j,k-1,idust) 
           if(speed<0.0d0) fz(idust)= speed*uin(l,i,j,k,idust)
           !Second order terms
           if(speed.ge.0.0d0) isl = k-1
           if(speed<0.0d0) isl = k
           
           call minmod_dust((uin(l,i,j,isl,idust)-uin(l,i,j,isl-1,idust))/dz,(uin(l,i,j,isl+1,idust)-uin(l,i,j,isl,idust))/dz,sigma)
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
  use hydro_parameters, only : upwind_dust

  implicit none
  real(dp)::a,b
  real(dp):: sigma
  if (abs(a).gt.abs(b)) sigma=b
  if (abs(b).ge.abs(a)) sigma=a
  if (a*b.le.0.0d0) sigma=0.0d0
  if (upwind_dust) sigma=0.0d0
end subroutine minmod_dust

subroutine vanleer(a,b,c,sigma)
  use amr_parameters
  use hydro_parameters, only : upwind_dust

  implicit none
  real(dp)::a,b,c
  real(dp)::sigma
  real(dp)::sigma2
  call minmod_dust(a,b,sigma2)
  call minmod_dust(sigma2,c,sigma)
  if (upwind_dust) sigma=0.0d0

end subroutine vanleer
