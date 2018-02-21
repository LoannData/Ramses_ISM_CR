subroutine check_subcycle_dust(uin,myflux,dx,dy,dz,dt,ngrid,ncycle,dust_cycle)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust+2)::uin
  real(dp)::dx_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fx
  real(dp) :: speed, sigma,dPdx,dPdy,dPdz,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  integer::i,j,k,l,isl,idust, idens, ipress
  integer::jlo,jhi,klo,khi,ihi,ilo
  integer:: ncycle
  logical:: dust_cycle
  idens=2*ndust+1
  ipress=2*ndust+2

  !x direction
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp
        dPdx= (uin(l,i,j,k,ipress)-uin(l,i-1,j,k,ipress))/dx
        do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i-1,j,k,idust)*uin(l,i-1,j,k,ndust+idust)/uin(l,i-1,j,k,idens)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)/uin(l,i,j,k,idens)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i-1,j,k,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        !Checks the stablity condition
        do idust=1,ndust
           !First order terms
           speed  = 0.5d0*(Tksright(idust)/uin(l,i,j,k,idens)+Tksleft(idust)/uin(l,i-1,j,k,idens))*dPdx
           if (speed .ne. 0.0d0) then           
           if (dt.gt. courant_factor * dx/abs(speed)) then
               !Check for diffusion approximation validity
              !if(uin(l,i,j,k,idust)*(Tksright(idust)+Tksleft(idust))/uin(l,i,j,k,idens).gt.dt) then
              !   write (*,*) 'Diffusion instable what have you done'
              !   stop
              !else
                 dust_cycle=.true.
                 ncycle=max(ncycle,floor(dt*abs(speed)/(dx*courant_factor))+1) !+1 is to be sure to subcycle
              !endif
              endif
              endif
       end do    
    enddo
  enddo
  enddo
  enddo

  !y direction
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)

  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp
        dPdy= (uin(l,i,j,k,ipress)-uin(l,i,j-1,k,ipress))/dy
        do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i,j-1,k,idust)*uin(l,i,j-1,k,ndust+idust)/uin(l,i,j-1,k,idens)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)/uin(l,i,j,k,idens)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i,j-1,k,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        !Checks the stablity condition
        do idust=1,ndust
           !First order terms
           speed  = 0.5d0*(Tksright(idust)/uin(l,i,j,k,idens)+Tksleft(idust)/uin(l,i,j-1,k,idens))*dPdy
           if (speed .ne. 0.0d0) then
           if (dt.gt. courant_factor * dy/abs(speed)) then
              ! !Check for diffusion approximation validity
              !if(uin(l,i,j,k,idust)*(Tksright(idust)+Tksleft(idust))/uin(l,i,j,k,idens).gt.dt) then
              !   write (*,*) 'Diffusion instable what have you done'
              !!   stop
              !else
                 dust_cycle=.true.
                 ncycle=max(ncycle,floor(dt*abs(speed)/(dy*courant_factor))+1) !+1 is for the residu
              !endif
              endif
           endif   
       end do    
    enddo
  enddo
  enddo
  enddo

  !z direction
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
        Tksleft     = 0.0_dp
        Tksright    = 0.0_dp
        Tksleft_tot = 0.0_dp
        Tksright_tot= 0.0_dp
        dPdz= (uin(l,i,j,k,ipress)-uin(l,i,j,k-1,ipress))/dz
        do idust= 1, ndust
           Tksleft_tot=Tksleft_tot-uin(l,i,j,k-1,idust)*uin(l,i,j,k-1,ndust+idust)/uin(l,i,j,k-1,idens)
           Tksright_tot=Tksright_tot-uin(l,i,j,k,idust)*uin(l,i,j,k,ndust+idust)/uin(l,i,j,k,idens)
        end do
        do idust= 1, ndust
           Tksleft(idust)=  uin(l,i,j,k-1,ndust+idust)+Tksleft_tot
           Tksright(idust)= uin(l,i,j,k,ndust+idust)+Tksright_tot
        end do
        !Checks the stablity condition
        do idust=1,ndust
           !First order terms
           speed  = 0.5d0*(Tksright(idust)/uin(l,i,j,k,idens)+Tksleft(idust)/uin(l,i,j,k-1,idens))*dPdz
           if (speed .ne. 0.0d0) then           
           if (dt.gt. courant_factor * dz/abs(speed)) then
              !Check for diffusion approximation validity
              !if(uin(l,i,j,k,idust)*(Tksright(idust)+Tksleft(idust))/uin(l,i,j,k,idens).gt.dt) then
              !   write (*,*) 'Diffusion instable what have you done'
              !   stop
              !else
                 dust_cycle=.true.
                 ncycle=max(ncycle,floor(dt*abs(speed)/(dz*courant_factor))+1) !+1 is for the residu
              !endif
              endif
              endif
       end do    
    enddo
  enddo
  enddo
  enddo
end subroutine check_subcycle_dust
