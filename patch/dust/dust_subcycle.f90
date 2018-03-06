!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine dust_cycle_fine(ilevel,d_cycle_ok,ncycle)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the dust diffusion scheme.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated.
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,idust
  integer,dimension(1:nvector),save::ind_grid
  logical:: d_cycle_ok
  integer :: icycle, ncycle,ncycle_all,info
  if(numbtot(1,ilevel)==0)return
  d_cycle_ok=.false.
  ncycle =0

   !Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call dustcycle1(ind_grid,ngrid,ilevel,d_cycle_ok,ncycle)
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ncycle,ncycle_all,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,info)
  ncycle=ncycle_all
#endif
  if(ncycle.gt.1) d_cycle_ok =.true.
  if(.not.d_cycle_ok) ncycle =1  

  if (myid==1) write(*,112) ilevel, ncycle
112 format('   Subcycling level ',i2, ' for dust with ncycle = ',i2)

end subroutine dust_cycle_fine
  !###########################################################
  !###########################################################
  !###########################################################
  !###########################################################    
subroutine dustcycle1(ind_grid,ncache,ilevel,d_cycle_ok,ncycle)
  use amr_commons
  use hydro_commons
  use radiation_parameters, only:mu_gas
  use poisson_commons
  use cooling_module
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  ! This routine check if subcycle is required in a patch

  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3)::u2
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust+2)::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::facdx
  integer ,dimension(1:nvector,1:threetondim     )::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       )::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         )::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         )::ind1
  integer ,dimension(1:nvector                   )::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer ::  ncycle
  logical :: d_cycle_ok

  integer::idust,ht
  integer::i,j,ivar,idim,irad,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::sum_dust,sum_dust_new,sum_dust_old
  real(dp)::d,u,v,w,A,B,C,enint,e_kin,e_mag,pressure,cs, temp
  real(dp)::rho_gas, pi, t_stop
  real(dp), dimension(1:ndust) ::d_grain,l_grain

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  !Saved variables set to 0
  u1   = 0.0d0
  u2   = 0.0d0
  uloc = 0.0d0
  facdx= 0.0d0
  oneontwotondim = 1.d0/dble(twotondim)
  
  ! Initialisation of dust related quantities
  
  pi =3.14159265358979323846_dp
  if(mrn.eqv..true.) then
     call size_dust(l_grain)
     do idust=1,ndust
       l_grain(idust) = l_grain(idust)/scale_l
       d_grain(idust)=grain_dens(idust)/scale_d
    end do
  else
     do idust=1,ndust
       d_grain(idust)=grain_dens(idust)/scale_d
       l_grain(idust)=grain_size(idust)/scale_l
    end do
 endif
  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
     
  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)
  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     !Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
           nbuffer=nbuffer+1
           ind_nexist(nbuffer)=i
           ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do
     !If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3 
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro(u1,ind1,u2,nbuffer)
     end if
     !Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        !Gather dust variables
        do idust=1,ndust
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,idust)=uold(ind_cell(i),firstindex_ndust+idust)
              facdx(ind_exist(i),i3,j3,k3)=1.0d0
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,idust)=u2(i,ind_son,firstindex_ndust+idust)
           end do
        end do
        do i=1,nbuffer
           if(interpol_type.gt.0)then
              facdx(ind_nexist(i),i3,j3,k3)=1.0d0
           else
              facdx(ind_nexist(i),i3,j3,k3)=1.5d0
           endif
        end do
        !Computes and stores (in uloc) the density(2ndust+1),
        !the pressure(2ndust+2), and the stopping time (between ndust and 2ndust)
        do i=1,nexist
           if (energy_fix) then
              d=uold(ind_cell(i),1)
              enint=uold(ind_cell(i),nvar)
           else
              d=uold(ind_cell(i),1)

              u=uold(ind_cell(i),2)/d
              v=uold(ind_cell(i),3)/d
              w=uold(ind_cell(i),4)/d

              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
              B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
              C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+uold(ind_cell(i),8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=uold(ind_cell(i),5)
              else   
                 enint=uold(ind_cell(i),5)-e_kin- e_mag
              endif
           endif

           sum_dust=0.0_dp
           do idust = 1, ndust
              sum_dust = sum_dust + uloc(ind_exist(i),i3,j3,k3,ndust)/d
           end do

           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pressure)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           if(dust_barr) pressure = (1.0_dp-sum_dust)*d*cs*cs

           !We fill uloc with the quantities required to check for subcycling (d, P, epsilon, ts/d)
           
           uloc(ind_exist(i),i3,j3,k3,2*ndust+1)=d
           uloc(ind_exist(i),i3,j3,k3,2*ndust+2)=pressure
           do idust = 1, ndust
              ! different prescriptions for t-stop
              t_stop = d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d!/d
              if(K_drag) t_stop = sum_dust*(1.0_dp-sum_dust)*d/K_dust(idust)
              if(dust_barr) t_stop = 0.1_dp
              uloc(ind_exist(i),i3,j3,k3,ndust+idust)= t_stop / (1.0_dp - sum_dust)
           end do   
        end do

        do i=1,nbuffer
           if (energy_fix) then
              d=u2(i,ind_son,1)
              enint=u2(i,ind_son,nvar)
           else
              d=u2(i,ind_son,1)
              u=u2(i,ind_son,2)/d
              v=u2(i,ind_son,3)/d
              w=u2(i,ind_son,4)/d
              e_mag=0.0_dp
#ifdef SOLVERmhd             
              A=0.5d0*(u2(i,ind_son,6)+u2(i,ind_son,nvar+1))
              B=0.5d0*(u2(i,ind_son,7)+u2(i,ind_son,nvar+2))
              C=0.5d0*(u2(i,ind_son,8)+u2(i,ind_son,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)

#if NENER>0
              do irad=1,nener
                 e_kin=e_kin+u2(i,ind_son,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=u2(i,ind_son,5)
              else   
                 enint=u2(i,ind_son,5)-e_kin-e_mag
              endif
           endif
 
           uloc(ind_nexist(i),i3,j3,k3,2*ndust+1)=d
           sum_dust=0.0_dp
           do idust = 1, ndust
              sum_dust = sum_dust + u2(i,ind_son,firstindex_ndust+idust)/u2(i,ind_son,1)
           end do

           call pressure_eos((1.0_dp-sum_dust)*d,enint,pressure)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint,cs) 

           if(dust_barr) cs = 1.0_dp
           if(dust_barr) pressure = (1.0_dp-sum_dust)*d*cs*cs

           uloc(ind_nexist(i),i3,j3,k3,2*ndust+1)=d
           uloc(ind_nexist(i),i3,j3,k3,2*ndust+2)=pressure
           do idust = 1, ndust
              ! different prescriptions for t-stop
              t_stop =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d
              if(K_drag)  t_stop = sum_dust*(1.0_dp-sum_dust)*d/K_dust(idust)
              if(dust_barr) t_stop = 0.1_dp              
              uloc(ind_nexist(i),i3,j3,k3,ndust+idust) = t_stop /(1.0_dp -sum_dust)
           enddo
        end do

     end do
     end do
     end do
     !End loop over cells
  end do
  end do
  end do
  !End loop over neighboring grids
  
  !-----------------------------------------------
  !Check if subcycling
  !-----------------------------------------------
  call check_subcycle_dust(uloc,dx,dx,dx,dtnew(ilevel),ncache,ncycle,d_cycle_ok)

end subroutine dustcycle1

subroutine check_subcycle_dust(uin,dx,dy,dz,dt,ngrid,ncycle,dust_cycle)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Output fluxes

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
#if NDIM>1
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
#endif
#if NDIM>2
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
#endif
end subroutine check_subcycle_dust
