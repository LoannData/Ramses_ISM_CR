!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_f_dust_new(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets the dust fluxes to 0 before calling
  ! the diffusion scheme. 
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,irad,idust
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Set fdust = 0 for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do idust=1,ndust
        do i=1,active(ilevel)%ngrid
           f_dust(active(ilevel)%igrid(i)+iskip,idust) = 0.0d0
        end do
     end do
  end do

  ! Set fdust to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do idust=1,ndust
        do i=1,reception(icpu,ilevel)%ngrid
          f_dust(reception(icpu,ilevel)%igrid(i)+iskip,idust)= 0.0d0
        end do
     end do
  end do
  end do

111 format('   Entering set_f_dust for level ',i2)

end subroutine set_f_dust_new
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew_dust(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew =uold for the quantities that have been updated because of dust diffusion.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,irad,idust

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           unew(active(ilevel)%igrid(i)+iskip,1) = uold(active(ilevel)%igrid(i)+iskip,1)
           unew(active(ilevel)%igrid(i)+iskip,2) = uold(active(ilevel)%igrid(i)+iskip,2)
           unew(active(ilevel)%igrid(i)+iskip,3) = uold(active(ilevel)%igrid(i)+iskip,3)
           unew(active(ilevel)%igrid(i)+iskip,4) = uold(active(ilevel)%igrid(i)+iskip,4)
           unew(active(ilevel)%igrid(i)+iskip,5) = uold(active(ilevel)%igrid(i)+iskip,5)
           do idust=1,ndust
              unew(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust) = uold(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust)
           end do
        end do
  end do

111 format('   Entering set_unew_dust for level ',i2)

end subroutine set_unew_dust
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold_dust(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold =unew for the quantities that have been updated because of dust diffusion.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,irad,idust

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
           uold(active(ilevel)%igrid(i)+iskip,1) = unew(active(ilevel)%igrid(i)+iskip,1)
           uold(active(ilevel)%igrid(i)+iskip,2) = unew(active(ilevel)%igrid(i)+iskip,2)
           uold(active(ilevel)%igrid(i)+iskip,3) = unew(active(ilevel)%igrid(i)+iskip,3)
           uold(active(ilevel)%igrid(i)+iskip,4) = unew(active(ilevel)%igrid(i)+iskip,4)
           uold(active(ilevel)%igrid(i)+iskip,5) = unew(active(ilevel)%igrid(i)+iskip,5)
           do idust=1,ndust
              uold(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust) = unew(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust)
           end do
        end do
  end do

111 format('   Entering set_uold_dust for level ',i2)

end subroutine set_uold_dust
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine dust_diffusion_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the dust diffusion scheme.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,idust
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel


  !Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call dustdifffine1(ind_grid,ngrid,ilevel)
  end do

  do idust=1,ndust
      call make_virtual_reverse_dp(f_dust(1,idust),ilevel)
  end do

111 format('   Entering diffusion_fine for level ',i2)

end subroutine dust_diffusion_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine dustdifffine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use radiation_parameters, only:mu_gas
  use poisson_commons
  use cooling_module
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !---------------------------------------------------------------------!
  ! This routine gathers first hydro variables from neighboring grids  -!
  ! to set initial conditions in a 6x6x6 grid. It interpolate from     -!
  ! coarser level missing grid variables. It then calls the            -!
  ! dust diffusion solver that compute the flux. This flux is zeroed at-!
  ! coarse-fine boundaries, since contribution from finer levels has   -!
  ! already been taken into account. Conservative variables are updated-!
  ! and stored in array unew(:), both at the current level and at the  -!
  ! coarser level if necessary.                                        -!
  !---------------------------------------------------------------------!

  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust+2),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust),save::uuloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:5),save::uuuloc

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::facdx
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust,1:ndim),save::flux

  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok

  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  integer ,dimension(1:nvector                   ),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::idust
  integer::i,j,ivar,idim,irad,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  
  real(dp)::dx,scale,oneontwotondim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::sum_dust,sum_dust_new,sum_dust_old
  real(dp)::d,u,v,w,A,B,C,enint,e_kin,e_mag,pressure,cs
  real(dp)::rho_gas,rho_grain_loc,size_grain_loc, pi
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  if(verbosed) write(*,*) 'dt_dust =', dtnew(ilevel)
  if(verbosed) write(*,1500) scale_l,scale_t,scale_d,scale_v
  
  1500 FORMAT(1pE16.8,1pE16.8,1pE16.8)
  !Saved variables set to 0
  u1   = 0.0d0
  u2   = 0.0d0
  flux = 0.0d0  
  uloc = 0.0d0
  uuloc= 0.0d0
  uuuloc= 0.0d0
  facdx= 0.0d0
  ok   = .false.
  oneontwotondim = 1.d0/dble(twotondim)
  ! Initialisation of dust related quantities
  pi =3.14159265358979323846_dp
  rho_grain_loc=rho_grain/scale_d
  size_grain_loc=size_grain/scale_l
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
           do ivar=1,ndust
              do i=1,nbuffer
                 u1(i,j,firstindex_ndust+idust)=uold(ibuffer_father(i,j),firstindex_ndust+idust)/uold(ibuffer_father(i,j),1)
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
        !Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do
        !Gather dust variables
        do idust=1,ndust
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,idust)=uold(ind_cell(i),firstindex_ndust+idust)/uold(ind_cell(i),1)
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
#ifdef solverMHD                           
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
              sum_dust = sum_dust + uold(ind_cell(i),firstindex_ndust+idust)/uold(ind_cell(i),1)
           end do
           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pressure)
           if(dust_barr) pressure = (1.0_dp-sum_dust)*d*kb*temp_dust/mu_gas/mh/scale_v**2
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = SQRT (gamma*kb*temp_dust/mu_gas/mh)/scale_v
            

           uuuloc(ind_exist(i),i3,j3,k3,2)=u
           uuuloc(ind_exist(i),i3,j3,k3,3)=v
           uuuloc(ind_exist(i),i3,j3,k3,4)=w
           uuuloc(ind_exist(i),i3,j3,k3,5)=enint+e_mag
           !We fill uloc with the quantities required to compute the fluxes (d, P, epsilon, ts/d)
           uloc(ind_exist(i),i3,j3,k3,2*ndust+1)=d
           uloc(ind_exist(i),i3,j3,k3,2*ndust+2)=pressure
           do idust = 1, ndust
              uloc(ind_exist(i),i3,j3,k3,ndust+idust)= rho_grain_loc*size_grain_loc*SQRT(pi*gamma/8.0_dp)/cs/d!/d
              if(sum_dust*uloc(ind_exist(i),i3,j3,k3,ndust+idust).gt. dtnew(ilevel))then
                 write (*,*) 'DUST DIFFUSION UNSTABLE WHAT HAVE YOU DONE?'
                 stop
              endif   
           !(price&laibe 2015)
           end do   
        if(verbosed)write(*,*)'dt, ts',sum_dust,d,enint,e_mag,e_kin,uold(ind_cell(i),5)
        if(verbosed)write(*,*) dtnew(ilevel), rho_grain_loc*size_grain_loc/d*SQRT(pi*gamma/8.0_dp)/cs
        if(verbosed)write(*,*)'cs_barren = ', cs*scale_v, ' cm.s-1'
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
#ifdef solverMHD             
              A=0.5d0*(u2(i,ind_son,6)+u2(i,ind_son,nvar+1))
              B=0.5d0*(u2(i,ind_son,7)+u2(i,ind_son,nvar+2))
              C=0.5d0*(u2(i,ind_son,8)+u2(i,ind_son,nvar+3))
              e_kin=0.5d0*d*(u**2+v**2+w**2)
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
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
              sum_dust = sum_dust + u2(i,ind_son,firstindex_ndust+idust)!/u2(i,ind_son,1)
           end do
           call pressure_eos((1.0_dp-sum_dust)*d,enint,pressure)
           if(dust_barr) pressure = (1.0_dp-sum_dust)*d*kb*temp_dust/mu_gas/mh/scale_v**2
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint,cs)
           if(dust_barr) cs = SQRT (gamma*kb*temp_dust/mu_gas/mh)/scale_v
           uloc(ind_nexist(i),i3,j3,k3,2*ndust+1)=d
           uloc(ind_nexist(i),i3,j3,k3,2*ndust+2)=pressure
           do idust = 1, ndust           
              uloc(ind_nexist(i),i3,j3,k3,ndust+idust) = rho_grain_loc*size_grain_loc*SQRT(pi*gamma/8.0_dp)/cs/d!/d
              if(sum_dust*uloc(ind_nexist(i),i3,j3,k3,ndust+idust).gt. dtnew(ilevel))then
                 write (*,*) 'DUST DIFFUSION UNSTABLE WHAT HAVE YOU DONE?'
                 stop
              endif    
              !(price&laibe 2015)
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
  ! Compute flux due to dust diffusion
  !-----------------------------------------------
  call dustdiff_split(uloc,flux,dx,dx,dx,dtnew(ilevel),ncache,facdx)
  !Reset fluxes at refined interfaces
  do idim=1,ndim
      i0=0; j0=0; k0=0
      if(idim==1)i0=1
      if(idim==2)j0=1
      if(idim==3)k0=1
      do k3=k3min,k3max+k0
      do j3=j3min,j3max+j0
      do i3=i3min,i3max+i0
         do idust=1,ndust
            do i=1,ncache
               if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                  flux(i,i3,j3,k3,idust,idim)=0.0d0
               end if
            end do
         end do
      end do
      end do
      end do
   end do
  !--------------------------------------------------------
  !Conservative  update at level ilevel for the dust fluxes
  !--------------------------------------------------------
  do idim = 1, ndim
  i0=0; j0=0; k0=0
  if(idim==1)i0=1
  if(idim==2)j0=1
  if(idim==3)k0=1
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     ind_son=1+i2+2*j2+4*k2
     iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        !Update epsilon
        do idust = 1, ndust
           do i=1,ncache
           if(son(ind_cell(i))==0)then        
                   uuloc(i,i3,j3,k3,idust)= uuloc(i,i3,j3,k3,idust)&
                   &+(flux(i,i3,j3,k3,idust,idim)&
                   &-flux(i,i3+i0,j3+j0,k3+k0,idust,idim))!/uold(ind_cell(i),1)
                   !print*,uuloc(i,i3,j3,k3,idust),flux(i,i3,j3,k3,idust,idim),flux(i,i3+i0,j3+j0,k3+k0,idust,idim),f_dust(ind_cell(i),idust)
             end if
             end do
        end do
     end do
     end do
     end do
  end do
  !Update conservative variables new state vector
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     ind_son=1+i2+2*j2+4*k2
     iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        do i=1,ncache
           if(son(ind_cell(i))==0)then
              sum_dust_old=0.0_dp
              do idust=1,ndust
                 !We compute sum_dust_old in order to retrieve rho_gas
                 sum_dust_old=sum_dust_old+uold(ind_cell(i),firstindex_ndust+idust)
              enddo
              !we deduce rho_gas 
              rho_gas = uold(ind_cell(i),1)-sum_dust_old
            do idust=1,ndust
                 !Update epsilon taking in account small fluxes from refined interfaces
                 unew(ind_cell(i),firstindex_ndust+idust)=uuloc(i,i3,j3,k3,idust)+(uold(ind_cell(i),firstindex_ndust+idust)& 
                                                                                &+f_dust(ind_cell(i),idust))!/uold(ind_cell(i),1)
            enddo   
              sum_dust_new=0.0_dp              
              do idust=1,ndust
               !We compute sum_dust_new
               sum_dust_new = sum_dust_new + unew(ind_cell(i),firstindex_ndust+idust)!*uold(ind_cell(i),1)
              enddo
              if(dust_barr.eqv. .false.) then
                 !Update all the quantities that depend on rho
                 uuuloc(i,i3,j3,k3,1)= rho_gas + sum_dust_new! rho_gas*(1.0_dp + sum_dust_new)
          
               else
                 !If we test barenblatt we only update P
                 unew(ind_cell(i),5)=(1.0_dp-sum_dust_new)*uold(ind_cell(i),1)*kb*temp_dust/mu_gas/mh/scale_v**2/(gamma-1.0_dp)
               endif
           end if
        end do
  end do
  end do
end do

  !Update conservative variables new state vector
  do k2=k2min,k2max
  do j2=j2min,j2max
  do i2=i2min,i2max
     ind_son=1+i2+2*j2+4*k2
     iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        do i=1,ncache
          if(son(ind_cell(i))==0)then
          if(dust_barr.eqv. .false.) then
                 d = uuuloc(i,i3,j3,k3,1)
                 u = uuuloc(i,i3,j3,k3,2)
                 v = uuuloc(i,i3,j3,k3,3)
                 w = uuuloc(i,i3,j3,k3,4)
                 unew(ind_cell(i),1) = d
                 unew(ind_cell(i),2) = d*u
                 unew(ind_cell(i),3) = d*v
                 unew(ind_cell(i),4) = d*w
                 unew(ind_cell(i),5) = uuuloc(i,i3,j3,k3,5) + 0.5d0*d*(u**2+v**2+w**2)
                 !do idust=1,ndust
                 !Update rho epsilon 
                 !unew(ind_cell(i),firstindex_ndust+idust)=d*unew(ind_cell(i),firstindex_ndust+idust)
                 !enddo   
              end if
              end if
        end do
  end do
  end do
  end do

  if(ilevel>levelmin)then
  !-----------------------------------------------------
  ! update at level ilevel-1
  !-----------------------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     !----------------------
     ! Left flux at boundary
     !----------------------     
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh)   = i
        end if
     end do
     ! Conservative update of the flux
     do idust=1,ndust
        ! Loop over boundary cells
        do k3=k3min,k3max-k0
        do j3=j3min,j3max-j0
        do i3=i3min,i3max-i0
           do i=1,nb_noneigh
              f_dust(ind_buffer(i),idust)=f_dust(ind_buffer(i),idust) &
                   &-flux(ind_cell(i),i3,j3,k3,idust,idim)*oneontwotondim            
           end do
        end do
        end do
        end do
     end do
     
     !-----------------------
     ! Right flux at boundary
     !-----------------------     
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh)   = i
        end if
     end do
     ! Conservative update of the flux
     do idust=1,ndust
        ! Loop over boundary cells
        do k3=k3min+k0,k3max
        do j3=j3min+j0,j3max
        do i3=i3min+i0,i3max
           do i=1,nb_noneigh
              f_dust(ind_buffer(i),idust)=f_dust(ind_buffer(i),idust) &
                   &+flux(ind_cell(i),i3+i0,j3+j0,k3+k0,idust,idim)*oneontwotondim
           end do
        end do
        end do
        end do
     end do
  end do
  ! End loop over dimensions
  end if
end subroutine dustdifffine1
