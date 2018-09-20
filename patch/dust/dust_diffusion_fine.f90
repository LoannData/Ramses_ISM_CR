subroutine set_dflux_dust_new(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets the dust fluxes to 0 before calling
  ! the diffusion scheme. 
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,irad,idust,idim
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  ! Set fdust = 0 for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do idust=1,ndust
        do i=1,active(ilevel)%ngrid
           dflux_dust(active(ilevel)%igrid(i)+iskip,idust) = 0.0d0
        end do
     end do
  end do
  ! Set fdust to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do idust=1,ndust
        do i=1,reception(icpu,ilevel)%ngrid
          dflux_dust(reception(icpu,ilevel)%igrid(i)+iskip,idust)= 0.0d0
        end do
     end do
  end do
end do
111 format('   Entering set_dflux_dust for level ',i2)

end subroutine set_dflux_dust_new
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
  ! This routine sets array unew =uold for the quantities that are updated because of dust diffusion.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,irad,idust

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
        do i=1,active(ilevel)%ngrid
            !unew(active(ilevel)%igrid(i)+iskip,5) = uold(active(ilevel)%igrid(i)+iskip,5)
           do idust=1,ndust
              unew(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust) = uold(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust)+dflux_dust(active(ilevel)%igrid(i)+iskip,idust)
           end do
        end do
     end do
  ! Set unew to 0 for virtual boundary cells
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do idust=1,ndust
      do i=1,reception(icpu,ilevel)%ngrid
         unew(reception(icpu,ilevel)%igrid(i)+iskip,firstindex_ndust+idust)=0.0d0
       end do
     end do
     do i=1,reception(icpu,ilevel)%ngrid
        !unew(reception(icpu,ilevel)%igrid(i)+iskip,5)=0.0d0
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
  integer,dimension(1:nvector)::ind_grid
  
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array uold =unew for the quantities that have been updated because of dust diffusion.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,irad,idust,ht,ncache
  real (dp):: rho_gas ,sum_dust_old,sum_dust_new, d, u ,v ,w, enint,temp, e_mag, e_kin, A,B,C
  
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
 
  !Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid
           d = uold(active(ilevel)%igrid(i)+iskip,1)
           sum_dust_old=0.0_dp
           do idust=1,ndust
               !We compute the old dust density
               sum_dust_old=sum_dust_old+uold(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust)/d
           enddo
           sum_dust_new=0.0_dp
           do idust=1,ndust
!              !We compute the old dust density
              sum_dust_new=sum_dust_new+unew(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust)/d
           enddo
           !Update all the quantities that depend on rho
           rho_gas = uold(active(ilevel)%igrid(i)+iskip,1)-sum_dust_old*d
           u = uold(active(ilevel)%igrid(i)+iskip,2)/d
           v = uold(active(ilevel)%igrid(i)+iskip,3)/d
           w = uold(active(ilevel)%igrid(i)+iskip,4)/d
           e_mag= 0.0_dp
#ifdef SOLVERmhd                           
           A=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,6)+uold(active(ilevel)%igrid(i)+iskip,nvar+1))
           B=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,7)+uold(active(ilevel)%igrid(i)+iskip,nvar+2))
           C=0.5d0*(uold(active(ilevel)%igrid(i)+iskip,8)+uold(active(ilevel)%igrid(i)+iskip,nvar+3))
           e_mag=0.5d0*(A**2+B**2+C**2)
#endif
           e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e_mag=e_mag+uold(active(ilevel)%igrid(i)+iskip,8+irad)
           end do
#endif
           enint=0.0d0
           call temperature_eos(rho_gas, uold(active(ilevel)%igrid(i)+iskip,5) -e_kin -e_mag , temp, ht,sum_dust_new )
           rho_gas =  uold(active(ilevel)%igrid(i)+iskip,1)-sum_dust_new*d
           call enerint_eos (rho_gas, temp , enint)
           unew(active(ilevel)%igrid(i)+iskip,5) = enint + e_kin +e_mag
           !If we test barenblatt we only update P
           if(static_gas) unew(active(ilevel)%igrid(i)+iskip,5)=(1.0_dp-sum_dust_new)*uold(active(ilevel)%igrid(i)+iskip,1)/(gamma-1.0_dp)
           uold(active(ilevel)%igrid(i)+iskip,5) = unew(active(ilevel)%igrid(i)+iskip,5)
        end do
     end do
  !Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid           
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
   
subroutine dust_diffusion_fine(ilevel,d_cycle_ok,ncycle,icycle)
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
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,idust,idim
  integer,dimension(1:nvector),save::ind_grid
  logical:: d_cycle_ok
  integer :: icycle, ncycle

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  do ivar =1, nvar
    call make_virtual_fine_dp(uold(1,ivar),ilevel)
  end do
  if (.not.mhd_dust)call set_vdust(ilevel)
  if (mhd_dust) call set_vdust_mhd(ilevel)
  call upload_fine(ilevel)

  do idim =1,ndim
     do idust=1,ndust
        call make_virtual_fine_dp(v_dust(1,idust,idim),ilevel)
     end do
  end do

  call set_unew_dust(ilevel)
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call dustdifffine1(ind_grid,ngrid,ilevel,d_cycle_ok,ncycle,icycle)
  end do
  do idust=1,ndust
     call make_virtual_reverse_dp(unew(1,firstindex_ndust+idust),ilevel)
  end do
  call set_uold_dust(ilevel)
  do idust=1,ndust
     call make_virtual_reverse_dp(dflux_dust(1,idust),ilevel)
  end do

  if (.not.mhd_dust)call set_vdust(ilevel)
  if (mhd_dust) call set_vdust_mhd(ilevel)
  call upload_fine(ilevel)
  do idust=1,ndust
     call make_virtual_fine_dp(uold(1,firstindex_ndust+idust),ilevel)
  end do
  call make_virtual_fine_dp(uold(1,5),ilevel)
  do idim =1,ndim
     do idust=1,ndust
       call make_virtual_fine_dp(v_dust(1,idust,idim),ilevel)
     end do
  end do
  if(simple_boundary)call make_boundary_hydro(ilevel)

111 format('   Entering dust_diffusion_fine for level ',i2)


end subroutine dust_diffusion_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine dustdifffine1(ind_grid,ncache,ilevel,d_cycle_ok,ncycle,icycle)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use cooling_module,ONLY:kB,mH
  use cloud_module
  use radiation_parameters
  use units_commons, only : scale_m


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
  real(dp),dimension(1:nvector,0:twondim  ,1:ndust*ndim),save::u1dust
  real(dp),dimension(1:nvector,1:twotondim,1:ndust*ndim),save::u2dust
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust+ndust*ndim),save::uloc
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust,1:ndim),save::flux
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  integer ,dimension(1:nvector                   ),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist
  real(dp),dimension(1:twotondim,1:3),save::xc
  real(dp),dimension(1:twotondim),save::rc

  integer::idust,ht
  integer::i,j,ivar,idim,irad,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  integer::  ncycle,icycle
  real(dp):: dt_dustcycle
  logical :: d_cycle_ok
  real(dp)::dx,scale,oneontwotondim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::sum_dust,sum_dust_new,sum_dust_old
  real(dp)::d,u,v,w,A,B,C,enint,e_kin,e_mag,pressure,cs, temp
  real(dp)::rho_gas, pi, t_stop,t_stop_floor,dens_floor,d0,r0
  real(dp), dimension(1:ndust) ::d_grain,l_grain
  real(dp):: epsilon_0
  real(dp),dimension(1:ndust):: dustMRN
  epsilon_0 = dust_ratio(1)  
  pi =3.14159265358979323846_dp
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  t_stop_floor =0.0d0
  sum_dust=0.0d0

  !Saved variables set to 0
  u1   = 0.0d0
  u2   = 0.0d0
  u1dust=0.0d0
  u2dust=0.0d0
  flux = 0.0d0  
  uloc = 0.0d0
  ok   = .false.
  oneontwotondim = 1.d0/dble(twotondim)
  
  ! Initialisation of dust related quantities
  
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
           do idim=1,ndim
              do idust=1,ndust
                 do i=1,nbuffer
                    u1dust(i,j,ndim*(idust-1)+idim)=v_dust(ibuffer_father(i,j),idust,idim)
                 end do
              end do
           end do
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
        call interpol_hydro_dust(u1dust,ind1,u2dust,nbuffer)

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

        do i=1,nexist
           do idust=1,ndust
              uloc(ind_exist(i),i3,j3,k3,idust)=uold(ind_cell(i),firstindex_ndust+idust)
            end do
              do idim= 1,ndim
                 do idust=1,ndust
                    uloc(ind_exist(i),i3,j3,k3,ndust+ndim*(idust-1)+idim)= v_dust(ind_cell(i),idust,idim)
                 end do
                 !print*, v_dust(ind_cell(i),idust,idim), ndust+ndim*(idust-1)+idim, idim
               end do
            end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,idust)=u2(i,ind_son,firstindex_ndust+idust)
              do idim= 1,ndim
                 do idust=1,ndust
                    uloc(ind_nexist(i),i3,j3,k3,ndust+ndim*(idust-1)+idim)= u2dust(i,ind_son,ndim*(idust-1)+idim)
                 end do
              end do
           end do
        end do

     end do
  end do
end do
     !End loop over cells
end do
end do
!end do
!End loop over neighboring grids
  
  !-----------------------------------------------
  ! Compute flux due to dust diffusion
  !-----------------------------------------------

  !call dustdiff_split(uloc,flux,dx,dx,dx,dtnew(ilevel),ncache)
  call dustdiff_predict(uloc,flux,dx,dx,dx,dtnew(ilevel),ncache)
  
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
  do idim=1,ndim
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
        do i=1,ncache
           if(son(ind_cell(i))==0)then
              do idust=1,ndust
                 !Update rhodust
                ! print *, unew(ind_cell(i),firstindex_ndust+idust)
                 unew(ind_cell(i),firstindex_ndust+idust)=unew(ind_cell(i),firstindex_ndust+idust) +(flux(i,i3,j3,k3,idust,idim)&
                      &-flux(i,i3+i0,j3+j0,k3+k0,idust,idim))
                 !print *, flux(i,i3,j3,k3,idust,idim)&
                 !     &,-flux(i,i3+i0,j3+j0,k3+k0,idust,idim), unew(ind_cell(i),firstindex_ndust+idust)
              enddo
           end if
     end do
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
              dflux_dust(ind_buffer(i),idust)=dflux_dust(ind_buffer(i),idust) &
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
              dflux_dust(ind_buffer(i),idust)=dflux_dust(ind_buffer(i),idust) &
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


