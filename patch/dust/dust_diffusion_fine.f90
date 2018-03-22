subroutine set_vdust(ilevel)
  use amr_commons
  use hydro_commons
  use units_commons
  use cloud_module
  use cooling_module,ONLY:kB,mH
  use radiation_parameters

  implicit none
  integer::ilevel
  integer::i,j,k,ivar,irad,ind,iskip,nx_loc,ind_cell1,idust
  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::scale,dx,dx_loc,d,u,v,w,eold,A,B,C,pressure

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim),save::igridn
  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim),save::dx_g,dx_d
  real(dp)::usquare,emag,erad_loc,ekin,eps,sum_dust,enint
  real(dp)::e_mag,e_kin,e_cons,e_prim,e_trunc,div,fact,e_r
  real(dp)::Pgdivu,u_square,d_loc,Tp_loc,Tr_loc,cal_Teg
  real(dp),dimension(1:nvector,1:ndim),save::Pleft,Pright
  real(dp),dimension(1:nvector,1:ndim)       ::gradP
  real(dp),dimension(1:ndust)  :: t_stop
  real(dp)  ::cs,pi,tstop_tot,t_stop_floor,dens_floor
  real(dp), dimension(1:ndust) ::d_grain,l_grain
  real(dp) :: dd,ee,cmp_Cv_eos,d0,r0
  integer  :: ht
  real(dp):: epsilon_0
  real(dp),dimension(1:ndust):: dustMRN
  epsilon_0 = dust_ratio(1)
 
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(numbtot(1,ilevel)==0)return
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale
  t_stop_floor=0.0d0
  sum_dust=0.0d0
  Pleft=0.0; Pright=0.0
  t_stop=0.0d0
  sum_dust=0.0d0
  pi =3.14159265358979323846_dp

#if NDUST>0
     do idust =1,ndust
        dustMRN(idust) = dust_ratio(idust)/(1.0d0+dust_ratio(idust))
     end do     
     if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
     do idust =1,ndust
           sum_dust = sum_dust + dustMRN(idust)
        end do   
#endif   
  r0=(alpha_dense_core*2.*6.67d-8*mass_c*scale_m*mu_gas*mH/(5.*kB*Tr_floor*(1.0d0-sum_dust)))/scale_l
  d0 = 3.0d0*mass_c/(4.0d0*pi*r0**3.)
  dens_floor=d0
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
  
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
   
     ! Gather nvector grids
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Gather neighboring grids
     do i=1,ngrid
        igridn(i,0)=ind_grid(i)
     end do
     do idim=1,ndim
        do i=1,ngrid
           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
           igridn(i,2*idim-1)=son(ind_left (i,idim))
           igridn(i,2*idim  )=son(ind_right(i,idim))
        end do
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Compute central cell index
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather all neighboring velocities
        do idim=1,ndim
           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
           ih1=ncoarse+(id1-1)*ngridmax
           do i=1,ngrid
              dx_d(i,idim)=dx_loc

              if(energy_fix)then
               eold=uold(ind_left(i,idim),nvar)
              else
              ! Gather left thermal energy
              d=max(uold(ind_left(i,idim),1),smallr)
              u=0.0; v=0.0; w=0.0
              if(ndim>0)u=uold(ind_left(i,idim),2)/d
              if(ndim>1)v=uold(ind_left(i,idim),3)/d
              if(ndim>2)w=uold(ind_left(i,idim),4)/d
              A=0.5d0*(uold(ind_left(i,idim),6)+uold(ind_left(i,idim),nvar+1))
              B=0.5d0*(uold(ind_left(i,idim),7)+uold(ind_left(i,idim),nvar+2))
              C=0.5d0*(uold(ind_left(i,idim),8)+uold(ind_left(i,idim),nvar+3))
              eold=uold(ind_left(i,idim),5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)
#if NENER>0
              do irad=1,nener
                 eold=eold-uold(ind_left(i,idim),8+irad)
              end do
#endif
              endif    
              sum_dust=0.0d0
              do idust = 1, Ndust
                 sum_dust=sum_dust+uold(ind_left(i,idim),firstindex_ndust+idust)/d
              end do
              call pressure_eos((1.0_dp-sum_dust)*d,eold,Pleft(i,idim))
      enddo
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid             
              dx_g(i,idim)=dx_loc

              if(energy_fix)then
              eold=uold(ind_left(i,idim),nvar)
              else
              ! Gather right thermal energy
              d=max(uold(ind_right(i,idim),1),smallr)
              u=0.0; v=0.0; w=0.0
              if(ndim>0)u=uold(ind_right(i,idim),2)/d
              if(ndim>1)v=uold(ind_right(i,idim),3)/d
              if(ndim>2)w=uold(ind_right(i,idim),4)/d
              A=0.5d0*(uold(ind_right(i,idim),6)+uold(ind_right(i,idim),nvar+1))
              B=0.5d0*(uold(ind_right(i,idim),7)+uold(ind_right(i,idim),nvar+2))
              C=0.5d0*(uold(ind_right(i,idim),8)+uold(ind_right(i,idim),nvar+3))
              eold=uold(ind_right(i,idim),5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)
#if NENER>0
              do irad=1,nener
                 eold=eold-uold(ind_right(i,idim),8+irad)
              end do
#endif
              endif                  
              sum_dust=0.0d0
              do idust = 1, Ndust
                 sum_dust=sum_dust+uold(ind_right(i,idim),firstindex_ndust+idust)/d
              end do
              call pressure_eos((1.0_dp-sum_dust)*d,eold,Pright(i,idim))
           enddo
        end do
        do i=1,ngrid
           do j=1,ndim
            gradP(i,j) = (Pright(i,j)-Pleft(i,j))/(dx_g(i,j)+dx_d(i,j))
           enddo
        end do
        
        do i=1,ngrid
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
              sum_dust = sum_dust + uold(ind_cell(i),firstindex_ndust+idust)/d
           end do
           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pressure)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           if(dust_barr) pressure = (1.0_dp-sum_dust)*d*cs*cs
            sum_dust=0.0d0
            do idust = 1, Ndust
               sum_dust=sum_dust+uold(ind_cell(i),firstindex_ndust+idust)/d
            end do
            tstop_tot=0.0d0            
            do idust = 1,ndust
               t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d/(1.0d0-sum_dust)
               if (d .le. dens_floor) t_stop(idust) =t_stop_floor 
               tstop_tot= tstop_tot-t_stop(idust)*(uold(ind_cell(i),firstindex_ndust+idust)/d)
            end do
            do idust = 1,ndust
               t_stop(idust) = t_stop(idust)+tstop_tot
            do idim=1,ndim
               v_dust(ind_cell(i),idust,idim) = t_stop(idust)* gradP(i,idim)/d
            end do
            end do
        end do
     enddo
     ! End loop over cells
  end do
  ! End loop over grids

end subroutine set_vdust

!###########################################################
!###########################################################
!###########################################################
!###########################################################
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
  do idim=1,ndim
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do idust=1,ndust
        do i=1,active(ilevel)%ngrid
           dflux_dust(active(ilevel)%igrid(i)+iskip,idust,idim) = 0.0d0
        end do
     end do
  end do
  end do
  ! Set fdust to 0 for virtual boundary cells
  do idim=1,ndim
  do icpu=1,ncpu
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do idust=1,ndust
        do i=1,reception(icpu,ilevel)%ngrid
          dflux_dust(reception(icpu,ilevel)%igrid(i)+iskip,idust,idim)= 0.0d0
        end do
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
              unew(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust) = uold(active(ilevel)%igrid(i)+iskip,firstindex_ndust+idust)
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
!           if(dust_barr.eqv. .false..and.isothermal) then
              !Update all the quantities that depend on rho
              rho_gas = uold(active(ilevel)%igrid(i)+iskip,1)-sum_dust_old*d
                 d = uold(active(ilevel)%igrid(i)+iskip,1)
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
                 call temperature_eos(rho_gas, uold(active(ilevel)%igrid(i)+iskip,5) -e_kin -e_mag , temp, ht )
                 rho_gas =  uold(active(ilevel)%igrid(i)+iskip,1)-sum_dust_new*d
                 call enerint_eos (rho_gas, temp , enint)
!                
                 unew(active(ilevel)%igrid(i)+iskip,5) = enint + e_kin +e_mag
!              else
!                 !If we test barenblatt we only update P
!                  unew(active(ilevel)%igrid(i)+iskip,5)=(1.0_dp-sum_dust_new)*uold(active(ilevel)%igrid(i)+iskip,1)/(gamma-1.0_dp)
!           endif
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
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,idust
  integer,dimension(1:nvector),save::ind_grid
  logical:: d_cycle_ok
  integer :: icycle, ncycle

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        call dustdifffine1(ind_grid,ngrid,ilevel,d_cycle_ok,ncycle,icycle)
     end do


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
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust+2),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::facdx
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
#if NDUST>0
     do idust =1,ndust
        dustMRN(idust) = dust_ratio(idust)/(1.0d0+dust_ratio(idust))
     end do     
     if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
     do idust =1,ndust
           sum_dust = sum_dust + dustMRN(idust)
        end do   
#endif
        
  r0=(alpha_dense_core*2.*6.67d-8*mass_c*scale_m*mu_gas*mH/(5.*kB*Tr_floor*(1.0d0-sum_dust)))/scale_l
  d0 = 3.0d0*mass_c/(4.0d0*pi*r0**3.)
  dens_floor=d0

  !Saved variables set to 0
  u1   = 0.0d0
  u2   = 0.0d0
  flux = 0.0d0  
  uloc = 0.0d0
  facdx= 0.0d0
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

           !We fill uloc with the quantities required to compute the fluxes (d, P, epsilon, ts/d)
           uloc(ind_exist(i),i3,j3,k3,2*ndust+1)=d
           uloc(ind_exist(i),i3,j3,k3,2*ndust+2)=pressure
           do idust = 1, ndust
              ! different prescriptions for t-stop
              t_stop = d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d!/d
              if(K_drag) t_stop = sum_dust*(1.0_dp-sum_dust)*d/K_dust(idust)
              if(dust_barr) t_stop = 0.1_dp
              if (d .le. dens_floor) t_stop =t_stop_floor 

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
              if (d .le. dens_floor) t_stop =t_stop_floor 
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
                 unew(ind_cell(i),firstindex_ndust+idust)=unew(ind_cell(i),firstindex_ndust+idust) +(flux(i,i3,j3,k3,idust,idim)&
                      &-flux(i,i3+i0,j3+j0,k3+k0,idust,idim))+dflux_dust(ind_cell(i),idust,idim)
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
              dflux_dust(ind_buffer(i),idust,idim)=dflux_dust(ind_buffer(i),idust,idim) &
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
              dflux_dust(ind_buffer(i),idust,idim)=dflux_dust(ind_buffer(i),idust,idim) &
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


