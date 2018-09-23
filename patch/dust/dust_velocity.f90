
subroutine vdust_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the dust velocity scheme.
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

  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call vdustfine1(ind_grid,ngrid,ilevel)
  end do

  call upload_fine(ilevel)
  do idim =1,ndim
     do idust=1,ndust
       call make_virtual_fine_dp(v_dust(1,idust,idim),ilevel)
     end do
  end do
  if(simple_boundary)call make_boundary_hydro(ilevel)

111 format('   Entering vdust_fine for level ',i2)


end subroutine vdust_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine vdustfine1(ind_grid,ncache,ilevel)
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
  ! vdust routine that compute the velocity.                           -!
  !---------------------------------------------------------------------!

  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim),save::udust
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::upress
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
  real(dp):: dt_dust
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
  !Saved variables set to 0
  u1   = 0.0d0
  u2   = 0.0d0

 
  ok   = .false.
  oneontwotondim = 1.d0/dble(twotondim)

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
           if (epsil_cons) then
           do idust=1,ndust 
              do i=1,nbuffer
                 u1(i,j,firstindex_ndust+idust)=uold(ibuffer_father(i,j),firstindex_ndust+idust)/uold(ibuffer_father(i,j),1)
              end do
           end do
           end if 
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

        do i=1,nexist
           do ivar=1,nvar+3
              upress(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           if(epsil_cons) then
          do idust=1,ndust
              upress(ind_exist(i),i3,j3,k3,firstindex_ndust+idust)=uold(ind_cell(i),firstindex_ndust+idust)/uold(ind_cell(i),1)
           end do
           end if
        end do
        do i=1,nbuffer
           do ivar =1, nvar+3
              upress(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)                  
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
  ! Compute dust velocity
!-----------------------------------------------
call v_dust1(upress,udust,dx,ncache,ind_grid,dtnew(ilevel))


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
                 v_dust(ind_cell(i),idust,idim) =  udust(i,i3,j3,k3,idust,idim)
              enddo
           end if
     end do
  end do
  end do
  end do
end do



end subroutine vdustfine1




!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine v_dust1(upress,udust,dx,ngrid,ind_grid,dt)
  use amr_parameters
  use const             
  use hydro_parameters
  use amr_commons
  use hydro_commons
  use units_commons,ONLY:scale_m
  use cloud_module
  use cooling_module,ONLY:kB,mH
  use radiation_parameters

  implicit none
  integer,dimension(1:nvector)::ind_grid
  integer ::ngrid
  real(dp)::dx,dy,dz
  real(dp)::scale,d,u,v,w,eold,A,B,C,pressure,e_mag,e_kin,enint,d0,r0
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist

  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::upress
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim)::udust
  integer ,dimension(1:nvector                   ),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist
  integer::ind_son,ind_father,iskip,nbuffer,ibuffer
  
  real(dp),dimension(1:ndust)  :: t_stop
  real(dp)  ::cs,pi,tstop_tot,t_stop_floor,dens_floor
  real(dp), dimension(1:ndust) ::d_grain,l_grain
  ! Output fluxes

  ! Local scalar variables
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2 
  integer::i,j,k,l,ivar, idust,idim
  integer::ilo,ihi,jlo,jhi,klo,khi
  real(dp):: epsilon_0, dt_dust,wnorm,pr,pl,sum_dust,dt
  real(dp),dimension(1:ndust):: dustMRN
  epsilon_0 = dust_ratio(1)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  t_stop_floor=0.0d0
  sum_dust=0.0d0
  t_stop=0.0d0
  sum_dust=0.0d0
  pi =3.14159265358979323846_dp
  dens_floor=0.0d0  
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
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

 ! Compute the dust velocity X direction
  do i=ilo,ihi         
  do j=jlo,jhi
  do k=klo,khi
     do l = 1, ngrid
            tstop_tot=0.0d0
            t_stop=0.0d0
              if(energy_fix)then
            !compute right pressure
              d=upress(l,i+1,j,k,1)
              enint=upress(l,i+1,j,k,nvar)
           else
              d=upress(l,i+1,j,k,1)
              u=upress(l,i+1,j,k,2)/d
              v=upress(l,i+1,j,k,3)/d
              w=upress(l,i+1,j,k,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i+1,j,k,6)+upress(l,i+1,j,k,nvar+1))
              B=0.5d0*(upress(l,i+1,j,k,7)+upress(l,i+1,j,k,nvar+2))
              C=0.5d0*(upress(l,i+1,j,k,8)+upress(l,i+1,j,k,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i+1,j,k,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i+1,j,k,5)
              else   
                 enint=upress(l,i+1,j,k,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
             if (.not.epsil_cons)  sum_dust = sum_dust + upress(l,i+1,j,k,firstindex_ndust+idust)/d
              if (epsil_cons) sum_dust = sum_dust + upress(l,i+1,j,k,firstindex_ndust+idust)
 
           end do
           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pr)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           if(dust_barr) pr = (1.0_dp-sum_dust)*d*cs*cs
            !compute left pressure
              if(energy_fix)then
              d=upress(l,i-1,j,k,1)
              enint=upress(l,i-1,j,k,nvar)
           else
              d=upress(l,i-1,j,k,1)
              u=upress(l,i-1,j,k,2)/d
              v=upress(l,i-1,j,k,3)/d
              w=upress(l,i-1,j,k,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i-1,j,k,6)+upress(l,i-1,j,k,nvar+1))
              B=0.5d0*(upress(l,i-1,j,k,7)+upress(l,i-1,j,k,nvar+2))
              C=0.5d0*(upress(l,i-1,j,k,8)+upress(l,i-1,j,k,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i-1,j,k,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i-1,j,k,5)
              else   
                 enint=upress(l,i-1,j,k,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
              if (.not.epsil_cons)sum_dust = sum_dust + upress(l,i-1,j,k,firstindex_ndust+idust)/d
              if (epsil_cons) sum_dust = sum_dust + upress(l,i-1,j,k,firstindex_ndust+idust)

           end do
           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pl)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           if(dust_barr) pl = (1.0_dp-sum_dust)*d*cs*cs
             !compute  central sound speed
              if(energy_fix)then
              d=upress(l,i,j,k,1)
              enint=upress(l,i,j,k,nvar)
           else
              d=upress(l,i,j,k,1)
              u=upress(l,i,j,k,2)/d
              v=upress(l,i,j,k,3)/d
              w=upress(l,i,j,k,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i,j,k,6)+upress(l,i,j,k,nvar+1))
              B=0.5d0*(upress(l,i,j,k,7)+upress(l,i,j,k,nvar+2))
              C=0.5d0*(upress(l,i,j,k,8)+upress(l,i,j,k,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i,j,k,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i,j,k,5)
              else   
                 enint=upress(l,i,j,k,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
              if (.not.epsil_cons)sum_dust = sum_dust + upress(l,i,j,k,firstindex_ndust+idust)/d
              if (epsil_cons) sum_dust = sum_dust + upress(l,i,j,k,firstindex_ndust+idust)

           end do
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           
            do idust = 1,ndust
               t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d/(1.0d0-sum_dust)
               if(K_drag)  t_stop(idust) = upress(l,i,j,k,firstindex_ndust+idust)/K_dust(idust)
               if(dust_barr) t_stop (idust)= 0.1_dp
                if (.not.epsil_cons) tstop_tot= tstop_tot-t_stop(idust)*(upress(l,i,j,k,firstindex_ndust+idust)/d)
               if (epsil_cons)  tstop_tot= tstop_tot-t_stop(idust)*upress(l,i,j,k,firstindex_ndust+idust)

            end do
            do idust = 1,ndust
               t_stop(idust) = t_stop(idust)+tstop_tot
               udust(l,i,j,k,idust,1)=t_stop(idust)*(pr-pl)/(2.0d0*dx)/d  

          end do
    enddo
  enddo
  enddo
  enddo
#if NDIM>1
 ! Compute the dust velocity y direction

  do i=ilo,ihi         
  do j=jlo,jhi
  do k=klo,khi
     do l = 1, ngrid

            tstop_tot=0.0d0
            t_stop=0.0d0
              if(energy_fix)then
            !compute right pressure
              d=upress(l,i,j+1,k,1)
              enint=upress(l,i,j+1,k,nvar)
           else
              d=upress(l,i,j+1,k,1)
              u=upress(l,i,j+1,k,2)/d
              v=upress(l,i,j+1,k,3)/d
              w=upress(l,i,j+1,k,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i,j+1,k,6)+upress(l,i,j+1,k,nvar+1))
              B=0.5d0*(upress(l,i,j+1,k,7)+upress(l,i,j+1,k,nvar+2))
              C=0.5d0*(upress(l,i,j+1,k,8)+upress(l,i,j+1,k,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i,j+1,k,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i,j+1,k,5)
              else   
                 enint=upress(l,i,j+1,k,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
              if (.not.epsil_cons)sum_dust = sum_dust + upress(l,i,j+1,k,firstindex_ndust+idust)/d
              if (epsil_cons) sum_dust = sum_dust + upress(l,i,j+1,k,firstindex_ndust+idust)
                            
           end do
           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pr)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           if(dust_barr) pr = (1.0_dp-sum_dust)*d*cs*cs
            !compute left pressure
              if(energy_fix)then
                 d=upress(l,i,j-1,k,1)
              enint=upress(l,i,j-1,k,nvar)
           else
              d=upress(l,i,j-1,k,1)
              u=upress(l,i,j-1,k,2)/d
              v=upress(l,i,j-1,k,3)/d
              w=upress(l,i,j-1,k,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i,j-1,k,6)+upress(l,i,j-1,k,nvar+1))
              B=0.5d0*(upress(l,i,j-1,k,7)+upress(l,i,j-1,k,nvar+2))
              C=0.5d0*(upress(l,i,j-1,k,8)+upress(l,i,j-1,k,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i,j-1,k,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i,j-1,k,5)
              else   
                 enint=upress(l,i,j-1,k,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
              if (.not.epsil_cons)sum_dust = sum_dust + upress(l,i,j-1,k,firstindex_ndust+idust)/d
                            
           end do
           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pl)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           if(dust_barr) pl = (1.0_dp-sum_dust)*d*cs*cs
             !compute  central sound speed
              if(energy_fix)then
              d=upress(l,i,j,k,1)
              enint=upress(l,i,j,k,nvar)
           else
              d=upress(l,i,j,k,1)
              u=upress(l,i,j,k,2)/d
              v=upress(l,i,j,k,3)/d
              w=upress(l,i,j,k,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i,j,k,6)+upress(l,i,j,k,nvar+1))
              B=0.5d0*(upress(l,i,j,k,7)+upress(l,i,j,k,nvar+2))
              C=0.5d0*(upress(l,i,j,k,8)+upress(l,i,j,k,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i,j,k,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i,j,k,5)
              else   
                 enint=upress(l,i,j,k,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
              if (.not.epsil_cons)sum_dust = sum_dust + upress(l,i,j,k,firstindex_ndust+idust)/d
              if (epsil_cons) sum_dust = sum_dust + upress(l,i,j,k,firstindex_ndust+idust)

           end do
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           
            do idust = 1,ndust
               t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d/(1.0d0-sum_dust)
               if(K_drag)  t_stop(idust) = upress(l,i,j,k,firstindex_ndust+idust)/K_dust(idust)
               if(dust_barr) t_stop (idust)= 0.1_dp
                if (.not.epsil_cons) tstop_tot= tstop_tot-t_stop(idust)*(upress(l,i,j,k,firstindex_ndust+idust)/d)
               if (epsil_cons)  tstop_tot= tstop_tot-t_stop(idust)*upress(l,i,j,k,firstindex_ndust+idust)
            end do
            do idust = 1,ndust
               t_stop(idust) = t_stop(idust)+tstop_tot
               if (.not.epsil_cons)udust(l,i,j,k,idust,2)=t_stop(idust)*(pr-pl)/(2.0d0*dx)/d
               if (epsil_cons)  udust(l,i,j,k,idust,2)=t_stop(idust)*(pr-pl)/(2.0d0*dx)

               !if(d .le. dens_floor)  udust(l,i,j,k,idust,2)=0.0d0
          end do
    enddo
  enddo
  enddo
  enddo
#endif 

#if NDIM>2
 ! Compute the dust velocity z direction

  do i=ilo,ihi         
  do j=jlo,jhi
  do k=klo,khi      
     do l = 1, ngrid

            tstop_tot=0.0d0
            t_stop=0.0d0
              if(energy_fix)then
            !compute right pressure
              d=upress(l,i,j,k+1,1)
              enint=upress(l,i,j,k+1,nvar)
           else
              d=upress(l,i,j,k+1,1)
              u=upress(l,i,j,k+1,2)/d
              v=upress(l,i,j,k+1,3)/d
              w=upress(l,i,j,k+1,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i,j,k+1,6)+upress(l,i,j,k+1,nvar+1))
              B=0.5d0*(upress(l,i,j,k+1,7)+upress(l,i,j,k+1,nvar+2))
              C=0.5d0*(upress(l,i,j,k+1,8)+upress(l,i,j,k+1,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i,j,k+1,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i,j,k+1,5)
              else   
                 enint=upress(l,i,j,k+1,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
              if (.not.epsil_cons)sum_dust = sum_dust + upress(l,i,j,k+1,firstindex_ndust+idust)/d
              if (epsil_cons) sum_dust = sum_dust + upress(l,i,j,k+1,firstindex_ndust+idust)

           end do
           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pr)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           if(dust_barr) pr = (1.0_dp-sum_dust)*d*cs*cs
            !compute left pressure
              if(energy_fix)then
                 d=upress(l,i,j,k-1,1)
              enint=upress(l,i,j,k-1,nvar)
           else
              d=upress(l,i,j,k-1,1)
              u=upress(l,i,j,k-1,2)/d
              v=upress(l,i,j,k-1,3)/d
              w=upress(l,i,j,k-1,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i,j,k-1,6)+upress(l,i,j,k-1,nvar+1))
              B=0.5d0*(upress(l,i,j,k-1,7)+upress(l,i,j,k-1,nvar+2))
              C=0.5d0*(upress(l,i,j,k-1,8)+upress(l,i,j,k-1,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i,j,k-1,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i,j,k-1,5)
              else   
                 enint=upress(l,i,j,k-1,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
              if (.not.epsil_cons)sum_dust = sum_dust + upress(l,i,j,k-1,firstindex_ndust+idust)/d
              if (epsil_cons) sum_dust = sum_dust + upress(l,i,j,k-1,firstindex_ndust+idust)

           end do
           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pl)
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           if(dust_barr) pl = (1.0_dp-sum_dust)*d*cs*cs
             !compute  central sound speed
              if(energy_fix)then
              d=upress(l,i,j,k,1)
              enint=upress(l,i,j,k,nvar)
           else
              d=upress(l,i,j,k,1)
              u=upress(l,i,j,k,2)/d
              v=upress(l,i,j,k,3)/d
              w=upress(l,i,j,k,4)/d
              e_mag= 0.0_dp
#ifdef SOLVERmhd                           
              A=0.5d0*(upress(l,i,j,k,6)+upress(l,i,j,k,nvar+1))
              B=0.5d0*(upress(l,i,j,k,7)+upress(l,i,j,k,nvar+2))
              C=0.5d0*(upress(l,i,j,k,8)+upress(l,i,j,k,nvar+3))
              e_mag=0.5d0*(A**2+B**2+C**2)
#endif
              e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                 e_mag=e_mag+upress(l,i,j,k,8+irad)
              end do
#endif
              if(dust_barr) then
                 enint=upress(l,i,j,k,5)
              else   
                 enint=upress(l,i,j,k,5)-e_kin- e_mag
              endif
           endif
           sum_dust=0.0_dp
           do idust = 1, ndust
              if (.not.epsil_cons) sum_dust = sum_dust + upress(l,i,j,k,firstindex_ndust+idust)/d
              if (epsil_cons) sum_dust = sum_dust + upress(l,i,j,k,firstindex_ndust+idust)

           end do
           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
           if(dust_barr)  cs = 1.0_dp
           
            do idust = 1,ndust
               t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d/(1.0d0-sum_dust)
               if(K_drag)  t_stop(idust) = upress(l,i,j,k,firstindex_ndust+idust)/K_dust(idust)
               if(dust_barr) t_stop (idust)= 0.1_dp
               if (.not.epsil_cons)  tstop_tot= tstop_tot-t_stop(idust)*(upress(l,i,j,k,firstindex_ndust+idust)/d)
               if (epsil_cons)  tstop_tot= tstop_tot-t_stop(idust)*upress(l,i,j,k,firstindex_ndust+idust)

            end do
            do idust = 1,ndust
               t_stop(idust) = t_stop(idust)+tstop_tot
               udust(l,i,j,k,idust,3)=t_stop(idust)*(pr-pl)/(2.0d0*dx)/d 
              
          end do
    enddo
  enddo
  enddo
  enddo
#endif

if (reduce_wdust .eqv. .true.) then
  do i=ilo,ihi         
  do j=jlo,jhi
  do k=klo,khi      
     do l = 1, ngrid
          
            do idust = 1,ndust
              if (NDIM.eq.1)  dt_dust   = courant_factor*dx/(abs( udust(l,i,j,k,idust,1)))
              if (NDIM.eq.2)  dt_dust  = courant_factor*dx/(abs( udust(l,i,j,k,idust,1))+abs( udust(l,i,j,k,idust,2)))
              if (NDIM.eq.3)  dt_dust  = courant_factor*dx/(abs(udust(l,i,j,k,idust,1))+abs( udust(l,i,j,k,idust,2))+abs(udust(l,i,j,k,idust,3)))
 
              if( dt_dust .le.dt.or.+upress(l,i,j,k,1) .lt. d0 .and. reduce_wdust .eqv. .true.) then   
	      if (NDIM.eq.1) wnorm =abs( udust(l,i,j,k,idust,1))
 	      if (NDIM.eq.2) wnorm =sqrt( udust(l,i,j,k,idust,1)**2.0+ udust(l,i,j,k,idust,2)**2.0)
	      if (NDIM.eq.3) wnorm =sqrt(udust(l,i,j,k,idust,1)**2.0+udust(l,i,j,k,idust,2)**2.0+udust(l,i,j,k,idust,3)**2.0)
              do idim=1,ndim       
             !udust(l,i,j,k,idust,idim) = sign(beta_dust*courant_factor*dx/dt*abs(udust(l,i,j,k,idust,idim))/wnorm,udust(l,i,j,k,idust,idim))
               end do   
            end if   
          end do
    enddo
  enddo
  enddo
  enddo
end if





                            


end subroutine v_dust1



