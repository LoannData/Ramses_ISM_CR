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
  if(verbose)write(*,111)ilevel
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
!  dens_floor=d0
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
              if(igridn(i,ig1)>0)then
                 dx_g(i,idim)=dx_loc              
              if(energy_fix)then
                 eold=uold(igridn(i,ig1)+ih1,nvar)
              else
                 ! Gather left thermal energy
                 d=max(uold(igridn(i,ig1)+ih1,1),smallr)
                 u=0.0; v=0.0; w=0.0
                 if(ndim>0)u=uold(igridn(i,ig1)+ih1,2)/d
                 if(ndim>1)v=uold(igridn(i,ig1)+ih1,3)/d
                 if(ndim>2)w=uold(igridn(i,ig1)+ih1,4)/d
                 A=0.5d0*(uold(igridn(i,ig1)+ih1,6)+uold(igridn(i,ig1)+ih1,nvar+1))
                 B=0.5d0*(uold(igridn(i,ig1)+ih1,7)+uold(igridn(i,ig1)+ih1,nvar+2))
                 C=0.5d0*(uold(igridn(i,ig1)+ih1,8)+uold(igridn(i,ig1)+ih1,nvar+3))
                 eold=uold(igridn(i,ig1)+ih1,5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)
#if NENER>0
                 do irad=1,nener
                    eold=eold-uold(igridn(i,ig1)+ih1,8+irad)
                 end do
#endif
              endif
              sum_dust=0.0d0
              do idust = 1, Ndust
                 sum_dust=sum_dust+uold(igridn(i,ig1)+ih1,firstindex_ndust+idust)/d
              end do
              call pressure_eos((1.0_dp-sum_dust)*d,eold,Pleft(i,idim))
           else
              dx_g(i,idim)=dx_loc*1.5d0         
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
         endif
      end do
   end do
   do idim=1,ndim
           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
           ih2=ncoarse+(id2-1)*ngridmax
           do i=1,ngrid
              if(igridn(i,ig2)>0)then
              dx_d(i,idim)=dx_loc              
              if(energy_fix)then
               eold=uold(igridn(i,ig2)+ih2,nvar)
              else
              ! Gather right thermal energy
              d=max(uold(igridn(i,ig2)+ih2,1),smallr)
              u=0.0; v=0.0; w=0.0
              if(ndim>0)u=uold(igridn(i,ig2)+ih2,2)/d
              if(ndim>1)v=uold(igridn(i,ig2)+ih2,3)/d
              if(ndim>2)w=uold(igridn(i,ig2)+ih2,4)/d
              A=0.5d0*(uold(igridn(i,ig2)+ih2,6)+uold(igridn(i,ig2)+ih2,nvar+1))
              B=0.5d0*(uold(igridn(i,ig2)+ih2,7)+uold(igridn(i,ig2)+ih2,nvar+2))
              C=0.5d0*(uold(igridn(i,ig2)+ih2,8)+uold(igridn(i,ig2)+ih2,nvar+3))
              eold=uold(igridn(i,ig2)+ih2,5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)
#if NENER>0
              do irad=1,nener
                 eold=eold-uold(igridn(i,ig2)+ih2,8+irad)
              end do
#endif
              endif    
              sum_dust=0.0d0
              do idust = 1, Ndust
                 sum_dust=sum_dust+uold(igridn(i,ig2)+ih2,firstindex_ndust+idust)/d
              end do
              call pressure_eos((1.0_dp-sum_dust)*d,eold,Pright(i,idim))
           else
              dx_d(i,idim)=dx_loc*1.5
              if(energy_fix)then
              eold=uold(ind_right(i,idim),nvar)
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
           endif
        end do
     end do
     do idim=1,ndim
           do i=1,ngrid
              gradP(i,idim) = (Pright(i,idim)-Pleft(i,idim))/(dx_g(i,idim)+dx_d(i,idim))
           end do
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
            do idust = 1, ndust
               sum_dust=sum_dust+uold(ind_cell(i),firstindex_ndust+idust)/d
            end do
            tstop_tot=0.0d0
            t_stop=0.0d0
            do idust = 1,ndust
               t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/d/(1.0d0-sum_dust)
               if(K_drag)  t_stop(idust) = uold(ind_cell(i),firstindex_ndust+idust)/K_dust(idust)
               if(dust_barr) t_stop (idust)= 0.1_dp
               if (d .le. dens_floor) t_stop(idust) =t_stop_floor 
               tstop_tot= tstop_tot-t_stop(idust)*(uold(ind_cell(i),firstindex_ndust+idust)/d)
            end do
            do idust = 1,ndust
               t_stop(idust) = t_stop(idust)+tstop_tot
               do idim=1,ndim
                  v_dust(ind_cell(i),idust,idim)=t_stop(idust)*gradP(i,idim)/d
               end do   
            end do
         end do
   end do
enddo

111 format('   Entering set_vdust for level ',i2)

end subroutine set_vdust


!###########################################################
!###########################################################
!###########################################################
!###########################################################
   

!subroutine regularize_dust(speedr,speed,dspeed,dx)
!  use amr_parameters
!  use hydro_parameters
!  implicit none
!  real(dp)::speedr,dx
!  real(dp)::speed,dspeed
!  if(visco_dust.eqv..true.) then
!     speedr= tanh(sign(dx,speed)/(eta_dust))*abs(speed)
!  else
!     speedr=speed
!  endif

!end subroutine regularize_dust
