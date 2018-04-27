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
       do i=1,ngrid
            do idust = 1,ndust
               v_dust(ind_cell(i),idust,1)=speedx
              IF(NDIM>1) v_dust(ind_cell(i),idust,2)=speedy
              IF (NDIM>2) v_dust(ind_cell(i),idust,3)=speedz
           end do
           end do
         end do
      end do

111 format('   Entering set_vdust for level ',i2)

end subroutine set_vdust


!###########################################################
!###########################################################
!###########################################################
!###########################################################
   

subroutine regularize_dust(speedr,speed,dspeed,dx)
  use amr_parameters
  use hydro_parameters
  implicit none
  real(dp)::speedr,dx
  real(dp)::speed,dspeed
  if(visco_dust.eqv..true.) then
     speedr= tanh(sign(dx,speed)/(eta_dust))*abs(speed)
  else
     speedr=speed
  endif

end subroutine regularize_dust
