subroutine relaxation_implementation(ilevel,nstp)
  use amr_commons
  use hydro_commons
  use units_commons
  use poisson_commons
   use cloud_module

  use cooling_module      , only : kb,mh
  use radiation_parameters
  implicit none

  integer::ilevel,nstp
  !--------------------------------------------------------------------------
  ! Enforcing relaxation damping to velocity components at early times,     !
  ! of in geometric regions. This now happens after set_uold, so do to uold !
  !--------------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::ix,iy,iz,id1,ig1,ih1,id2,ig2,ih2
  integer::ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp) ::A,B,C,e_mag,e_kin1,e_kin2,enint,d,u,v,w,e
  real(dp)::dx,scale,x0,xn,yn,zn, radius,rho0,cs0,cs,omega_kep,radiusin, radiusout,smoothing
  real(dp)::relax,rin,rr,dx_loc,rout,emass,relaxinit,H,csdr,drrho,radiusdr,rrdr,rrr,rrrm,csback,csmax,cloc,hsmooth
  real(dp):: sinthetai, sintheta,v_kep,cs_iso,old_Temp,new_Temp,ht,tauc,r0

  real(dp) :: Teq
  real(dp),dimension(ndim)::vel

  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp):: epsilon_0,sum_dust
  integer::idust

  real(dp) :: sfive,ssix,trel,sone
  real(dp) :: dref, pref,uref,vref,wref
   real(dp):: sigmaHayash, THayash, rHayash,pi
 
#if NDUST>0
  real(dp),dimension(1:ndust):: dustMRN
#endif
    integer,dimension(1:3,1:2,1:8)::iii,jjj

#if NDUST>0
  epsilon_0 = dust_ratio(1)
#endif
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
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  x0 = boxlen/2.0d0
  
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
    ngrid=MIN(nvector,ncache-igrid+1)
    do i=1,ngrid
      ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
    end do

    ! Loop over cells
    do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
         ind_cell(i)=ind_grid(i)+iskip
      end do
          ! Gather cell centre positions
      do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
         end do
      end do
      ! Rescale position from code units to user units
      do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
         end do
      end do
      do i=1,ngrid
         xn=(xx(i,1)-x0)
         yn=(xx(i,2)-x0)
         zn=(xx(i,3)-x0)
         RR= sqrt(xn**2.0+yn**2.0+yn**2.0)
      
      
         if(RR.ge.r0) v_dust(ind_cell(i),:,:) =0.0d0
  
      enddo

   enddo
enddo
111 format('   Entering relaxation_implementation for level ',I2)

end subroutine relaxation_implementation
