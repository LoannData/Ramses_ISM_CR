subroutine relaxation_implementation(ilevel,nstp)
  use amr_commons
  use hydro_commons
  use units_commons
  use poisson_commons
  use cooling_module      , only : kb,mh
  use radiation_parameters,only:mu_gas
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
    real(dp):: sinthetai, sintheta,alpha_disk,k_corona,v_kep,cs_iso,n_disk,buffer_H

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
    if(hayashi) then
     rhayash=1.0d0
  endif
#if NDUST>0
  epsilon_0 = dust_ratio(1)
#endif
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  relaxinit=1000.0
  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  x0 = boxlen/2.0d0
  
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  rout= 5.0d0!2.0*4.0
  rin =rd_factor
  rho0=rhocen/scale_d
  emass=rin! Softening length
  H=HoverR*rout
  Cs0 =  sqrt(gamma*kb*Tp0/(mu_gas*mH))/scale_v
  trel= trelax/scale_t*3.154d7
    csback= sqrt(gamma*kb*Tpback/(mu_gas*mH))/scale_v
  csmax=sqrt(gamma*kb*4000/(mu_gas*mH))/scale_v
  ! Set position of cell centers relative to grid center
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
         RR = sqrt(xn**2.0+yn**2.0+rsmooth**2.0)
         RRR = sqrt(xn**2.0+yn**2.0)
         RRdR= sqrt(xn**2.0+yn**2.0+rsmooth**2.0)+dx
         radius = sqrt(RR**2.0+zn**2.0)
         sintheta=zn/sqrt(zn**2.0+RR**2.0)
         sinthetai= 1.0 + log10(10e-3)*hoverr**2.
     
         alpha_disk=-1.
         n_disk=-2.
         k_corona= 6.
         radiusdr= sqrt(RRdr**2.0+zn**2.0)

 
          if(iso_smooth) then
             if (rrr<rin.or.abs(zn)>Hsmooth) then
                if(Gressel)cs= cs0/sqrt(RR/rout)*sfive(rrr/rsmooth)+csback
                if(Hayashi)THayash= 280.0d0*(rr/rhayash)**(-1./2.)
                if(Hayashi) cs =  sfive(rr/rsmooth)*sqrt(gamma*kb*THayash/(mu_gas*mH))/scale_v+csback
                H=hoverr*RRR
                 buffer_H=0.1
                 if(bethune) then
                      alpha_disk=-1.
                      n_disk=-2.
                      k_corona= 6.
                      THayash= 300.0*(rr/rhayash)**(alpha_disk/2.0)
                      cs_iso=sqrt(gamma*kb*THayash/(mu_gas*mH))/scale_v
                 if(abs(zn).lt.3.72*H-H*buffer_H) then
                    cs = cs_iso
                 else
                    if (abs(zn).gt.3.72*H) then
                       cs=k_corona*cs_iso
                    else
                       cs= cs_iso*(k_corona+(1.0-k_corona)*(abs(zn)-3.72*H+H*buffer_H)/(H*buffer_H))
                    endif
                 endif
               endif
                d=max(uold(ind_cell(i),1),smallr)
                u=uold(ind_cell(i),2)/d
                v=uold(ind_cell(i),3)/d
                w=uold(ind_cell(i),4)/d
                A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
                B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
                C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
                e=0.5d0*d*(u**2+v**2+w**2)+0.5d0*(A**2+B**2+C**2)
                sum_dust=0.0d0
#if NDUST>0
                do idust =1,ndust
                   sum_dust = sum_dust + uold(ind_cell(i), firstindex_ndust+idust)/d
                end do
#endif 
                uold(ind_cell(i),5)=e+(1.0d0-sum_dust)*d*cs**2.0/(gamma-1.0d0)
             endif
                if(Gressel)cs= cs0/sqrt(RR/rout)*sfive(rrr/rsmooth)+csback
                if(Hayashi)THayash= 280.0d0*(rr/rhayash)**(-1./2.)
                if(Hayashi) cs =  sfive(rr/rsmooth)*sqrt(gamma*kb*THayash/(mu_gas*mH))/scale_v+csback
                 H=hoverr*RRR
                 buffer_H=0.1
                 if(bethune) then
                      alpha_disk=-1.
                      n_disk=-2.
                      k_corona= 6.
                      THayash= 300.0*(rr/rhayash)**(alpha_disk/2.0)
                      cs_iso=sqrt(gamma*kb*THayash/(mu_gas*mH))/scale_v                    
                 if(abs(zn).lt.3.72*H-H*buffer_H) then
                    cs = cs_iso
                 else
                    if (abs(zn).gt.3.72*H) then
                       cs=k_corona*cs_iso
                    else
                       cs= cs_iso*(k_corona+(1.0-k_corona)*(abs(zn)-3.72*H+H*buffer_H)/(H*buffer_H))
                    endif
                 endif
               endif  
                 d=max(uold(ind_cell(i),1),smallr)
             smoothing=sfive(abs(0.5*boxlen-Hsmooth)/abs(zn))*sfive(rrr/rsmooth)
             if(damp)uold(ind_cell(i),2)=uold(ind_cell(i),2)/abs(uold(ind_cell(i),2))*min(abs(uold(ind_cell(i),2)),1e8/scale_v)
             if(damp)uold(ind_cell(i),3)=uold(ind_cell(i),3)/abs(uold(ind_cell(i),3))*min(abs(uold(ind_cell(i),3)),1e8/scale_v)
             if(damp)uold(ind_cell(i),4)=uold(ind_cell(i),4)/abs(uold(ind_cell(i),4))*min(abs(uold(ind_cell(i),4)),1e8/scale_v)
             
             u=uold(ind_cell(i),2)/d
             v=uold(ind_cell(i),3)/d
             w=uold(ind_cell(i),4)/d
             A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
             B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
             C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
             e=0.5d0*d*(u**2+v**2+w**2)+0.5d0*(A**2+B**2+C**2)
             cloc= sqrt(gamma*(uold(ind_cell(i),5)-e)*(gamma-1.0)/d)
             if(cloc>csmax) uold(ind_cell(i),5)=e+d*cs**2.0/(gamma-1.0d0)
             
          endif
         
  
      enddo

   enddo
enddo
111 format('   Entering relaxation_implementation for level ',I2)

end subroutine relaxation_implementation
function Sone(x)
  use amr_commons
  implicit none
  real(dp) :: x
  real(dp) :: sone,sf
  
  sf = -2*x**3+3*x**2
  if(x.le.0.0d0) then
     sone= 0.0d0
  else if(x<1.0d0)then
     sone= sf
  else if (x.ge.1.0d0) then
     sone= 1.0d0
  endif
  
end function Sone


function Sfive(r)
  use amr_commons
  implicit none
  real(dp) :: r
  real(dp) :: sfive,sf
  
  sf = r**6.0*(462.0d0 -1980.0d0*r+3465.0d0*r**2.0-3080.0d0*r**3.0+1386.0d0*r**4.0-252.0d0*r**5.0)
  
  if(r.le.0.0d0) then
     sfive= 0.0d0
  else if(r<1.0d0)then
     sfive= sf
  else if (r.ge.1.0d0) then
     sfive= 1.0d0
  endif
  
end function Sfive


function Ssix(x)
  use amr_commons
  implicit none
  real(dp) :: x
  real(dp) :: ssix,sf
  
  sf = 924.*x**13-6006.*x**12+16380.*x**11-24024.*x**10+20020.*x**9-9009.*x**8+1716.*x**7
  if(x.le.0.0d0) then
     ssix= 0.0d0
  else if(x<1.0d0)then
     ssix= sf
  else if (x.ge.1.0d0) then
     ssix= 1.0d0
  endif
  
end function Ssix

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold_relax(ilevel)
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
           do idust=1,ndust
!!$              uold(active(ilevel)%igrid(i)+iskip,1) = unew(active(ilevel)%igrid(i)+iskip,1)
!!$              uold(active(ilevel)%igrid(i)+iskip,2) = unew(active(ilevel)%igrid(i)+iskip,2)
!!$              uold(active(ilevel)%igrid(i)+iskip,3) = unew(active(ilevel)%igrid(i)+iskip,3)
!!$              uold(active(ilevel)%igrid(i)+iskip,4) = unew(active(ilevel)%igrid(i)+iskip,4)
!!$              uold(active(ilevel)%igrid(i)+iskip,5) = unew(active(ilevel)%igrid(i)+iskip,5)

           end do
        end do
     end do

111 format('   Entering set_uold_dust for level ',i2)

end subroutine set_uold_relax
!###########################################################
!###########################################################
!###########################################################
!###########################################################
