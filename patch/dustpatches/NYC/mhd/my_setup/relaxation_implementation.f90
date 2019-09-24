subroutine relaxation(ilevel)
  use amr_commons
  use hydro_commons
  use units_commons
  use poisson_commons
  implicit none

  integer::ilevel
  !--------------------------------------------------------------------------
  ! Enforcing relaxation damping to velocity components at early times,     !
  ! of in geometric regions. This now happens after set_uold, so do to uold !
  !--------------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::ix,iy,iz
  integer::ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,scale,v2,rad2,relax_f,x0,xn,yn
  real(dp)::relax,rin,rr
  real(dp),dimension(ndim)::vel

  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  relax=1000.0

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  rin= 0.5d0
  x0 = boxlen/2.0d0


  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx-x0
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx-x0
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx-x0
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
         RR = sqrt(xn**2.0+yn**2.0)

         if (RR<rin) then
         !   uold(ind_cell(i),2)=0.0d0
         !   uold(ind_cell(i),3)=0.0d0
         !   uold(ind_cell(i),4)=0.0d0
         endif
      enddo

   enddo
    enddo
write (*,*) 'relaxation is on'
111 format('   Entering relaxation_implementation for level ',I2)

end subroutine relaxation
