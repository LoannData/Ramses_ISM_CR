subroutine interpol_hydro_dust(u1,ind1,u2,nn)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:ndust*ndim)::u1
  integer ,dimension(1:nvector,0:twondim)           ::ind1
  real(dp),dimension(1:nvector,1:twotondim,1:ndust*ndim)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are:
  ! the dust velocities
  real(dp)::oneover_twotondim
  integer::i,j,ivar,jdim,idust,idim,ind,ix,iy,iz
#if NENER>0
  integer::irad
#endif
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim),save::a
  real(dp),dimension(1:nvector,1:ndim),save::w

  ! volume fraction of a fine cell relative to a coarse cell
  oneover_twotondim=1.D0/dble(twotondim)

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  !------------------------------------------------
  ! Loop over cell-centered interpolation variables
  !------------------------------------------------
  do idim=1,ndim
     do idust=1,ndust
     ! Load father variable
     do j=0,twondim
        do i=1,nn
           a(i,j)=u1(i,j,ndim*(idust-1)+idim)
         !  print *,a(i,j)
        end do
     end do
     ! Reset gradient
     w(1:nn,1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type==1)call compute_limiter_minmod(a,w,nn)
     if(interpol_type==2)call compute_limiter_central(a,w,nn)
     if(interpol_type==3)call compute_central(a,w,nn)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(1:nn,ind,ndim*(idust-1)+idim)=a(1:nn,0)
        do jdim=1,ndim
           do i=1,nn
              u2(i,ind,ndim*(idust-1)+idim)=u2(i,ind,ndim*(idust-1)+idim)+w(i,jdim)*xc(ind,jdim)
           end do
        end do
     end do

end do
enddo
  ! End loop over cell-centered variables

end subroutine interpol_hydro_dust
