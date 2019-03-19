! ------------------------------------------------------------------------------
!  DUSTDIFF_SPLIT  This routine solves the dust flux after the operator splitting  
!              
!
!  inputs/outputs
!  uin         => (const)  input state
!  iu1,iu2     => (const)  first and last index of input array,
!  ju1,ju2     => (const)  cell centered,    
!  ku1,ku2     => (const)  including buffer cells.
!  flux       <=  (modify) return fluxes in the 3 coord directions
!  if1,if2     => (const)  first and last index of output array,
!  jf1,jf2     => (const)  edge centered,
!  kf1,kf2     => (const)  for active cells only.
!  dx,dy,dz    => (const)  (dx,dy,dz)
!  dt          => (const)  time step
!  ngrid       => (const)  number of sub-grids
!  ndim        => (const)  number of dimensions
!
!  uin = (\rho eps_1, .... \rho eps_ndust,vdust ....,\rho,P)
! --------------------------------------------------------------------------
subroutine dustdiff_predict_eint(uin,flux,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none
  
  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim),save::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim),save::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim),save::qp
  
  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndim)::flux
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::Xflux,Yflux,Zflux

  ! Local scalar variables
  integer::i,j,k,l,ivar, idust
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  !flux=0.0d0
  Xflux=0.0d0
  Yflux=0.0d0
  Zflux=0.0d0
  qm=0.0d0
  qp= 0.0d0
  dq=0.0d0
  ! Compute TVD slopes

  call uslope_dust2(uin,dq,dx,dt,ngrid)
  ! Compute 3D traced-states in all three directions
#if NDIM==1
     call trace1d_dust_eint(uin,dq,qm,qp,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call trace2d_dust_eint(uin,dq,qm,qp,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call trace3d_dust_eint(uin,dq,qm,qp,dx,dy,dz,dt,ngrid)
#endif
  
 ! Solve for 1D flux in X direction
  call cmpflxmdust_eint(qm,iu1+1,iu2+1,ju1  ,ju2  ,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          if1  ,if2  ,jlo  ,jhi  ,klo  ,khi  , 2,3,4,Xflux,ngrid)
  ! Save flux in output array
  do i=if1,if2
  do j=jlo,jhi
  do k=klo,khi
     do idust=1,ndust
        do l=1,ngrid
           flux(l,i,j,k,1)=Xflux(l,i,j,k)*dt/dx
        end do
     end do
  end do
  end do
  end do

  ! Solve for 1D flux in Y direction
#if NDIM>1
  call cmpflxmdust_eint(qm,iu1  ,iu2  ,ju1+1,ju2+1,ku1  ,ku2  , &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jf1  ,jf2  ,klo  ,khi  , 3,2,4,Yflux,ngrid)
 
  ! Save flux in output array
  do i=ilo,ihi
  do j=jf1,jf2
  do k=klo,khi
     do idust=1,ndust
        do l=1,ngrid
           flux(l,i,j,k,2)=Yflux(l,i,j,k)*dt/dy
        end do
     end do
  end do
  end do
  end do
#endif

  ! Solve for 1D flux in Z direction
#if NDIM>2
  call cmpflxmdust_eint(qm,iu1  ,iu2  ,ju1  ,ju2  ,ku1+1,ku2+1, &
       &       qp,iu1  ,iu2  ,ju1  ,ju2  ,ku1  ,ku2  , &
       &          ilo  ,ihi  ,jlo  ,jhi  ,kf1  ,kf2  , 4,2,3,zflux,ngrid)

  ! Save flux in output array
  do i=ilo,ihi
  do j=jlo,jhi
  do k=kf1,kf2
     do idust=1,ndust
        do l=1,ngrid
           flux(l,i,j,k,3)=Zflux(l,i,j,k)*dt/dz
        end do
     end do
  end do
end do
end do
#endif
end subroutine dustdiff_predict_eint

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope_dust2(q,dq,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::dq
  ! local arrays
  integer::i, j, k, l, n,ndvar
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop
#if NDIM==2
  real(dp)::dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr
#endif
#if NDIM==3
  real(dp)::dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl
  real(dp)::dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm
  real(dp)::dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr
  real(dp)::dfz
#endif
#if NDIM>1
  real(dp)::vmin,vmax,dfx,dfy,dff
#endif
  integer::ilo,ihi,jlo,jhi,klo,khi
  ndvar =1 +ndim
  
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_dust==0)then
     dq=zero
     return
  end if
#if NDIM==1
  do n = 1, ndvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_dust==1.or.slope_dust==2.or.slope_dust==3)then  ! minmod or average
                 do l = 1, ngrid
                    dlft = MIN(slope_dust,2)*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = MIN(slope_dust,2)*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/MIN(slope_dust,2)
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else if(slope_dust==4)then ! superbee
                 do l = 1, ngrid
                    dcen = q(l,i,j,k,2)*dt/dx
                    dlft = two/(one+dcen)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                    drgt = two/(one-dcen)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                    dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dsgn = sign(one, dlft)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                 end do
              else if(slope_dust==5)then ! ultrabee
                 if(n==1)then
                    do l = 1, ngrid
                       dcen = q(l,i,j,k,2)*dt/dx
                       if(dcen>=0)then
                          dlft = two/(zero+dcen+1d-10)*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                          drgt = two/(one -dcen      )*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       else
                          dlft = two/(one +dcen      )*(q(l,i,j,k,n)-q(l,i-1,j,k,n))
                          drgt = two/(zero-dcen+1d-10)*(q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       endif
                       dsgn = sign(one, dlft)
                       slop = min(abs(dlft),abs(drgt))
                       dlim = slop
                       dcen = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                       if((dlft*drgt)<=zero)dlim=zero
                       dq(l,i,j,k,n,1) = dsgn*dlim !min(dlim,abs(dcen))
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0.0
                    end do
                 end if
              else if(slope_dust==6)then ! unstable
                 if(n==1)then
                    do l = 1, ngrid
                       dlft = (q(l,i,j,k,n)-q(l,i-1,j,k,n))
                       drgt = (q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       slop = 0.5d0*(dlft+drgt)
                       dlim = slop
                       dq(l,i,j,k,n,1) = dlim
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0.0
                    end do
                 end if
              else if(slope_dust==7)then ! van Leer
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2.0d0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              else if(slope_dust==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)



                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else
                 write(*,*)'Unknown slope type',dx,dt
                 stop
              end if
           end do
        end do
     end do
  end do
#endif

#if NDIM==2
  if(slope_dust==1.or.slope_dust==2)then  ! minmod or average
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_dust*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_dust*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_dust
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_dust*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_dust*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_dust
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_dust==3)then ! positivity preserving 2d unsplit slope
     do n = 1, nvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dfll = q(l,i-1,j-1,k,n)-q(l,i,j,k,n)
                    dflm = q(l,i-1,j  ,k,n)-q(l,i,j,k,n)
                    dflr = q(l,i-1,j+1,k,n)-q(l,i,j,k,n)
                    dfml = q(l,i  ,j-1,k,n)-q(l,i,j,k,n)
                    dfmm = q(l,i  ,j  ,k,n)-q(l,i,j,k,n)
                    dfmr = q(l,i  ,j+1,k,n)-q(l,i,j,k,n)
                    dfrl = q(l,i+1,j-1,k,n)-q(l,i,j,k,n)
                    dfrm = q(l,i+1,j  ,k,n)-q(l,i,j,k,n)
                    dfrr = q(l,i+1,j+1,k,n)-q(l,i,j,k,n)

                    vmin = min(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)
                    vmax = max(dfll,dflm,dflr,dfml,dfmm,dfmr,dfrl,dfrm,dfrr)

                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dff  = half*(abs(dfx)+abs(dfy))

                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif

                    dlim = slop

                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy

                 end do
              end do
           end do
        end do
     end do
  else if(slope_dust==7)then ! van Leer
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2.0*dlft*drgt/(dlft+drgt))
                       end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2)=zero
                    else
                       dq(l,i,j,k,n,2)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_dust==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type',dx,dt
     stop
  endif
#endif

#if NDIM==3
  if(slope_dust==1)then  ! minmod
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i  ,j,k,n) - q(l,i-1,j,k,n)
                    drgt = q(l,i+1,j,k,n) - q(l,i  ,j,k,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,1) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,1) = max(dlft,drgt)
                    end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i,j  ,k,n) - q(l,i,j-1,k,n)
                    drgt = q(l,i,j+1,k,n) - q(l,i,j  ,k,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,2) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,2) = max(dlft,drgt)
                    end if
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = q(l,i,j,k  ,n) - q(l,i,j,k-1,n)
                    drgt = q(l,i,j,k+1,n) - q(l,i,j,k  ,n)
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,3) = zero
                    else if(dlft>0) then
                       dq(l,i,j,k,n,3) = min(dlft,drgt)
                    else
                       dq(l,i,j,k,n,3) = max(dlft,drgt)
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_dust==2)then ! moncen
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_dust*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_dust*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_dust
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_dust*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_dust*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_dust
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = slope_dust*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = slope_dust*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)/slope_dust
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else if(slope_dust==3)then ! positivity preserving 3d unsplit slope
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 do l = 1, ngrid
                    dflll = q(l,i-1,j-1,k-1,n)-q(l,i,j,k,n)
                    dflml = q(l,i-1,j  ,k-1,n)-q(l,i,j,k,n)
                    dflrl = q(l,i-1,j+1,k-1,n)-q(l,i,j,k,n)
                    dfmll = q(l,i  ,j-1,k-1,n)-q(l,i,j,k,n)
                    dfmml = q(l,i  ,j  ,k-1,n)-q(l,i,j,k,n)
                    dfmrl = q(l,i  ,j+1,k-1,n)-q(l,i,j,k,n)
                    dfrll = q(l,i+1,j-1,k-1,n)-q(l,i,j,k,n)
                    dfrml = q(l,i+1,j  ,k-1,n)-q(l,i,j,k,n)
                    dfrrl = q(l,i+1,j+1,k-1,n)-q(l,i,j,k,n)

                    dfllm = q(l,i-1,j-1,k  ,n)-q(l,i,j,k,n)
                    dflmm = q(l,i-1,j  ,k  ,n)-q(l,i,j,k,n)
                    dflrm = q(l,i-1,j+1,k  ,n)-q(l,i,j,k,n)
                    dfmlm = q(l,i  ,j-1,k  ,n)-q(l,i,j,k,n)
                    dfmmm = q(l,i  ,j  ,k  ,n)-q(l,i,j,k,n)
                    dfmrm = q(l,i  ,j+1,k  ,n)-q(l,i,j,k,n)
                    dfrlm = q(l,i+1,j-1,k  ,n)-q(l,i,j,k,n)
                    dfrmm = q(l,i+1,j  ,k  ,n)-q(l,i,j,k,n)
                    dfrrm = q(l,i+1,j+1,k  ,n)-q(l,i,j,k,n)

                    dfllr = q(l,i-1,j-1,k+1,n)-q(l,i,j,k,n)
                    dflmr = q(l,i-1,j  ,k+1,n)-q(l,i,j,k,n)
                    dflrr = q(l,i-1,j+1,k+1,n)-q(l,i,j,k,n)
                    dfmlr = q(l,i  ,j-1,k+1,n)-q(l,i,j,k,n)
                    dfmmr = q(l,i  ,j  ,k+1,n)-q(l,i,j,k,n)
                    dfmrr = q(l,i  ,j+1,k+1,n)-q(l,i,j,k,n)
                    dfrlr = q(l,i+1,j-1,k+1,n)-q(l,i,j,k,n)
                    dfrmr = q(l,i+1,j  ,k+1,n)-q(l,i,j,k,n)
                    dfrrr = q(l,i+1,j+1,k+1,n)-q(l,i,j,k,n)

                    vmin = min(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)
                    vmax = max(dflll,dflml,dflrl,dfmll,dfmml,dfmrl,dfrll,dfrml,dfrrl, &
                         &     dfllm,dflmm,dflrm,dfmlm,dfmmm,dfmrm,dfrlm,dfrmm,dfrrm, &
                         &     dfllr,dflmr,dflrr,dfmlr,dfmmr,dfmrr,dfrlr,dfrmr,dfrrr)

                    dfx  = half*(q(l,i+1,j,k,n)-q(l,i-1,j,k,n))
                    dfy  = half*(q(l,i,j+1,k,n)-q(l,i,j-1,k,n))
                    dfz  = half*(q(l,i,j,k+1,n)-q(l,i,j,k-1,n))
                    dff  = half*(abs(dfx)+abs(dfy)+abs(dfz))

                    if(dff>zero)then
                       slop = min(one,min(abs(vmin),abs(vmax))/dff)
                    else
                       slop = one
                    endif

                    dlim = slop

                    dq(l,i,j,k,n,1) = dlim*dfx
                    dq(l,i,j,k,n,2) = dlim*dfy
                    dq(l,i,j,k,n,3) = dlim*dfz

                 end do
              end do
           end do
        end do
     end do
  else if(slope_dust==7)then ! van Leer
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,2)=zero
                    else
                       dq(l,i,j,k,n,2)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = (q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,3)=zero
                    else
                       dq(l,i,j,k,n,3)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              end do
           end do
        end do
     end do
  else if(slope_dust==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one, dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = (q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = (q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = (q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)
                    dsgn = sign(one,dcen)
                    slop = min(slope_theta*abs(dlft),slope_theta*abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,3) = dsgn*min(dlim,abs(dcen))
                 end do
              end do
           end do
        end do
     end do
  else
     write(*,*)'Unknown slope type',dx,dt
     stop
  endif
#endif

end subroutine uslope_dust2

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmpflxmdust_eint(qm,im1,im2,jm1,jm2,km1,km2, &
     &             qp,ip1,ip2,jp1,jp2,kp1,kp2, &
     &                ilo,ihi,jlo,jhi,klo,khi, ln,lt1,lt2, &
     &            flx,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  integer ::ln,lt1,lt2
  integer ::im1,im2,jm1,jm2,km1,km2
  integer ::ip1,ip2,jp1,jp2,kp1,kp2
  integer ::ilo,ihi,jlo,jhi,klo,khi
  real(dp),dimension(1:nvector,im1:im2,jm1:jm2,km1:km2,1:1+ndim,1:ndim)::qm
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2,1:1+ndim,1:ndim)::qp
  real(dp),dimension(1:nvector,ip1:ip2,jp1:jp2,kp1:kp2)::flx

  ! local variables
  integer ::i, j, k, l, xdim,idust,i0,j0,k0
  real(dp)::entho
  real(dp)::qleft,qright
  real(dp)::uleft,uright


  xdim=ln-1
 
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              qleft = qm(l,i,j,k,1,xdim)
              qright = qp(l,i,j,k,1,xdim)
              ! dust velocity
              uleft = qm(l,i,j,k,1+xdim,xdim)
              uright = qp(l,i,j,k,1+xdim,xdim)
              ! Compute fluxes
              call upwind_dust(uleft,uright,qleft,qright,flx(l,i,j,k))
              
           end do
     end do
  end do
end do
end subroutine cmpflxmdust_eint

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d_dust_eint(q,dq,qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dt
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::qp

  ! Local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::idust
  real(dp)::dtdx
  real(dp)::rhod,ud
  real(dp):: drhox,dux,dux1,dux2
  real(dp)::srho0, su0
  integer :: Ndvar
  Ndvar = 2*ndust*ndim

  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              ! Cell centered values
              rhod =  q(l,i,j,k,1)
              ud   =  q(l,i,j,k,2)
              ! TVD slopes in all 3 directions
              drhox = dq(l,i,j,k,1,1)
              dux= dq(l,i,j,k,2,1)
              ! Source terms (including transverse derivatives)
              srho0 =  -ud*drhox - dux*rhod
              qp(l,i,j,k,1,1) = rhod      + half*srho0*dtdx - half* drhox
              qp(l,i,j,k,2,1) = ud  - half* dux 
              ! Left state at left interface
              qm(l,i,j,k,1,1) = rhod        + half*srho0*dtdx+ half*drhox
              qm(l,i,j,k,2,1) = ud + half*dux
           end do
        end do
     end do
  end do


end subroutine trace1d_dust_eint
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d_dust_eint(q,dq,qm,qp,dx,dy,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::qp
  ! declare local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::idust
  real(dp)::dtdx, dtdy
  real(dp)::rhod,ud,vd
  real(dp):: drhox,drhoy,dux,duy,dvx,dvy
  real(dp):: drhox2,drhoy2,dux2,duy2,dvx2,dvy2
  
  real(dp)::srho0, su0, sv0,srhox,srhoy
  integer :: Ndvar
  Ndvar = 1+ndim

  dtdx = dt/dx
  dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              ! Cell centered values
              rhod =  q(l,i,j,k,1)
              ud   =  q(l,i,j,k,2)
              vd   =  q(l,i,j,k,3)
              ! TVD slopes in all 3 directions
              drhox = dq(l,i,j,k,1,1)
              dux =   dq(l,i,j,k,2,1)
              dvx =   dq(l,i,j,k,3,1)
              drhoy = dq(l,i,j,k,1,2)
              duy =   dq(l,i,j,k,2,2)
              dvy =   dq(l,i,j,k,3,2)

               
              ! Source terms (including transverse derivatives)
              srho0 = -ud*drhox-vd*drhoy- (dux+dvy)*rhod



              ! Right state at left interface
              qp(l,i,j,k,1,1) = rhod       - half*drhox + srho0*dtdx*half
              qp(l,i,j,k,2,1) = ud     - half*dux   !+ su0*dtdx*half
              qp(l,i,j,k,3,1) = vd     - half*dvx  ! + sv0*dtdx*half
              ! Left state at right interface
              qm(l,i,j,k,1,1) = rhod       + half*drhox + srho0*dtdx*half
              qm(l,i,j,k,2,1) = ud     + half*dux   !+ su0*dtdx*half
              qm(l,i,j,k,3,1) = vd     + half*dvx   !+ sv0*dtdx*half
              ! Top state at bottom interface
              qp(l,i,j,k,1,2) = rhod       - half*drhoy + srho0*dtdy*half
              qp(l,i,j,k,2,2) = ud     - half*duy   !+ su0*dtdy*half
              qp(l,i,j,k,3,2) = vd     - half*dvy   !+ sv0*dtdy*half
              ! Bottom state at top interface
              qm(l,i,j,k,1,2) = rhod       + half*drhoy + srho0*dtdy*half
              qm(l,i,j,k,2,2) = ud     + half*duy   !+ su0*dtdy*half
              qm(l,i,j,k,3,2) = vd     + half*dvy   !+ sv0*dtdy*half
           end do
        end do
     end do
  end do


end subroutine trace2d_dust_eint
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d_dust_eint(q,dq,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:1+ndim,1:ndim)::qp

  ! declare local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::idust
  real(dp)::dtdx, dtdy, dtdz
  real(dp)::rhod,ud,vd,wd
  real(dp):: drhox,drhoy,drhoz,dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz
  real(dp):: drhox2,drhoy2,drhoz2,dux2,duy2,duz2,dvx2,dvy2,dvz2,dwx2,dwy2,dwz2
  
  real(dp)::srho0, su0, sv0, sw0, srhox,srhoy,srhoz
  integer :: Ndvar
  Ndvar =1+ndim

  dtdx = dt/dx
  dtdy = dt/dy
  dtdz = dt/dz
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              do idust =1,ndust
              ! Cell centered values
              rhod =  q(l,i,j,k,1)
              ud   =  q(l,i,j,k,2)
              vd   =  q(l,i,j,k,3)
              wd   =  q(l,i,j,k,4)
              ! TVD slopes in all 3 directions
              drhox = dq(l,i,j,k,1,1)
              dux = dq(l,i,j,k,2,1)
              dvx = dq(l,i,j,k,3,1)
              dwx = dq(l,i,j,k,4,1)
              drhoy = dq(l,i,j,k,1,2)
              duy = dq(l,i,j,k,2,2)
              dvy = dq(l,i,j,k,3,2)
              dwy = dq(l,i,j,k,4,2)
              drhoz = dq(l,i,j,k,1,3)
              duz = dq(l,i,j,k,2,3)
              dvz = dq(l,i,j,k,3,3)
              dwz = dq(l,i,j,k,4,3)


              
              ! Source terms (including transverse derivatives)
              srho0 = -ud*drhox-vd*drhoy-wd*drhoz - (dux+dvy+dwz)*rhod


              !su0 = -ud*dux-vd*duy-wd*duz
              !sv0 = -ud*dvx-vd*dvy-wd*dvz 
              !sw0 = -ud*dwx-vd*dwy-wd*dwz 
              ! Right state at left interface
              qp(l,i,j,k,1,1) = rhod       - half*drhox + srho0*dtdx*half
              qp(l,i,j,k,2,1) = ud - half*dux   !+ su0*dtdx*half
              qp(l,i,j,k,3,1) = vd - half*dvx  ! + sv0*dtdx*half
              qp(l,i,j,k,4,1) = wd - half*dwx   !+ sw0*dtdx*half
              ! Left state at left interface
              qm(l,i,j,k,1,1) = rhod       + half*drhox + srho0*dtdx*half
              qm(l,i,j,k,2,1) = ud + half*dux  ! + su0*dtdx*half
              qm(l,i,j,k,3,1) = vd + half*dvx  ! + sv0*dtdx*half
              qm(l,i,j,k,4,1) = wd + half*dwx  ! + sw0*dtdx*half
              ! Top state at bottom interface
              qp(l,i,j,k,1,2) = rhod       - half*drhoy + srho0*dtdy*half
              qp(l,i,j,k,2,2) = ud - half*duy  ! + su0*dtdy*half
              qp(l,i,j,k,3,2) = vd - half*dvy  ! + sv0*dtdy*half
              qp(l,i,j,k,4,2) = wd - half*dwy  ! + sw0*dtdy*half
              ! Bottom state at top interface
              qm(l,i,j,k,1,2) = rhod       + half*drhoy + srho0*dtdy*half
              qm(l,i,j,k,2,2) = ud + half*duy  ! + su0*dtdy*half
              qm(l,i,j,k,3,2) = vd + half*dvy  ! + sv0*dtdy*half
              qm(l,i,j,k,4,2) = wd + half*dwy   !+ sw0*dtdy*half
              ! Back state at front interface
              qp(l,i,j,k,1,3) = rhod       - half*drhoz + srho0*dtdz*half
              qp(l,i,j,k,2,3) = ud - half*duz   !+ su0*dtdz*half
              qp(l,i,j,k,3,3) = vd - half*dvz  ! + sv0*dtdz*half
              qp(l,i,j,k,4,3) = wd - half*dwz  ! + sw0*dtdz*half
              ! Front state at back interface
              qm(l,i,j,k,1,3) = rhod       + half*drhoz + srho0*dtdz*half
              qm(l,i,j,k,2,3) = ud + half*duz  ! + su0*dtdz*half
              qm(l,i,j,k,3,3) = vd + half*dvz  ! + sv0*dtdz*half
              qm(l,i,j,k,4,3) = wd + half*dwz  ! + sw0*dtdz*half
              end do
           end do
        end do
     end do
  end do


end subroutine trace3d_dust_eint
#endif
