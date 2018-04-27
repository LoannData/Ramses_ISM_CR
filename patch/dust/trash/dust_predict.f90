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
!  uin = (\rho eps_1, .... \rho eps_ndust, tstop1, tstopndust,\rho,P)

!
!  This routine was adapted from Yohan Dubois & Benoit Commerçon's
!  heat conduction routine by Ugo Lebreuilly
! --------------------------------------------------------------------------
subroutine dustdiff_predict(uin,flux,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none
  
  integer ::ngrid
  real(dp)::dx,dy,dz,dt

  ! Input states
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::uin

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::Xflux,Yflux,Zflux

  ! Local scalar variables
  integer::i,j,k,l,ivar, idust
  integer::ilo,ihi,jlo,jhi,klo,khi

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)

  flux=0.0d0
  Xflux=0.0d0
  Yflux=0.0d0
  Zflux=0.0d0
   ! Translate to primitive variables, compute sound speeds
  call ctoprim_dust(uin,qin,cin,gravin,dt,ngrid)

  ! Compute TVD slopes
  call uslope_dust(qin,dq,dx,dt,ngrid)

  ! Compute 3D traced-states in all three directions
#if NDIM==1
     call trace1d_dust(qin,dq,qm,qp,dx      ,dt,ngrid)
#endif
#if NDIM==2
     call trace2d_dust(qin,dq,qm,qp,dx,dy   ,dt,ngrid)
#endif
#if NDIM==3
     call trace3d_dust(qin,dq,qm,qp,dx,dy,dz,dt,ngrid)
#endif
 
 ! Compute the dust flux in X direction
  call dustXflxP(uin,Xflux,dx,dt,ngrid)

  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        do idust = 1, ndust
           flux(l,i,j,k,idust,1)=Xflux(l,i,j,k,idust)
        end do
    enddo
  enddo
  enddo
  enddo
#if NDIM>1
  ! Compute dust flux in Y direction
  call dustYflxP(uin,Yflux,dy,dt,ngrid)
  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        do idust = 1, ndust         
           flux(l,i,j,k,idust,2)=Yflux(l,i,j,k,idust)
        end do
     enddo
  enddo
  enddo
  enddo
#endif
  
#if NDIM>2
  ! Compute the dust flux in Z direction
  call dustZflxP(uin,Zflux,dz,dt,ngrid)
  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
        do idust = 1, ndust         
          flux(l,i,j,k,idust,3)=Zflux(l,i,j,k,idust)
        end do
     enddo
  enddo
  enddo
  enddo
#endif




end subroutine dustdiff_predict

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

subroutine dustXflxP(uin,myflux,dx,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dx,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::uin

  real(dp)::dx_loc,sum_dust
  real(dp),dimension(1:ndust)::fdust
  real(dp),dimension(1:ndust)::fx
  real(dp) :: speed, ur,ul, sigmas,sigmau,speedr,speedl,dpdx
  integer::i,j,k,l,isl,idust, idens, ipress
  integer::jlo,jhi,klo,khi

  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)
  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  do k=klo,khi
  do j=jlo,jhi
  do i=if1,if2
     do l = 1, ngrid
        do idust=1,ndust
           !First order terms
           speed =  uin(l,i-1,j,k,ndust+idust)
           !left state
           if(slope_dust.eq.1)sigmau=0.0d0
           if(slope_dust.eq.2)call minmod_dust((uin(l,i-1,j,k,idust)-uin(l,i-2,j,k,idust))/dx,(uin(l,i,j,k,idust)-uin(l,i-1,j,k,idust))/dx,sigmau)
           if(slope_dust.eq.3)call vanleer(0.5d0*(uin(l,i,j,k,idust)-uin(l,i-2,j,k,idust))/dx,2.0d0*(uin(l,i-1,j,k,idust)-uin(l,i-2,j,k,idust))/dx,&
                 &2.0d0*(uin(l,i,j,k,idust)-uin(l,i-1,j,k,idust))/dx,sigmau)             
           if(slope_dust.eq.1)sigmas=0.0d0
           if(slope_dust.eq.2)call minmod_dust((uin(l,i-1,j,k,ndust+idust)-uin(l,i-2,j,k,ndust+idust))/dx,(uin(l,i,j,k,ndust+idust)-uin(l,i-1,j,k,ndust+idust))/dx,sigmas)
           if(slope_dust.eq.3)call vanleer(0.5d0*(uin(l,i,j,k,ndust+idust)-uin(l,i-2,j,k,ndust+idust))/dx,2.0d0*(uin(l,i-1,j,k,ndust+idust)-uin(l,i-2,j,k,ndust+idust))/dx,&
                 &2.0d0*(uin(l,i,j,k,ndust+idust)-uin(l,i-1,j,k,ndust+idust))/dx,sigmas)   
           dpdx= (uin(l,i+1,j,k,2*ndust*ndim+2)- uin(l,i-1,j,k,2*ndust*ndim+2))/dx
           call regularize_dust(speedl,speed+0.5d0*sigmas*dx,dpdx)

           ul=uin(l,i-1,j,k,idust)+0.5d0*dt*sigmau*speedl+0.5d0*sigmas*uin(l,i-1,j,k,idust)*dt+0.5d0*sigmau*dx
           !right state
           speed =  uin(l,i,j,k,ndust+idust)
           if(slope_dust.eq.1)sigmau=0.0d0
           if(slope_dust.eq.2)call minmod_dust((uin(l,i,j,k,idust)-uin(l,i-1,j,k,idust))/dx,(uin(l,i+1,j,k,idust)-uin(l,i,j,k,idust))/dx,sigmau)
           if(slope_dust.eq.3)call vanleer(0.5d0*(uin(l,i+1,j,k,idust)-uin(l,i-1,j,k,idust))/dx,2.0d0*(uin(l,i,j,k,idust)-uin(l,i-1,j,k,idust))/dx,&
                 &2.0d0*(uin(l,i+1,j,k,idust)-uin(l,i,j,k,idust))/dx,sigmau)             
           if(slope_dust.eq.1)sigmas=0.0d0
           if(slope_dust.eq.2)call minmod_dust((uin(l,i,j,k,ndust+idust)-uin(l,i-1,j,k,ndust+idust))/dx,(uin(l,i+1,j,k,ndust+idust)-uin(l,i,j,k,ndust+idust))/dx,sigmas)
           if(slope_dust.eq.3)call vanleer(0.5d0*(uin(l,i+1,j,k,ndust+idust)-uin(l,i-1,j,k,ndust+idust))/dx,2.0d0*(uin(l,i,j,k,ndust+idust)-uin(l,i-1,j,k,ndust+idust))/dx,&
                 &2.0d0*(uin(l,i+1,j,k,ndust+idust)-uin(l,i,j,k,ndust+idust))/dx,sigmas)   
           ur=uin(l,i,j,k,idust)+0.5d0*dt*sigmau*speedr+0.5d0*sigmas*uin(l,i,j,k,idust)*dt-0.5d0*sigmau*dx
           call regularize_dust(speedr,speed-0.5d0*sigmas*dx,dpdx)

           call hlldust(speedr,speedl,ul,ur,fx(idust))
           !fx(idust) = fx(idust) + 0.5d0*abs(speed)*(dx-abs(speed)*dt)*sigma
        end do
        do idust= 1, ndust
           myflux(l,i,j,k,idust)=fx(idust)*dt/dx
        end do
    enddo
  enddo
  enddo
  enddo
  
end subroutine dustXflxP

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

subroutine dustYflxP(uin,myflux,dy,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dy,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::uin

  real(dp)::dPdyl,dPdyr,dy_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fy
  !Slopes and advection velocity
  real(dp) :: speed,dspeed,speedtempr,speedtempl, sigma,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp):: t_dyn,entho,pi
  ! Local scalar variables
  integer::i,j,k,l,ivar, idust, idens, ipress,isl
  integer::ilo,ihi,klo,khi
  entho= one /(gamma -one)
  pi =3.14159265358979323846_dp

  klo=MIN(1,ku1+2); khi=MAX(1,ku2-2)
  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  do k=klo,khi
  do j=jf1,jf2
  do i=ilo,ihi
     do l = 1, ngrid
        do idust=1,ndust
  
        end do
        do idust= 1, ndust
           myflux(l,i,j,k,idust)=fy(idust)*dt/dy
        end do  
    enddo
  enddo
  enddo
  enddo

end subroutine dustYflxP

!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################

subroutine dustZflxP(uin,myflux,dz,dt,ngrid)
  use amr_parameters
  use const             
  use hydro_parameters
  implicit none 

  integer ::ngrid
  real(dp)::dz,dt

  ! Output fluxes
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:ndust)::myflux

  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::uin
  real(dp)::dPdzl,dPdzr,dz_loc,sum_dust,Tksleft_tot,Tksright_tot
  real(dp),dimension(1:ndust)::fdust, Tksleft, Tksright
  real(dp),dimension(1:ndust)::fz

  !Slopes and advection velocity
  real(dp) :: speed,speedtempr,speedtempl,dspeed, sigma ,scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp):: t_dyn ,entho, pi
  ! Local scalar variables
  integer::i,j,k,l,ivar, idust, idens, ipress, isl
  integer::ilo,ihi,jlo,jhi,klo,khi
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  entho= one /(gamma -one)
  pi =3.14159265358979323846_dp

  ilo=MIN(1,iu1+2); ihi=MAX(1,iu2-2)
  jlo=MIN(1,ju1+2); jhi=MAX(1,ju2-2)

  do k=kf1,kf2
  do j=jlo,jhi
  do i=ilo,ihi
     do l = 1, ngrid
        do idust=1,ndust

        end do
        do idust= 1, ndust
           myflux(l,i,j,k,idust)=fz(idust)*dt/dz
        end do  
    enddo
  enddo
  enddo
  enddo

end subroutine dustZflxP






!###########################################################
!###########################################################
!###########################################################
!###########################################################
SUBROUTINE hlldust(vleft,vright,uleft,uright,fgdnv)
  USE amr_parameters
  USE const
  ! 1D HLL Riemann solver
  IMPLICIT NONE
  REAL(dp)::fleft,fright,fgdnv
  REAL(dp)::uleft,uright
  REAL(dp):: vleft,vright,bx_mean,SL,SR


  SL=min(min(vleft,vright),zero)
  SR=max(max(vleft,vright),zero)
  fleft= uleft*vleft
  fright= uright*vright
 
  ! the HLL flux
  fgdnv = (SR*fleft-SL*fright+SR*SL*(uright-uleft))/(SR-SL)

END SUBROUTINE hlldust

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine trace1d_dust(q,dq,qm,qp,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dt
  
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::qp

  ! Local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, ip
  real(dp)::dtdx
  real(dp)::rhod,ud
  real(dp):: drhodx,drhody,dux,duy
  real(dp)::srho0, su0
  integer :: Ndvar
  Ndvar = 2*ndust*ndim+2

  dtdx = dt/dx

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; ip=3

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              ! Cell centered values
              rhod =  q(l,i,j,k,idust)
              ud   =  q(l,i,j,k,ndust +idust)
              ! TVD slopes in all 3 directions
              drhodx = dq(l,i,j,k,idust,1)
              dux = dq(l,i,j,k,iu,ndust +idust,1)
              ! Source terms (including transverse derivatives)
              srho0 = -ud*drhox-vd*drhoy- (dux)*rhod
              su0 = -ud*dux
              ! Right state at left interface
              qp(l,i,j,k,idust,1) = rhod       - half*drhox + srho0*dtdx*half
              qp(l,i,j,k,ndust+idust,1) = ud   - half*dux   + su0*dtdx*half
              ! Left state at left interface
              qm(l,i,j,k,idust,1) = rhod       + half*drhox + srho0*dtdx*half
              qm(l,i,j,k,ndust+idust,1) = ud   + half*dux   + su0*dtdx*half

           end do
        end do
     end do
  end do


end subroutine trace1d_dust
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>1
subroutine trace2d_dust(q,dq,qm,qp,dx,dy,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::qp
  ! declare local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, ip
  real(dp)::dtdx, dtdy
  real(dp)::rhod,ud,vd
  real(dp):: drhodx,drhody,dux,duy,dvx,dvy
  real(dp)::srho0, su0, sv0
  integer :: Ndvar
  Ndvar = 2*ndust*ndim+2

  dtdx = dt/dx
  dtdy = dt/dy
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  ir=1; iu=2; iv=3; ip=4

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              do idust=1,ndust
              ! Cell centered values
              rhod =  q(l,i,j,k,idust)
              ud   =  q(l,i,j,k,ndust +idust)
              vd   =  q(l,i,j,k,ndust +idust*2)
              ! TVD slopes in all 3 directions
              drhodx = dq(l,i,j,k,idust,1)
              dux = dq(l,i,j,k,iu,ndust +idust,1)
              dvx = dq(l,i,j,k,iv,ndust+2*idust,1)
              drhody = dq(l,i,j,k,idust,2)
              duy = dq(l,i,j,k,iu,ndust +idust,2)
              dvy = dq(l,i,j,k,iv,ndust+2*idust,2)
              ! Source terms (including transverse derivatives)
              srho0 = -ud*drhox-vd*drhoy- (dux+dvy)*rhod
              su0 = -ud*dux-vd*duy
              sv0 = -ud*dvx-vd*dvy
              ! Right state at left interface
              qp(l,i,j,k,idust,1) = rhod       - half*drhox + srho0*dtdx*half
              qp(l,i,j,k,ndust+idust,1) = ud   - half*dux   + su0*dtdx*half
              qp(l,i,j,k,ndust+idust*2,1) = vd - half*dvx   + sv0*dtdx*half
              ! Left state at left interface
              qm(l,i,j,k,idust,1) = rhod       + half*drhox + srho0*dtdx*half
              qm(l,i,j,k,ndust+idust,1) = ud   + half*dux   + su0*dtdx*half
              qm(l,i,j,k,ndust+idust*2,1) = vd + half*dvx   + sv0*dtdx*half
              ! Top state at bottom interface
              qp(l,i,j,k,idust,2) = rhod       - half*drhoy + srho0*dtdy*half
              qp(l,i,j,k,ndust+idust,2) = ud   - half*duy   + su0*dtdy*half
              qp(l,i,j,k,ndust+idust*2,2) = vd - half*dvy   + sv0*dtdy*half
              ! Bottom state at top interface
              qm(l,i,j,k,idust,2) = rhod       + half*drhoy + srho0*dtdy*half
              qm(l,i,j,k,ndust+idust,2) = ud   + half*duy   + su0*dtdy*half
              qm(l,i,j,k,ndust+idust*2,2) = vd + half*dvy   + sv0*dtdy*half
              end do
           end do
        end do
     end do
  end do


end subroutine trace2d_dust
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
#if NDIM>2
subroutine trace3d_dust(q,dq,qm,qp,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::dq
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::qm
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::qp

  ! declare local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::ir, iu, iv, iw, ip
  real(dp)::dtdx, dtdy, dtdz
  real(dp)::rhod,ud,vd,wd
  real(dp):: drhodx,drhody,drhodz,dux,duy,duz,dvx,dvy,dvz,dwx,dwy,dwz
  real(dp)::srho0, su0, sv0, sw0
  integer :: Ndvar
  Ndvar = 2*ndust*ndim+2

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
              rhod =  q(l,i,j,k,idust)
              ud   =  q(l,i,j,k,ndust +idust)
              vd   =  q(l,i,j,k,ndust +idust*2)
              wd   =  q(l,i,j,k,ndust +idust*3)
              ! TVD slopes in all 3 directions
              drhodx = dq(l,i,j,k,idust,1)
              dux = dq(l,i,j,k,iu,ndust +idust,1)
              dvx = dq(l,i,j,k,iv,ndust+2*idust,1)
              dwx = dq(l,i,j,k,iw,ndust +3*idust,1)
              drhody = dq(l,i,j,k,idust,2)
              duy = dq(l,i,j,k,iu,ndust +idust,2)
              dvy = dq(l,i,j,k,iv,ndust+2*idust,2)
              dwy = dq(l,i,j,k,iw,ndust +3*idust,2)
              drhodz = dq(l,i,j,k,idust,3)
              duz = dq(l,i,j,k,iu,ndust +idust,3)
              dvz = dq(l,i,j,k,iv,ndust+2*idust,3
              dwz = dq(l,i,j,k,iw,ndust +3*idust,3)
              ! Source terms (including transverse derivatives)
              srho0 = -ud*drhox-vd*drhoy-wd*drhoz - (dux+dvy+dwz)*rhod
              su0 = -ud*dux-vd*duy-wd*duz
              sv0 = -ud*dvx-vd*dvy-wd*dvz 
              sw0 = -ud*dwx-vd*dwy-wd*dwz 
              ! Right state at left interface
              qp(l,i,j,k,idust,1) = rhod       - half*drhox + srho0*dtdx*half
              qp(l,i,j,k,ndust+idust,1) = ud   - half*dux   + su0*dtdx*half
              qp(l,i,j,k,ndust+idust*2,1) = vd - half*dvx   + sv0*dtdx*half
              qp(l,i,j,k,ndust+idust*3,1) = wd - half*dwx   + sw0*dtdx*half
              ! Left state at left interface
              qm(l,i,j,k,idust,1) = rhod       + half*drhox + srho0*dtdx*half
              qm(l,i,j,k,ndust+idust,1) = ud   + half*dux   + su0*dtdx*half
              qm(l,i,j,k,ndust+idust*2,1) = vd + half*dvx   + sv0*dtdx*half
              qm(l,i,j,k,ndust+idust*3,1) = wd + half*dwx   + sw0*dtdx*half
              ! Top state at bottom interface
              qp(l,i,j,k,idust,2) = rhod       - half*drhoy + srho0*dtdy*half
              qp(l,i,j,k,ndust+idust,2) = ud   - half*duy   + su0*dtdy*half
              qp(l,i,j,k,ndust+idust*2,2) = vd - half*dvy   + sv0*dtdy*half
              qp(l,i,j,k,ndust+idust*3,2) = wd - half*dwy   + sw0*dtdy*half
              ! Bottom state at top interface
              qm(l,i,j,k,idust,2) = rhod       + half*drhoy + srho0*dtdy*half
              qm(l,i,j,k,ndust+idust,2) = ud   + half*duy   + su0*dtdy*half
              qm(l,i,j,k,ndust+idust*2,2) = vd + half*dvy   + sv0*dtdy*half
              qm(l,i,j,k,ndust+idust*3,2) = wd + half*dwy   + sw0*dtdy*half
              ! Back state at front interface
              qp(l,i,j,k,idust,3) = rhod       - half*drhoz + srho0*dtdz*half
              qp(l,i,j,k,ndust+idust,3) = ud   - half*duz   + su0*dtdz*half
              qp(l,i,j,k,ndust+idust*2,3) = vd - half*dvz   + sv0*dtdz*half
              qp(l,i,j,k,ndust+idust*3,3) = wd - half*dwz   + sw0*dtdz*half
              ! Front state at back interface
              qm(l,i,j,k,idust,3) = rhod       + half*drhoz + srho0*dtdz*half
              qm(l,i,j,k,ndust+idust,3) = ud   + half*duz   + su0*dtdz*half
              qm(l,i,j,k,ndust+idust*2,3) = vd + half*dvz   + sv0*dtdz*half
              qm(l,i,j,k,ndust+idust*3,3) = wd + half*dwz   + sw0*dtdz*half
              end do
           end do
        end do
     end do
  end do


end subroutine trace3d_dust
#endif

!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine uslope_dust(q,dq,dx,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  implicit none

  integer::ngrid
  real(dp)::dx,dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2)::q
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2*ndust*ndim+2,1:ndim)::dq
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
  ndvar = 2*ndust*ndim+2

  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)

  if(slope_type==0)then
     dq=zero
     return
  end if

#if NDIM==1
  do n = 1, ndvar
     do k = klo, khi
        do j = jlo, jhi
           do i = ilo, ihi
              if(slope_type==1.or.slope_type==2.or.slope_type==3)then  ! minmod or average
                 do l = 1, ngrid
                    dlft = MIN(slope_type,2)*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = MIN(slope_type,2)*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/MIN(slope_type,2)
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
              else if(slope_type==4)then ! superbee
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
              else if(slope_type==5)then ! ultrabee
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
              else if(slope_type==6)then ! unstable
                 if(n==1)then
                    do l = 1, ngrid
                       dlft = (q(l,i,j,k,n)-q(l,i-1,j,k,n))
                       drgt = (q(l,i+1,j,k,n)-q(l,i,j,k,n))
                       slop = 0.5*(dlft+drgt)
                       dlim = slop
                       dq(l,i,j,k,n,1) = dlim
                    end do
                 else
                    do l = 1, ngrid
                       dq(l,i,j,k,n,1) = 0.0
                    end do
                 end if
              else if(slope_type==7)then ! van Leer
                 do l = 1, ngrid
                    dlft = (q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = (q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    if((dlft*drgt)<=zero) then
                       dq(l,i,j,k,n,1)=zero
                    else
                       dq(l,i,j,k,n,1)=(2.0*dlft*drgt/(dlft+drgt))
                    end if
                 end do
              else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
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
  if(slope_type==1.or.slope_type==2)then  ! minmod or average
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
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
  else if(slope_type==3)then ! positivity preserving 2d unsplit slope
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
  else if(slope_type==7)then ! van Leer
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
  else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
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
  if(slope_type==1)then  ! minmod
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
  else if(slope_type==2)then ! moncen
     do n = 1, ndvar
        do k = klo, khi
           do j = jlo, jhi
              do i = ilo, ihi
                 ! slopes in first coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i  ,j,k,n) - q(l,i-1,j,k,n))
                    drgt = slope_type*(q(l,i+1,j,k,n) - q(l,i  ,j,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one, dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,1) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in second coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j  ,k,n) - q(l,i,j-1,k,n))
                    drgt = slope_type*(q(l,i,j+1,k,n) - q(l,i,j  ,k,n))
                    dcen = half*(dlft+drgt)/slope_type
                    dsgn = sign(one,dcen)
                    slop = min(abs(dlft),abs(drgt))
                    dlim = slop
                    if((dlft*drgt)<=zero)dlim=zero
                    dq(l,i,j,k,n,2) = dsgn*min(dlim,abs(dcen))
                 end do
                 ! slopes in third coordinate direction
                 do l = 1, ngrid
                    dlft = slope_type*(q(l,i,j,k  ,n) - q(l,i,j,k-1,n))
                    drgt = slope_type*(q(l,i,j,k+1,n) - q(l,i,j,k  ,n))
                    dcen = half*(dlft+drgt)/slope_type
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
  else if(slope_type==3)then ! positivity preserving 3d unsplit slope
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
  else if(slope_type==7)then ! van Leer
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
  else if(slope_type==8)then ! generalized moncen/minmod parameterisation (van Leer 1979)
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

end subroutine uslope_dust
