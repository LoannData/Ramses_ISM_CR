subroutine cmpdtdust(uu,gg,dx,vdu,dt,ncell)
  use amr_parameters
  use hydro_parameters
  use const
  use hydro_commons
  implicit none
  integer::ncell,ht
  real(dp)::dx,dt
  real(dp),dimension(1:nvector,1:nvar+3)::uu
  real(dp),dimension(1:nvector,1:ndim)::gg
  real(dp),dimension(1:nvector)::vdu

  real(dp),dimension(1:nvector),save::a2,B2,rho,ctot
  real(dp)::dtcell,smallp,cf,cc,bc,bn
  integer::k,idim
  real(dp)::sum_dust,dt_dust 
  integer::idust
#if NENER>0
  integer::irad
#endif
  smallp = smallr*smallc**2/gamma
  print *, smallc
  ! Convert to primitive variables
  do k = 1,ncell
     uu(k,1)=max(uu(k,1),smallr)
     rho(k)=uu(k,1)
    
  end do
  do idim = 1,3
     do k = 1, ncell
        uu(k,idim+1) = uu(k,idim+1)/rho(k)
     end do
  end do

  do k = 1,ncell
     B2(k)=zero
  end do
  do idim = 1,3
     do k = 1, ncell
        Bc = half*(uu(k,5+idim)+uu(k,nvar+idim))
        B2(k)=B2(k)+Bc**2
        uu(k,5) = uu(k,5)-half*uu(k,1)*uu(k,idim+1)**2-half*Bc**2
     end do
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        uu(k,5) = uu(k,5)-uu(k,8+irad)
     end do
  end do
#endif

  ! Compute thermal sound speed
     do k = 1, ncell
     sum_dust= 0.0d0
#if NDUST>0
     do idust=1,ndust
        sum_dust= sum_dust + uu(k,firstindex_ndust+idust)/rho(k)
     enddo
#endif
     
     uu(k,5) = max((gamma-one)*uu(k,5),smallp)
     a2(k)=gamma*uu(k,5)/uu(k,1)/(1.0d0-sum_dust)
  end do
#if NENER>0
  do irad = 1,nener
     do k = 1, ncell
        a2(k) = a2(k) + gamma_rad(irad)*(gamma_rad(irad)-1.0d0)*uu(k,8+irad)/uu(k,1)
     end do
  end do
#endif

  ! Compute maximum wave speed (fast magnetosonic)
  do k = 1, ncell
     ctot(k)=zero
  end do
  if(ischeme.eq.1)then
     do idim = 1,ndim   ! WARNING: ndim instead of 3
        do k = 1, ncell
           ctot(k)=ctot(k)+abs(uu(k,idim+1))
#if NDUST>0           
           ctot(k)=ctot(k)+vdu(k)
#endif              
        end do
     end do
  else
     do idim = 1,ndim   ! WARNING: ndim instead of 3
        do k = 1, ncell
           cc=half*(B2(k)/rho(k)+a2(k))
           BN=half*(uu(k,5+idim)+uu(k,nvar+idim))
           cf=sqrt(cc+sqrt(cc**2-a2(k)*BN**2/rho(k)))
           ctot(k)=ctot(k)+abs(uu(k,idim+1))+cf
#if NDUST>0
           ctot(k)=ctot(k)+vdu(k)
#endif            
        end do
     end do
  endif

  ! Compute gravity strength ratio
  do k = 1, ncell
     rho(k)=zero
  end do
  do idim = 1,ndim
     do k = 1, ncell
        rho(k)=rho(k)+abs(gg(k,idim))
     end do
  end do
  do k = 1, ncell
     rho(k)=rho(k)*dx/ctot(k)**2
     rho(k)=MAX(rho(k),0.0001_dp)
  end do

  ! Compute maximum time step for each authorized cell
  dt = courant_factor*dx/smallc
  do k = 1,ncell
           dtcell=dx/ctot(k)*(sqrt(one+two*courant_factor*rho(k))-one)/rho(k)
           dt = min(dt,dtcell,dx/vdu(k))
  end do

end subroutine cmpdtdust
