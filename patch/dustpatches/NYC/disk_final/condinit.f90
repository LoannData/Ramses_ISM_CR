!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use hydro_commons
  use cooling_module      , only : kb,mh
  use radiation_parameters,only:mu_gas
   use units_commons

  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar, idust, i
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::xn,x0,sum_dust, RR,rrr, rin, rout,yn,zn, radius,rho0,cs0,cs,omega_kep,radiusin, radiusout,emass,H,drRho
  real(dp)::rrdr,radiusdr,csdr,sfive,Bz,bzdr,rrdmr,bzdmr,csdmr,radiusdmr,csback,Hsmooth,delta_rho
  real(dp):: sinthetai, sintheta,alpha_disk,k_corona,v_kep,cs_iso,n_disk,buffer_H,rho_surf
!hayashi's params
  real(dp):: sigmaHayash, THayash, rHayash,pi
  real(dp),dimension(1:ndim) :: vtur
#if NDUST>0
  real(dp),dimension(1:ndust):: dustMRN
#endif
  real(dp):: epsilon_0
#if NDUST>0
  epsilon_0 = dust_ratio(1)
#endif

  do ivar=9,nvar
     q(1:nn,ivar)=0.0d0
  end do
  rin =rd_factor
  rout= 5.0d0!2.0*4.0
  H=HoverR*rout
  Cs0 =  sqrt(gamma*kb*Tp0/(mu_gas*mH))/scale_v
  csback= sqrt(gamma*kb*Tpback/(mu_gas*mH))/scale_v
  pi =3.14159265358979
  rho0=rhocen/scale_d
  x0=boxlen/2.0
  
  delta_rho=0.1
  rhayash=1.0
  if(hayashi) then
     rhayash=1.0d0
  endif
  do i=1,nn
     ! Compute position in normalized coordinates
     xn=(x(i,1)-x0)
     yn=(x(i,2)-x0)
     zn=(x(i,3)-x0)
     !cylindrical radius
     !spherical radius
     RR = sqrt(xn**2.0+yn**2.0+rsmooth**2.0)
     RRR = sqrt(xn**2.0+yn**2.0)
     
     radius = sqrt(RR**2.0+zn**2.0)

     RRdR=RR+dx
     RRdmR=abs(RR-dx)

     radiusdr= sqrt(RRdr**2.0+zn**2.0)
     radiusdmr=sqrt(RRdmr**2.0+zn**2.0)

!Gressel model     
     if (Gressel) then
        
     cs= sfive(rrr/rsmooth)*cs0/sqrt(RR/rout)+csback
     csdr=sfive(rrr/rsmooth)*Cs0/sqrt(RRdr/rout)+csback
     csdmr=sqrt(RRdmr**2.0+zn**2.0)     
     Hsmooth=5.*hoverr*rrr
      sum_dust=0.0d0

      q(i,1)=max(rho0*(RR/rout)**(-2.5)*exp(-Mstar_cen/cs**2.0*(1/RR-1/radius)),rhoext/scale_d)&
           &*(1.0+delta_rho*cos(2.*xn/rrr))
      
      
#if NDUST>0
     if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
        do idust =1,ndust
           q(i, firstindex_ndust+idust)= dust_ratio(idust)/(1.0d0+dust_ratio(idust))
           !if(q(i,1).eq.rhoext/scale_d) q(i, firstindex_ndust+idust)=1.d-18

           if(mrn) q(i, firstindex_ndust+idust) = dustMRN(idust)
           
           sum_dust = sum_dust + q(i, firstindex_ndust+idust)
        end do   
#endif    
        drRho= (csdr**2.0*rho0*(RRdr/rout)**(-2.5)*exp(-Mstar_cen/csdr**2.0*(1/RRdr-1/radiusdr))-cs**2.0*rho0*(RR/rout)**(-2.5)*exp(-Mstar_cen/cs**2.0*(1/RR-1/radius)))/q(i,1)/RR/dx




     omega_kep=sqrt(Mstar_cen/rr**3.0+drRho)!-3.5*Bz*(RR/rout)**(-1.0)/q(i,1)/RR)
     q(i,5)=cs**2.0*q(i,1)
     Bz=sfive(rr/rsmooth)*sqrt(2.0d0*(cs**2.0*rho0*(RR/rout)**(-2.5))/beta_mag)

     call get_vturb(turb_perc,cs,vtur)
     !call get_dturb(turb_perc,q(i,1),delro)
     q(i,2)=-omega_kep*rr*yn/rrr+vtur(1)
     q(i,3)=omega_kep*rr*xn/rrr+vtur(2)
     q(i,4)=vtur(3)

     q(i,6)     = 0.d0
     q(i,7)     = 0.d0
     q(i,8)     = bz
     
     ! Right B fields. Div * B is zero with linear operator (Bxr - Bxl)/dx ...
     q(i,nvar+1)= 0.d0
     q(i,nvar+2)= 0.d0
     q(i,nvar+3)= q(i,8)
  endif

  !MMSN model Hayashi
  if (Hayashi) then
     omega_kep=sqrt(Mstar_cen/rr**3.0)

     sigmaHayash= 1700.d0/scale_m*scale_l**2.0*(rr/rhayash)**(-3./2.)
   
     THayash= 280.0d0*(rr/rhayash)**(-1./2.)
     cs =  sfive(rr/rsmooth)*sqrt(gamma*kb*THayash/(mu_gas*mH))/scale_v+csback
     H=cs/omega_kep
#if NDUST>0
     if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
        do idust =1,ndust
           q(i, firstindex_ndust+idust)= dust_ratio(idust)/(1.0d0+dust_ratio(idust))
           !if(q(i,1).eq.rhoext/scale_d) q(i, firstindex_ndust+idust)=1.d-18

           if(mrn) q(i, firstindex_ndust+idust) = dustMRN(idust)
           
           sum_dust = sum_dust + q(i, firstindex_ndust+idust)
        end do   
#endif
     q(i,1)= max(sigmaHayash/sqrt(2.0*pi*H**2.0)*exp(-zn**2.0/(2.0*H**2.0)),rhoext/scale_d)
     q(i,5)=cs**2.0*q(i,1)
     Bz=sqrt(2.0d0*(cs**2.0*sigmaHayash/sqrt(2.0*pi*H**2.0))/beta_mag)
     call get_vturb(turb_perc,cs,vtur)
     !call get_dturb(turb_perc,q(i,1),delro)
     q(i,2)=-omega_kep*rr*yn/rrr+vtur(1)
     q(i,3)=omega_kep*rr*xn/rrr+vtur(2)
     q(i,4)=vtur(3)

     q(i,6)     = 0.d0
     q(i,7)     = 0.d0
     q(i,8)     = bz
     
     ! Right B fields. Div * B is zero with linear operator (Bxr - Bxl)/dx ...
     q(i,nvar+1)= 0.d0
     q(i,nvar+2)= 0.d0
     q(i,nvar+3)= q(i,8)
  endif

  if(bethune) then
     sintheta=abs(zn)/sqrt(zn**2.0+RRR**2.0)
     sinthetai= 1.0 + log10(1.e-3)*hoverr**2.0
#if NDUST>0
     if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
        do idust =1,ndust
           q(i, firstindex_ndust+idust)= dust_ratio(idust)/(1.0d0+dust_ratio(idust))
           !if(q(i,1).eq.rhoext/scale_d) q(i, firstindex_ndust+idust)=1.d-18

           if(mrn) q(i, firstindex_ndust+idust) = dustMRN(idust)
           
           sum_dust = sum_dust + q(i, firstindex_ndust+idust)
        end do   
#endif
     
     alpha_disk=-1.
     n_disk=-3./2.
     k_corona= 6.
     THayash= 300.0*(rr/rhayash)**(alpha_disk/2.0)
     cs_iso=sfive(rr/rsmooth)*sqrt(gamma*kb*THayash/(mu_gas*mH))/scale_v+csback
     H=hoverr*RRR
     buffer_H=2.5
     cs = cs_iso
     rho_surf=rho0*(RR/rhayash)**n_disk*exp(-Mstar_cen/cs**2.0*(1/RR-1/sqrt(RR**2.0+(3.72*H)**2.0)))
     if(abs(zn).lt.3.72*H) then
     cs = cs_iso
     q(i,1)=rho0*(RR/rhayash)**n_disk*exp(-Mstar_cen/cs**2.0*(1/RR-1/radius))
     else
     if (abs(zn).gt.3.72*H+H*buffer_H) then
     cs=k_corona*cs_iso
     q(i,1)=rho_surf*exp(-Mstar_cen/cs**2.0*(1/RR-1/radius))
     else
     cs= cs_iso*(k_corona-(1.0-k_corona)*(abs(zn)-3.72*H-H*buffer_H)/(H*buffer_H))
     q(i,1)=rho_surf*exp(-Mstar_cen/cs**2.0*(1/RR-1/radius))
     endif
     endif

     omega_kep=sqrt(Mstar_cen/rr**3.0)!-3.5*Bz*(RR/rout)**(-1.0)/q(i,1)/RR)
     call get_vturb(turb_perc,cs,vtur)
     !call get_dturb(turb_perc,q(i,1),delro)
     q(i,2)=-omega_kep*rr*yn/rrr+vtur(1)
     q(i,3)=omega_kep*rr*xn/rrr+vtur(2)
     q(i,5)=q(i,1)*cs**2.0

     Bz=sqrt(2.0d0*(cs_iso**2.0*rho0*(RR/rhayash)**n_disk/beta_mag))
     q(i,6)     = 0.d0
     q(i,7)     = 0.d0
     q(i,8)     = bz
     
     ! Right B fields. Div * B is zero with linear operator (Bxr - Bxl)/dx ...
     q(i,nvar+1)= 0.d0
     q(i,nvar+2)= 0.d0
     q(i,nvar+3)= q(i,8)
  endif

  
     end do

     ! Convert primitive to conservative variables
     ! density -> density
     u(1:nn,1)=q(1:nn,1)
     ! velocity -> momentum
     u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
     u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
     u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
     ! kinetic energy
     u(1:nn,5)=0.0d0
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
     u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
     ! pressure -> total fluid energy
     u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
     ! magnetic energy -> total fluid energy
     u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
     u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
     u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
     u(1:nn,6:8)=q(1:nn,6:8)
     u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
#if NENER>0
     ! non-thermal pressure -> non-thermal energy
     ! non-thermal energy   -> total fluid energy
     do ivar=1,nener-ngrp
        u(1:nn,8+ivar)=q(1:nn,8+ivar)/(gamma_rad(ivar)-1.0d0)
        u(1:nn,5)=u(1:nn,5)+u(1:nn,8+ivar)
     enddo
    ! Radiative transfer
#if NGRP>0
     ! radiative energy   -> total fluid energy
     do ivar=1,ngrp
        u(1:nn,firstindex_er+ivar)= q(1:nn,firstindex_er+ivar)
        u(1:nn,5)=u(1:nn,5)+ u(1:nn,firstindex_er+ivar)
     enddo
#if USE_M_1==1
     ! radiative flux
     do ivar=1,ndim*ngrp
        do i=1,ncache
           u(1:nn,fisrtindex_fr+ivar)=q(1:nn,firstindex+ivar)
        end do
        write(ilun)xdp
     end do
#endif
#endif
#endif
#if NEXTINCT>0
     ! Extinction
     if(extinction)u(1:nn,firstindex_extinct+nextinct)=1.0D0
#endif
#if NPSCAL>0
     ! passive scalars
     do ivar=1,npscal
        u(1:nn,firstindex_pscal+ivar)=q(1:nn,1)*q(1:nn,firstindex_pscal+ivar)
     end do
     ! Internal energy
     u(1:nn,nvar)=q(1:nn,5)/(gamma-1.0d0)
#endif
#if NDUST>0
     ! dust
     do ivar=1,ndust
        u(1:nn,firstindex_ndust+ivar)=q(1:nn,1)*q(1:nn,firstindex_ndust+ivar)
     end do
#endif

   end subroutine condinit

   !================================================================
!================================================================
!================================================================
!================================================================
subroutine columninit(x,sigm,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use cooling_module      , only : kb,mh
  use radiation_parameters,only:mu_gas
   use units_commons

  implicit none
  integer ::nn                              ! Number of cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  real(dp),dimension(1:nvector,1:nstore_disk),save::sig1   ! column dens etc
  real(dp),dimension(1:nvector,1:nstore_disk)::sigm   ! column dens etc

  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,6:8): Bleft,
  ! U(i,nvar+1:nvar+3): Bright
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,6:8): Bleft,
  ! Q(i,nvar+1:nvar+3): Bright
  ! If nvar > 8, remaining variables (9:nvar) are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::  ivar, idust, i
  real(dp):: xn,x0,sum_dust, RR,yn,zn
  real(dp):: rhayash,pi,thayash,cs,omega_kep,H,csback,sig0,sfive,rout,rin,delta_rho
  real(dp):: dens,Bz
  real(dp):: sigmacr,sigmafuv,zetafuv,zetarad,zetacr,xfuv,sigmasc,sigmaab,alph,bet,gammai,sigmave,sigmahayash
  real(dp):: eel, me, clum, alpha_dr, xe, zeta,zeta1
  real(dp) :: scale_dcol
  real(dp):: betaambd,etaohmic,eta_cap,beta_cap
  real(dp):: sinthetai, sintheta,alpha_disk,k_corona,v_kep,cs_iso,n_disk,buffer_H,radius,rho0,rrr

  !hayashi's params
  rin =rd_factor
  rout= 5.0d0!2.0*4.0
  H=HoverR*rout
   rho0=rhocen/scale_d
 
  csback= sqrt(gamma*kb*Tpback/(mu_gas*mH))/scale_v

  x0=boxlen/2.0
  pi=3.141592653589

  scale_dcol= scale_m/scale_l**2.0
  !resistivity params
  sigmacr=96./scale_dcol
  sigmafuv=0.03/scale_dcol
  zetafuv=2e-5
  zetarad=1e-19
  zetacr=1e-16
  sigmasc=7e23*mu_gas*mH/scale_dcol
  sigmaab=1e21*mu_gas*mH/scale_dcol
  alph=0.65
  bet=0.4
  gammai=3.5e13*scale_d*scale_t
  eel=4.8e-10
  me=9.1e-28
  clum=2.997e10
  if(hayashi.or.bethune) then
     rhayash=1.0d0
  endif
  do i=1,nn
     ! Compute position in normalized coordinates
     xn=(x(i,1)-x0)
     yn=(x(i,2)-x0)
     zn=(x(i,3)-x0)
     !cylindrical radius
     !spherical radius
     RR = sqrt(xn**2.0+yn**2.0+rsmooth**2.0)
     radius = sqrt(RR**2.0+zn**2.0)
     RRR = sqrt(xn**2.0+yn**2.0)


     if (bethune) then     
     alpha_disk=-1.
     n_disk=-2.
     k_corona= 6.
     THayash= 300.0*(rr/rhayash)**(alpha_disk/2.0)
     cs_iso=sqrt(gamma*kb*THayash/(mu_gas*mH))/scale_v
     H=hoverr*RRR
     buffer_H=0.1
     if(abs(zn).lt.3.72*H-H*buffer_H) then
     cs = cs_iso
     dens=rho0*(RR/rhayash)**n_disk*exp(-Mstar_cen/cs**2.0*(1/RR-1/radius))
     else
     if (abs(zn).gt.3.72*H) then
     cs=k_corona*cs_iso
     dens=rho0*(RR/rhayash)**n_disk*exp(-Mstar_cen/cs**2.0*(1/RR-1/radius))
     else
     cs= cs_iso*(k_corona+(1.0-k_corona)*(abs(zn)-3.72*H+H*buffer_H)/(H*buffer_H))
     dens=rho0*(RR/rhayash)**n_disk*exp(-Mstar_cen/cs**2.0*(1/RR-1/radius))
     endif
     endif

     omega_kep=sqrt(Mstar_cen/rr**3.0)

     Bz=sqrt(2.0d0*(cs_iso**2.0*rho0*(RR/rhayash)**n_disk/beta_mag))



     sig0=1700/sqrt(2.0*pi)*(rr/rhayash)**(-3./2.)/scale_m*scale_l**2.0
     sig1(i,1)= min(1700.d0/scale_m*scale_l**2.0*(rr/rhayash)**(-3./2.),sig0*(1.0d0-erf(abs(zn)/(sqrt(2.0)*H))))
     zeta1=(rr/rhayash)**(-2.2)*(1e-15*(exp(-(sig1(i,1)/sigmasc)**alph))+6e-12*(exp(-(sig1(i,1)/sigmaab)**bet)))
     zeta=zetacr*exp(-sig1(i,1)/sigmacr)+zetarad+zeta1
     xfuv=zetafuv*exp(-(sig1(i,1)/sigmafuv)**4.0)
     alpha_dr=3e-6/sqrt(Thayash)
     xe= sqrt(zeta/alpha_dr/(H2_fraction*dens*scale_d/(mu_gas*mH)))+xfuv
     sigmave=8.28e-9*sqrt(Thayash/100.)
     eta_cap=10*omega_kep*H**2.0
     beta_cap=eta_cap/bz**2.0
     etaohmic=min(clum**2.0*me/(4.0*pi*eel**2.0)/xe*sigmave/scale_l**2*scale_t,eta_cap)
     betaambd=min(1./gammai/dens**2/xe/29,beta_cap)
     sig1(i,2)=etaohmic
     sig1(i,3)=betaambd
!!$     sig1(i,4)=xe
!!$     sig1(i,5)=1./dens/betaambd/omega_kep
!!$     sig1(i,6)=Bz**2.0/dens/etaohmic/omega_kep
     ! print *, sig(i,2),sig(i,3)
  endif
  
  if (Hayashi) then
     omega_kep=sqrt(Mstar_cen/rr**3.0)

     sigmaHayash= 1700.d0/scale_m*scale_l**2.0*(rr/rhayash)**(-3./2.)
   
     THayash= 280.0d0*(rr/rhayash)**(-1./2.)
     cs =  sfive(rr/rsmooth)*sqrt(gamma*kb*THayash/(mu_gas*mH))/scale_v+csback
     H=cs/omega_kep

     dens= max(sigmaHayash/sqrt(2.0*pi*H**2.0)*exp(-zn**2.0/(2.0*H**2.0)),rhoext/scale_d)
     Bz=sqrt(2.0d0*(cs**2.0*sigmaHayash/sqrt(2.0*pi*H**2.0))/beta_mag)

     sig0=1700/sqrt(2.0*pi)*(rr/rhayash)**(-3./2.)/scale_m*scale_l**2.0
     sig1(i,1)= min(1700.d0/scale_m*scale_l**2.0*(rr/rhayash)**(-3./2.),sig0*(1.0d0-erf(abs(zn)/(sqrt(2.0)*H))))
     zeta1=(rr/rhayash)**(-2.2)*(1e-15*(exp(-(sig1(i,1)/sigmasc)**alph))+6e-12*(exp(-(sig1(i,1)/sigmaab)**bet)))
     zeta=zetacr*exp(-sig1(i,1)/sigmacr)+zetarad+zeta1
     xfuv=zetafuv*exp(-(sig1(i,1)/sigmafuv)**4.0)
     alpha_dr=3e-6/sqrt(Thayash)
     xe= sqrt(zeta/alpha_dr/(H2_fraction*dens*scale_d/(mu_gas*mH)))+xfuv
     sigmave=8.28e-9*sqrt(Thayash/100.)
     eta_cap=10*omega_kep*H**2.0
     beta_cap=eta_cap/bz**2.0
     etaohmic=min(clum**2.0*me/(4.0*pi*eel**2.0)/xe*sigmave/scale_l**2*scale_t,eta_cap)
     betaambd=min(1./gammai/dens**2/xe/29,beta_cap)
     sig1(i,2)=etaohmic
     sig1(i,3)=betaambd
!!$     sig1(i,4)=xe
!!$     sig1(i,5)=1./dens/betaambd/omega_kep
!!$     sig1(i,6)=Bz**2.0/dens/etaohmic/omega_kep
     ! print *, sig(i,2),sig(i,3)
  endif
  
  end do
sigm(1:nn,1)=sig1(1:nn,1)
sigm(1:nn,2)=sig1(1:nn,2)
sigm(1:nn,3)=sig1(1:nn,3)
!!$sig(1:nn,4)=sig1(1:nn,4)
!!$sig(1:nn,5)=sig1(1:nn,5)
!!$sig(1:nn,6)=sig1(1:nn,6)

end subroutine columninit
!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

!!$     xx=x(i,1)
!!$#if NDIM > 1
!!$     yy=x(i,2)
!!$#endif
!!$#if NDIM > 2
!!$     zz=x(i,3)
!!$#endif
!!$     ! ABC
!!$     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
!!$     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
!!$     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega

!!$     v(i,1)=vx
!!$#if NDIM > 1
!!$     v(i,2)=vy
!!$#endif
!!$#if NDIM > 2
!!$     v(i,3)=vz
!!$#endif

  end do


end subroutine velana

subroutine get_vturb(vrms,cs,v_turb)

  use amr_commons, only:myid,ncpu
  use pm_commons, only:localseed,IRandNumSize,iseed
  use random

  implicit none

  integer :: i
  integer ,dimension(1:ncpu,1:IRandNumSize)    :: allseed
  double precision ::  vrms, cs, vel
  double precision, dimension(3) :: v_turb

  double precision :: u1, v1, u2, v2, theta, phi, x, y, z

#ifdef DEBUGRANDOM
  logical, save :: first=.true.
#endif

  if (localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if

  ! magnitude --> Gressel is v_rms = 0.01 * cs
  if (vrms .lt. 0.d0) then
    call gaussdev(localseed,vel)
    vel = vel * abs(vrms)*cs
  else
    call gaussdev(localseed,vel)
    vel = vel * vrms
  endif

  call ranf(localseed,v_turb(1))
  call ranf(localseed,v_turb(2))
  call ranf(localseed,v_turb(3))

  v_turb = (v_turb - 0.5d0)*2.d0
  v_turb = v_turb/sqrt(sum((v_turb(1:3))**2)) * vel

  ! NEED TO HAVE Cs OR Vrms in code units.

#ifdef DEBUGRANDOM
  if (myid == 1 .and. first) then
    open(42,file="test_gauss.dat",status='unknown')
    do i=1,10000
      call gaussdev(localseed,vel)
      write(42,*) vel
    enddo
    close(42)
    first = .false.
  endif
#endif

end subroutine get_vturb


subroutine get_dturb(pert,rho0,delrho)

  use amr_commons, only:myid,ncpu
  use pm_commons, only:localseed,IRandNumSize,iseed
  use random

  implicit none

  integer :: i
  integer ,dimension(1:ncpu,1:IRandNumSize)    :: allseed
  double precision ::  pert,rho0,delrho,randno


#ifdef DEBUGRANDOM
  logical, save :: first=.true.
#endif

  if (localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if

  call ranf(localseed,randno)

  delrho = randno*pert*rho0


end subroutine get_dturb
