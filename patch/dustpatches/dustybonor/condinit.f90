!================================================================
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_commons
  use radiation_parameters  
  use poisson_parameters
  use const
  use cooling_module,ONLY:kB,mH,clight
  use cloud_module
  
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  
  integer::i,dim,j,id,iu,iv,iw,ip,iEr,ie
  real(dp):: pi,xc,yc,zc,rs,rc,x0,y0,z0,r0,d0,B0,p0,omega_c,x1,y1,z1
  real(dp)::Mass,temp,vit_son,vit_son_norm
  real(dp)::G,mu,psic,xi2psipb,G_norm,xi2psipc,psib
  real(dp)::betab, megayear,pc , tff,alpha,xic
  real(dp)::dens_bord,dens_centre,rayon_phys,rayon_bord,rayon_norm
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  character(LEN=60)::file_bonnor
  real(dp),allocatable,dimension(:)::xi,psi,phi,xi2psip,xi2empsi,dif2,rb,dif
  integer,dimension(1)::mindif
  real(dp)::a,b,c,d,e, xi2psip_nouv,rbb,rtr, phi_p, delta_p, a_R,mnuage,msubH
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
  integer::ivar
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp):: sum_dust
#if NDUST>0
  integer:: idust
  real(dp):: epsilon_0
  real(dp),dimension(1:ndust):: dustMRN
  epsilon_0 = dust_ratio(1)
#endif
  
  ! Call built-in initial condition generator
!  call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions
  ! ........
  id=1; iu=2; iv=3; iw=4; ip=5; iEr=9; ie=nvar
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  x0=0.5*boxlen
  y0=0.5*boxlen
  z0=0.5*boxlen
  
  pc=3.08d18
  pi=ACOS(-1.0d0)
  G=6.67d-8
  a_r=7.5657d-15
  msubH=mu_gas*mH

  mnuage=mass_c*msun
  temp=T2_star*mu_gas
 
!  xic=7.0d0
  betab= kB*temp/(4.0d0*pi*G*msubH)


  vit_son = sqrt(kB *temp / msubH )


  file_bonnor='isodata.asc'
  open(10,file=file_bonnor,status='old')
  
  i=1
  do 
     read(10,66,end=100)a,b,c,d,e
     i=i+1
  enddo

100 dim=i-1
  
  rewind(10)
  allocate(xi(1:dim),psi(1:dim),phi(1:dim),xi2psip(1:dim),xi2empsi(1:dim), &
       & rb(1:dim), dif(1:dim))
  
 ! mindif=100.
  do i=1,dim
     read(10,66)xi(i),psi(i),phi(i),xi2psip(i),xi2empsi(i)
     dif(i) =  abs(xi(i) - xic)
   !  if(dif .lt. dif_min) mindif=dif
     rb(i)=xi(i)*(betab)**0.5
  
  enddo


  mindif =  minloc(dif)

  xi2psipb = xi2psip(mindif(1))  ! les parametres au bord de la sphere isotherme

  psib = psi(mindif(1))
 
  rbb= rb(mindif(1))
  !parametre de la sphere isotherme
  
  dens_centre = ((vit_son**6.)  * (xi2psipb**2.)) &
       & /((4. * pi)* (G**3.) * (mnuage**2.))
 ! print*,dens_centre,mnuage,xi2psipb,xic
  dens_bord = dens_centre * exp(-psib)

  P0= temp*dens_centre *kB /msubH 


! brot=90.

  tff= sqrt( 3. * pi / ( 32. * dens_centre / scale_d))
  omega_c = omega_machida * sqrt(4.*pi)!2. * pi / (brot * tff)
 
  rayon_phys = G * mnuage / (vit_son**2.) !/ scale_l
  rayon_bord = rayon_phys / xi2psipb * xic
  rayon_norm= rayon_bord / scale_l
  alpha = P0/( G * (dens_centre**2.)*8.*pi* rayon_bord **2./15.)
 
  do i=1,nn
     
     xc=x(i,1)-x0
     yc=x(i,2)-y0
     zc=x(i,3)-z0
     rs = sqrt(xc**2. + yc**2. + zc**2.)
     rc = sqrt(xc**2. + yc**2.)
 
     B0 = sqrt(4.*pi/5.)/0.53*crit

     !Bx component
     q(i,6     ) = 0.
     q(i,nvar+1) = 0.

     !By component
     q(i,7     ) = 0.
     q(i,nvar+2) = 0.

     !Bz component
     q(i,8     ) = B0
     q(i,nvar+3) = B0


     if( rs .lt. rayon_norm ) then 

        rtr=rs / rayon_norm * rbb
        dif(1:dim) =  abs(rb(1:dim) - rtr)
        mindif =  minloc(dif)

       sum_dust = 0.0d0
#if NDUST>0
        if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
        do idust =1,ndust
           q(i, firstindex_ndust+idust)= dust_ratio(idust)/(1.0d0+dust_ratio(idust))
           if(mrn) q(i, firstindex_ndust+idust)= dustMRN(idust)

           sum_dust = sum_dust + q(i, firstindex_ndust+idust)
        end do   
#endif
        phi_p      = atan(yc/xc)
        delta_p  = delta_rho*cos(2.0d0*phi_p)*(rs/rayon_norm)**2.
        q(i,id)  = dens_centre * exp(-psi(mindif(1))) / scale_d
        q(i,id)  = q(i,id) * f_dens * (1.0 + delta_p) 
        q(i,iu)  = omega_c * yc 
        q(i,iv)  = -omega_c * xc
        q(i,iw)  = 0.0
        q(i,ip)  = temp * q(i,id) *(1.0-sum_dust)* kB /msubH /( scale_l**2. / scale_t**2.)
        q(i,iEr) = a_R*Temp**4./(scale_d*scale_v**2)
        q(i,ie)  = q(i,ip)/(gamma-1.0d0)/q(i,id)

     ELSE
       sum_dust = 0.0d0
#if NDUST>0
        if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
        do idust =1,ndust
           q(i, firstindex_ndust+idust)= dust_ratio(idust)/(1.0d0+dust_ratio(idust))
           if(mrn) q(i, firstindex_ndust+idust)= dustMRN(idust)

           sum_dust = sum_dust + q(i, firstindex_ndust+idust)
        end do   
#endif        
        q(i,id)  = dens_bord / scale_d !/ 100.
        q(i,iu)  = omega_c * yc 
        q(i,iv)  = -omega_c * xc
        q(i,iw)  = 0.0
        q(i,ip)  = temp * q(i,id) *(1.0-sum_dust)*  kB /msubH /( scale_l**2. / scale_t**2.)
        q(i,iEr) = a_R*Temp**4./(scale_d*scale_v**2)
        q(i,ie)  = q(i,ip)/(gamma-1.0d0)/q(i,id)
     ENDIF 

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
  !kinetic + magnetic energy
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125*(q(1:nn,8)+q(1:nn,nvar+3))**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic field 
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
  u(1:nn,9)=q(1:nn,9)
  u(1:nn,5)=u(1:nn,5)+u(1:nn,9)
  ! passive scalars
  do ivar=9+nrad,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#if NDUST>0
     ! dust
     do idust=1,ndust
        u(1:nn,firstindex_ndust+idust)=q(1:nn,1)*q(1:nn,firstindex_ndust+idust)
     end do
#endif  
  close(10)
66    format (12(e13.6,3x))
end subroutine condinit
!=============================================================================
!=============================================================================
!=============================================================================
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
  twopi = 2.0*acos(-1.0d0)
  
  do i=1,ncell

     xx=x(i,1)
     yy=x(i,2)
     zz=x(i,3)

     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

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
     
     v(i,1)=vx
     v(i,2)=vy
     v(i,3)=vz

  end do


end subroutine velana
