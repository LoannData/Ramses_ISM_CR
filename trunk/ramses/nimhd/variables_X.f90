module variables_X
   use amr_parameters, ONLY:dp
   implicit none

   !integer,parameter   :: alpha=10           ! number of bins
   integer   :: alpha           ! number of bins
   !integer,parameter   :: Nvarchimie=7+(3*alpha)   ! nombre d'especes ioniques
   integer   :: Nvarchimie   ! nombre d'especes ioniques
   !integer,parameter   :: Ntot=Nvar+6  ! nombre total d'especes considerees

   !!!! Parametre d'integration
   !real(dp)  :: nH    ! densite en cm-3 d'atomes d'hydrogene

   !!!! Variables
   real(dp), allocatable, dimension(:) :: x_g        ! abondance des grains par bins
   real(dp), allocatable, dimension(:) :: r_g        ! rayon des grains par bins
   real(dp), allocatable, dimension(:) :: m_g        ! masse des grains par bins


   !!!! Pour les collisions avec les grains
   !real(dp), parameter :: pi=acos(-1.)
   real(dp), parameter :: pi=3.1415927410125732422_dp  
   real(dp), parameter :: rg=0.1d-4   ! rayon des grains (en cm)
   real(dp), parameter :: mp=1.6726d-24 ! masse du proton, en g
   real(dp), parameter :: me=9.1094d-28 ! masse de l'electron, en g
   real(dp), parameter :: mg=1.2566d-14 ! masse des grains, en g
   real(dp), parameter :: e=4.803204d-10 ! charge de l'electron, en cgs
   real(dp), parameter :: mol_ion=29._dp*mp    ! molecular ion, en masse de proton
   real(dp), parameter :: Met_ion=23.5_dp*mp    ! atomic ion, en masse de proton
   real(dp), parameter :: Kb = 1.3807d-16    ! erg/K
   real(dp), allocatable, dimension(:) :: q           ! charge 
   real(dp), allocatable, dimension(:) :: m           ! masse

   !!!! bins de grains
   real(dp), parameter :: rho_s=2.3_dp        ! g.cc
   real(dp), parameter :: rho_n_tot=1.17d-21        ! g.cc
   real(dp), parameter :: a_0=0.0375d-4     ! cm
   real(dp), parameter :: a_min=0.0181d-4   ! cm
   real(dp), parameter :: a_max=0.9049d-4   ! cm
   real(dp), parameter :: zeta=a_min/a_max  ! a_min/a_max

   real(dp)            :: rho_gtot          
   !real(dp), parameter :: P_e=0.6d0
   !real(dp), parameter :: P_I=1.0d0


   ! resistivites (cf Kunz & Mouschovias 2009)
   real(dp), allocatable, dimension(:)  :: sigma         ! sigma_s
   real(dp), allocatable, dimension(:)  :: zetas
   real(dp), allocatable, dimension(:)  :: phi
   real(dp), allocatable, dimension(:)  :: tau_sn,tau_gp,tau_gm
   real(dp), allocatable, dimension(:,:)  :: tau_inel
   real(dp), allocatable, dimension(:)  :: gamma_zeta,gamma_omega
   real(dp), allocatable, dimension(:)  :: omega,omega_bar
   real(dp), parameter        :: c_l=299792458.d2     ! c en cm/s

   integer :: cmp_X

   ! for RAMSES
   integer :: tchimie   ! tchmie steps in temperature
   real(dp) :: dtchimie   ! dtchmie step in temperature
   real(dp) :: tminchimie   ! min tchmie 
   integer :: nchimie   ! nchmie steps in density 
   real(dp) :: dnchimie   ! dnchmie step in density 
   real(dp) :: nminchimie   ! min nchmie  
   integer :: bchimie   ! bchmie steps in B field
   real(dp) :: dbchimie   ! dbchmie step in B field
   real(dp) :: bminchimie   ! min bchmie 
   real :: alpha_real
   real(dp),allocatable,dimension(:,:,:,:)::resistivite_chimie ! resistivites chimie

end module variables_X

subroutine rq    ! this routine is not used in the code, just use it once if you change the grains distribution to update the variables q(:) and r_g(:)
   use variables_X
   implicit none
   integer :: i

   !rho_gtot=0.01d0*(a_min/a_0)*zeta**(-0.5)*rho_n_tot

   !cmp_X=0

   alpha_real=(nvarchimie-7)/3
   alpha=floor(alpha_real)

   if (alpha_real.ne.real(alpha)) then
      print*, 'issue in number of species'
      stop
   endif

   allocate(x_g(alpha))
   allocate(r_g(alpha))
   allocate(m_g(alpha))
   allocate(q(nvarchimie))
   allocate(m(nvarchimie))
   allocate(sigma(nvarchimie))
   allocate(zetas(nvarchimie))
   allocate(phi(nvarchimie))
   allocate(tau_sn(nvarchimie))
   allocate(tau_gp(nvarchimie))
   allocate(tau_gm(nvarchimie))
   allocate(gamma_zeta(nvarchimie))
   allocate(gamma_omega(nvarchimie))
   allocate(omega(nvarchimie))
   allocate(omega_bar(nvarchimie))
   allocate(tau_inel(nvarchimie,nvarchimie))
   

   if  (alpha==1) then
      !    rho_gtot=0.01d0*rho_n_tot
      !    xg=1.d0/nH*(rho_gtot/(4.d0/3.d0*pi*rho_s*a_0**3))
      !    x(8)=xg/3.d0
      !    x(9)=xg/3.d0
      !    x(10)=xg/3.d0
      r_g(1)=a_0
      !xg=8.9d-13
      !x(8)=2.d-17
      !x(9)=8.d-13
      !x(10)=xg-x(8)-x(9)
      !r_g(1)=1.d-5
   else
      !xg=1.d0/3.d2*(rho_gtot/(4.d0/3.d0*rho_s*pi*a_min**3))*&
      ! & (1.d0/5.d0*(1.d0-zeta**2.5)/(1.d0-zeta**0.5))*zeta**0.5
      do  i=1,alpha    ! cf Kunz & Mouschovias ou mes notes
         !    x_g(i)=xg*zeta**(2.5*(i-1.d0)/alpha)*((1.d0-zeta**(2.5/alpha))/(1.d0-zeta**2.5))
         r_g(i)=a_min*zeta**(-(i-1.d0)/alpha)*(5.d0*(1.d0-zeta**(0.5/alpha))/(1.d0-zeta**(2.5/alpha)))**0.5
         !    x(8+3*(i-1))=x_g(i)/3.d0
         !    x(9+3*(i-1))=x_g(i)/3.d0
         !    x(10+3*(i-1))=x_g(i)/3.d0
      end do
   end if



   !!!! pour les grains
   ! (i) particule chargée + grain neutre
   q(:)=1.d0*e    ! cations
   q(1)=-1.d0*e   ! electron
   ! i>7
   do  i=8,Nvarchimie
      if (mod(i-8,3)==0) q(i)=1.d0*e   ! g+
      if (mod(i-8,3)==1) q(i)=-1.d0*e  ! g-
      if (mod(i-8,3)==2) q(i)=0.d0     ! g0
   end do
   m(:)=0.d0
   m(1) = me           ! e-
   m(2) = 23.5d0*mp    ! ions metalliques
   m(3) = 29.d0*mp     ! ions moleculaires
   m(4) = 3*mp         ! H3+
   m(5) = mp           ! H+
   m(6) = 12.d0*mp     ! C+
   m(7) = 4.d0*mp      ! He+
   do  i=1,alpha       ! masse des grains
      m_g(i)=4.d0/3.d0*pi*r_g(i)**3*rho_s
      !print*, m_g(i)
      m(8+3*(i-1))=m_g(i)
      m(9+3*(i-1))=m_g(i)
      m(10+3*(i-1))=m_g(i)
   end do


end subroutine rq


subroutine nimhd_3dtable
   use variables_X
   use hydro_commons, only:resistivite_chimie_x
   use amr_commons, only : myid
   implicit none

   integer  :: iB,iH,iT,i
   real(dp) :: B,bmaxchimie,nH,T,sigH,sigO,sigP
   real(dp), dimension(nvarchimie) :: x

   if(myid==1) write(*,*) 'Computing 3D resistivities table'
   
   ! values for Btable
   bminchimie=1.d-10
   bmaxchimie=1.d5               ! ok for first core in nimhd. maybe not enough for second core.
   bchimie=100
   dbchimie=(log10(bmaxchimie)-log10(bminchimie))/real((bchimie-1),dp)

   allocate(resistivite_chimie(0:3,1:nchimie,1:tchimie,1:bchimie))

   tau_sn      = 0.0_dp
   omega       = 0.0_dp
   sigma       = 0.0_dp
   phi         = 0.0_dp
   zetas       = 0.0_dp
   gamma_zeta  = 0.0_dp
   gamma_omega = 0.0_dp
   omega_bar   = 0.0_dp
      
   do  iB=1,bchimie
   do  iT=1,tchimie
   do  iH=1,nchimie

      nh=resistivite_chimie_x(0,iH,iT)  ! density (.cc) of current point
      B =10.0d0**(log10(bminchimie)+dble(iB-1)*dbchimie)
      T =resistivite_chimie_x(-1,iH,iT)
      
!      write(*,*) ih,it,ib,nh,b,t
!      read(*,*)
      !inp=nh/2.d0/H2_fraction     ! inp is neutrals.cc, to fit densionbis
      !print *, q(1),m(1)
      !stop
      !print*, nH,'dummy',barotrop1D(nH*scale_d)
      do  i=1,nvarchimie-3*alpha
         if  (i==1) then  ! electron
            x(i) = resistivite_chimie_x(i,iH,iT)
            tau_sn(i) = 1.d0/1.16d0*(m(i)+2.d0*mp)/(2.d0*mp)*1.d0/(nH/2.d0*1.3d-9)
            omega(i) = q(i)*B/(m(i)*c_l)
            sigma(i) = x(i)*nH*(q(i))**2*tau_sn(i)/m(i)
            phi(i) = 0.d0
            zetas(i) = 0.d0
            gamma_zeta(i) = 1.d0
            gamma_omega(i) = 1.d0
            omega_bar(i) = 0.d0
         else if (i>=2 .and. i<=7) then ! ions
            x(i) = resistivite_chimie_x(i,iH,iT)
            tau_sn(i) = 1.d0/1.14d0*(m(i)+2.d0*mp)/(2.d0*mp)*1.d0/(nH/2.d0*1.69d-9)
            omega(i) = q(i)*B/(m(i)*c_l)
            sigma(i) = x(i)*nH*(q(i))**2*tau_sn(i)/m(i)
            phi(i) = 0.d0
            zetas(i) = 0.d0
            gamma_zeta(i) = 1.d0
            gamma_omega(i) = 1.d0
            omega_bar(i) = 0.d0
         end if
      end do
      do  i=1,alpha   ! grains
         x(8+3*(i-1)) = resistivite_chimie_x(8+3*(i-1),iH,iT)
         tau_sn(8+3*(i-1))=1.d0/1.28d0*(m_g(i)+2.d0*mp)/(2.d0*mp)*1.d0/(nH/2.d0*(pi*r_g(i)**2*(8.d0*Kb*T/(pi*2.d0*mp))**0.5))    ! g+
         omega(8+3*(i-1)) = q(8+3*(i-1))*B/(m_g(i)*c_l)
         sigma(8+3*(i-1)) = x(8+3*(i-1))*nH*(q(8+3*(i-1)))**2*tau_sn(8+3*(i-1))/m_g(i)

         x(9+3*(i-1)) = resistivite_chimie_x(9+3*(i-1),iH,iT)
         tau_sn(9+3*(i-1))=tau_sn(8+3*(i-1))          ! g-
         omega(9+3*(i-1)) = q(9+3*(i-1))*B/(m_g(i)*c_l)
         sigma(9+3*(i-1)) = x(9+3*(i-1))*nH*(q(9+3*(i-1)))**2*tau_sn(9+3*(i-1))/m_g(i)

         x(10+3*(i-1)) = resistivite_chimie_x(10+3*(i-1),iH,iT)
         tau_sn(10+3*(i-1))=tau_sn(8+3*(i-1))          ! g0
         omega(10+3*(i-1)) = q(10+3*(i-1))*B/(m_g(i)*c_l)
         sigma(10+3*(i-1)) = x(10+3*(i-1))*nH*(q(10+3*(i-1)))**2*tau_sn(10+3*(i-1))/m_g(i)

      end do

      sigP=0.d0
      sigO=0.d0
      sigH=0.d0

      do i=1,nvarchimie
         sigP=sigP+sigma(i)
         sigO=sigO+sigma(i)/(1.d0+(omega(i)*tau_sn(i))**2)
         sigH=sigH-sigma(i)*omega(i)*tau_sn(i)/(1.d0+(omega(i)*tau_sn(i))**2)
      end do
      resistivite_chimie(1,iH,iT,iB)=log10(sigP)
      resistivite_chimie(2,iH,iT,iB)=log10(sigO)
      resistivite_chimie(3,iH,iT,iB)=log10(abs(sigH))
      resistivite_chimie(0,iH,iT,iB)=sign(1.0d0,sigH)
   end do
end do
end do
   
   if(myid==1) write(*,*) '3D resistivities table complete'

end subroutine nimhd_3dtable

