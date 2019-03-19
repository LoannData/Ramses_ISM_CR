module hydro_commons
  use amr_parameters
  use hydro_parameters
  real(dp),allocatable,dimension(:,:)::uold,unew ! State vector and its update
  real(dp),allocatable,dimension(:)  :: pstarold,pstarnew ! Stellar momentum and its update
  real(dp),allocatable,dimension(:)  ::divu,enew ! Non conservative variables
  real(dp),allocatable,dimension(:,:)::Rho_eos,Ener_eos,Temp_eos,P_eos,Cs_eos,S_eos,eint_eos
  real(dp),allocatable,dimension(:,:)::xH_eos, xH2_eos, xHe_eos,xHep_eos,Cv_eos,Dc_eos
  real(dp),allocatable,dimension(:,:)::resistivite_chimie_res ! resistivites chimie
  real(dp),allocatable,dimension(:,:,:,:)::resistivite_chimie_x ! resistivites chimie
#if NDUST>0
  real(dp),allocatable,dimension(:):: dflux_dust_eint ! Dust flux from refined interfaces for enint

  real(dp),allocatable,dimension(:,:):: dflux_dust ! Dust flux from refined interfaces
  real(dp),allocatable,dimension(:,:,:):: v_dust ! Differential velocity of the dust phases
  integer::jdust
#endif
#if MC>0  
  real(dp),allocatable,dimension(:,:)::fluxes ! Mass flux on the faces of the cells
#endif  
  real(dp),allocatable,dimension(:)::rho_barotrop,temp_barotrop
  real(dp)::mass_tot=0.0D0,mass_tot_0=0.0D0
  real(dp)::ana_xmi,ana_xma,ana_ymi,ana_yma,ana_zmi,ana_zma
  integer::nbins
  integer,allocatable,dimension(:)::liste_ind
  integer::nb_ind
  real(dp)   ::dt_imp                            ! Implicit timestep   
end module hydro_commons

module const
  use amr_parameters
  real(dp)::bigreal = 1.0e+30
  real(dp)::zero    =  0.0d0
  real(dp)::one     =  1.0d0
  real(dp)::two     =  2.0d0
  real(dp)::three   =  3.0d0
  real(dp)::four    =  4.0d0
  real(dp)::five    =  5.0d0
  real(dp)::six     =  6.0d0
  real(dp)::seven   =  7.0d0
  real(dp)::eight   =  8.0d0
  real(dp)::nine    =  9.0d0
  real(dp)::ten     = 10.0d0
  real(dp)::two3rd  = 0.6666666666666667d0
  real(dp)::half    = 0.5d0
  real(dp)::third   = 0.33333333333333333d0
  real(dp)::forth   = 0.25d0
  real(dp)::sixth   = 0.16666666666666667d0
end module const

! Units
module units_commons
  use amr_parameters, only : dp
  real(dp):: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp):: scale_E0,scale_kappa,scale_m,P_cal,C_cal
end module units_commons
