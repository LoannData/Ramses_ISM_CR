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
  real(dp),allocatable,dimension(:,:):: dflux_dust ! Dust flux from refined interfaces
  real(dp),allocatable,dimension(:,:,:):: v_dust ! Differential velocity of the dust phases
  !real(dp),allocatable,dimension(:,:,:):: v_dust_0 ! Initial differential velocity of the dust phases, to compute decay

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
  real(dp)::bigreal = 1.0e+30_dp
  real(dp)::zero    =  0.0_dp
  real(dp)::one     =  1.0_dp
  real(dp)::two     =  2.0_dp
  real(dp)::three   =  3.0_dp
  real(dp)::four    =  4.0_dp
  real(dp)::five    =  5.0_dp
  real(dp)::six     =  6.0_dp
  real(dp)::seven   =  7.0_dp
  real(dp)::eight   =  8.0_dp
  real(dp)::nine    =  9.0_dp
  real(dp)::ten     = 10.0_dp
  real(dp)::two3rd  = 0.6666666666666667_dp
  real(dp)::half    = 0.5_dp
  real(dp)::third   = 0.33333333333333333_dp
  real(dp)::forth   = 0.25_dp
  real(dp)::sixth   = 0.16666666666666667_dp
end module const

! Units
module units_commons
  use amr_parameters, only : dp
  real(dp):: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp):: scale_E0,scale_kappa,scale_m,P_cal,C_cal
end module units_commons
