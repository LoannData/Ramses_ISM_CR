subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  !use hydro_parameters, only : Ndust
  use radiation_parameters,only:mu_gas
  use cooling_module
  use units_commons, only : scale_kappa,scale_m

  implicit none

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l,pc,au,pi, msol
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------
      pc = 3.08e18_dp
      au = 1.5d13
      Msol =2d33
      pi=  3.14159256358


      scale_d = 1d-10
      scale_T2=1.0_dp
      scale_l = au
      scale_t = 3600.0d0
      scale_v = scale_l/scale_t
      scale_nH =1.0_dp
      scale_kappa = 1.0_dp/scale_l
      scale_m = scale_d*scale_l**3




end subroutine units
