subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  use units_commons, only : scale_kappa,scale_m
  use radiation_parameters, only : mu_gas

  implicit none

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  real(dp)::au, Msol,pi

  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  !code units are: G=1, rho cm^-3, x pc
  au = 1.5d13
  Msol =2d33
  pi=  3.14159256358d0
  scale_l = au
  scale_m = Msol
  scale_d = scale_m/(scale_l**3.0d0)
  scale_t = sqrt(scale_l**3.0d0/grav/scale_m)
  scale_v = scale_l / scale_t
  scale_T2 = 1.0d0
  scale_nH = 1.0_dp ! X/(mH) * scale_d
  scale_kappa = 1.0_dp/scale_l




end subroutine units
