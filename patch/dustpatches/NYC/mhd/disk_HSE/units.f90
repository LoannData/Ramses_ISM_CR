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
    !  scale_d = mu_gas * mH
    !  scale_T2=1.0_dp
    !  scale_l = pc
    !  scale_t = 3.0_dp*3.14159265358979323846_dp/SQRT (32.0_dp*Grav*scale_d)
    !  scale_v = scale_l/ scale_t
    !  scale_nH = X/mH * scale_d
    !  scale_kappa = 1.0_dp/scale_l
    !  scale_m = scale_d*scale_l**3


      scale_d = 1.0_dp
      scale_T2=1.0_dp
      scale_l = 1.0_dp
      scale_t = 1.0_dp
      scale_v = 1.0_dp
      scale_nH =1.0_dp
      scale_kappa = 1.0_dp/scale_l
      scale_m = scale_d*scale_l**3

      au = 1.5d13
      Msol =2d33
      pi=  3.14159256358
      scale_l = au
      scale_m = Msol
      scale_d = scale_m/(scale_l**3)
      scale_t = 2.0*pi*(300*au)**(3.0/2.0)/sqrt(Grav*Msol)
      scale_v = scale_l / scale_t
      scale_T2 = mu_gas**2 * mH**2 * au**2 * Grav / kb
      scale_nH = 1.0_dp ! X/(mH) * scale_d
      scale_kappa = 1.0_dp/scale_l


end subroutine units
