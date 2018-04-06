subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  use units_commons, only : scale_kappa,scale_m
  use radiation_parameters, only : mu_gas

  implicit none

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  real(dp)::pc
  
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  !code units are: G=1, rho cm^-3, x pc
  pc = 3.08e18_dp

  ! scale_d converts mass density from user units into g/cc
  !scale_d = units_density
  !if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3
  !calculate the initial density
  scale_d = 1.d-20

  ! scale_t converts time from user units into seconds
  !scale_t = units_time
  !if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)
  ! scale_t converts time from user units into seconds
  scale_t = 1.0_dp/sqrt(Grav*scale_d)

  ! scale_l converts distance from user units into cm
  !scale_l = units_length
  !if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)
  !calculate the initial cloud radius
  scale_l = 1.496d13

  ! scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  !kT = 10. * kB  / G / pc^2 / mpart^2
  scale_T2 = mu_gas**2 * mH**2 * pc**2 * Grav / kb

  ! scale_nH converts rho in user units into nH in H/cc
!  scale_nH = X/(mH*mu) * scale_d
  scale_nH = 1.0_dp ! X/(mH) * scale_d
  
  scale_kappa = 1.0_dp/scale_l
  
  scale_m = scale_d*scale_l**3
  
#if NIMHD==1
  ! modif nimhd
  if(ntestDADM.eq.1) then
      scale_d = 1.0_dp
      scale_l = 1.0_dp
      scale_v = 1.d0
      scale_t = 1.0_dp
  end if
  ! fin modif nimhd
#endif
  
end subroutine units

