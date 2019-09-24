subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  use units_commons, only : scale_kappa,scale_m

  implicit none

  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  ! scale_d converts mass density from user units into g/cc
  scale_d = units_density
  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3

  ! scale_t converts time from user units into seconds
  scale_t = units_time
  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)

  ! scale_l converts distance from user units into cm
  scale_l = units_length
  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)

  ! scale_v converts velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  !scale_T2 = mH/kB * scale_v**2
  scale_T2 = 1.0_dp ! T0 !scale_v**2 * mu_gas

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

  scale_kappa = 1.0_dp/scale_l

  scale_m = scale_d*scale_l**3
  
#if NIMHD==1
  ! modif nimhd
  if(ntestDADM.eq.1) then
      scale_d = 1.d0
      scale_l = 1.d0
      scale_v = 1.d0
      scale_t = 1.d0
  end if
  ! fin modif nimhd
#endif

end subroutine units
