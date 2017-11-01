subroutine init_dust_distrib(dust_zone,epsilon_loc,sum_dust,idust)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp) :: dust_zone
  real(dp) :: epsilon_loc, sum_dust
  !Distrib related quantities, normalisation coeff
  integer ::idust  
  real :: dust_bin, size_min_bin,size_max_bin, norm_dust
if(ndust.eq. 1) then 
   !epsilon_loc = dust_region
else
   norm_dust=(powerlaw_dust-1.0_dp)/(size_grain**(1.0_dp-powerlaw_dust)-size_grain_max**(1.0_dp-powerlaw_dust))
   dust_bin = (size_grain_max- size_grain)/DBLE(ndust)
   size_min_bin = size_grain + DBLE(idust)*dust_bin
   size_max_bin = size_min_bin + dust_bin
   !epsilon_loc  = dust_region*norm_dust*&
    !&(powerlaw_dust-1.0_dp)/(size_min_bin**(1.0_dp-powerlaw_dust)-size_max_bin**(1.0_dp-powerlaw_dust))
endif 
!sum_dust    = sum_dust+ dust_region
   
end subroutine init_dust_distrib
