subroutine init_dust_ratio(dustratio,epsilondust)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp) :: dustratio
  real(dp), dimension(1:ndust):: epsilondust
  real(dp), dimension(1:ndust+1):: sdust
  real(dp) :: epsilon_0,Anorm
  integer  :: idust
  epsilon_0 = dustratio/(1.0d0+dustratio)
  Anorm = 1.0d0/(size_max**(4.0d0-mrn_index)-size_min**(4.0d0-mrn_index))
  do idust =1,ndust+1
     sdust(idust) = size_min+(size_max-size_min)*DBLE(idust-1)/DBLE(Ndust)
  enddo
  do idust=1,ndust
     epsilondust(idust)= epsilon_0*Anorm*(sdust(idust+1)**(4.0d0-mrn_index)-sdust(idust)**(4.0d0-mrn_index))
  enddo
end subroutine init_dust_ratio

subroutine size_dust(sdust)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp), dimension(1:ndust):: sdust
  real(dp), dimension(1:ndust+1):: sdust_interval
  real(dp) :: Bnorm
  integer  :: idust
  Bnorm = 1.0d0/(size_max**(1.0d0-mrn_index)-size_min**(1.0d0-mrn_index))
  do idust =1,ndust+1
     sdust_interval(idust) = size_min+(size_max-size_min)*DBLE(idust-1)/DBLE(Ndust)
  enddo
   !We compute the average dust size in the bin to get the correct stopping time 
   do idust =1,ndust
      sdust(idust) = (0.5d0*(sdust_interval(idust)+sdust_interval(idust+1)))**(-1.0d0/mrn_index)*Bnorm**(1.0/mrn_index)
   enddo
end subroutine size_dust


