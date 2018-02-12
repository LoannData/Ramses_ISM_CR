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
  sdust(1)=size_min
  do idust =1,ndust
     sdust(idust+1) = sdust(idust)+10**(log10(size_max/size_min)/dble(ndust))
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
  integer  :: idust
  do idust =1,ndust
     sdust(idust+1) = sdust(idust)+10**(log10(size_max/size_min)/dble(ndust))
  enddo
   !We compute the average dust size in the bin to get the correct stopping time 
   do idust =1,ndust
      sdust(idust) = sqrt(sdust_interval(idust)*sdust_interval(idust+1))
   enddo
end subroutine size_dust


