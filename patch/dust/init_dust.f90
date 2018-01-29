subroutine init_dust_ratio(dustratio,epsilondust)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp) :: dustratio
  real(dp), dimension(1:ndust):: epsilondust
  real(dp), dimension(1:ndust+1):: sdustbords

  real(dp) :: epsilon_0
  real(dp) :: smin,smax
  real(dp) :: mrn_index
  integer  :: idust
  epsilon_0 = dustratio/(1.0d0+dustratio)
  mrn_index = 3.5d0
  smin = 5.0d-7
  smax = 2.5d-5
  
  !We bin the distribution
  do idust =1,ndust+1
     sdustbords(idust) = 10.0D0**(log10(smax/smin)*log10(dble(idust))/log10(dble(ndust))+log10(smin))
  enddo
  ! We initialise it
  do idust=1,ndust
     epsilondust(idust)= epsilon_0*(sdustbords(idust+1)**(1.0d0-mrn_index)-sdustbords(idust)**(1.0d0-mrn_index))/ (smax**(1.0d0-mrn_index)-smin**(1.0d0-mrn_index))
  enddo
end subroutine init_dust_ratio

subroutine size_dust(sdust)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp), dimension(1:ndust):: sdust
  real(dp), dimension(1:ndust+1):: sdustbords
  real(dp) :: smin,smax
  real(dp) :: mrn_index
 
  integer  :: idust
  mrn_index = 3.5d0
  smin = 5.0d-7
  smax = 2.5d-5 
  !We bin the distribution
  do idust =1,ndust+1
     sdustbords(idust) = 10.0D0**(log10(smax/smin)*log10(dble(idust))/log10(dble(ndust))+log10(smin))
  enddo
  do idust=1,ndust
     sdust(idust)= sqrt(sdustbords(idust+1)* sdustbords(idust))
  enddo
end subroutine size_dust
