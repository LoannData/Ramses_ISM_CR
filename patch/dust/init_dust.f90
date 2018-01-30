subroutine init_dust_ratio(dustratio,epsilondust)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp) :: dustratio
  real(dp), dimension(1:ndust):: epsilondust
  real(dp), dimension(1:ndust):: sdust
  real(dp) :: epsilon_0,Anorm,Bnorm
  integer  :: idust
  epsilon_0 = dustratio/(1.0d0+dustratio)
  Anorm =  (1.0d0-mrn_index)/(size_max**(1.0d0-mrn_index)-size_min**(1.0d0-mrn_index))
  Bnorm =(1.0D0-mrn_index)/DBLE(ndust)/Anorm
  sdust(1)=size_min
  !We bin the distribution
  do idust =1,ndust-1
     sdust(idust+1) = (sdust(idust)**(1.0d0-mrn_index)+Bnorm)**(1.0/(1.0-mrn_index))
  enddo
  ! We initialise it
  do idust=1,ndust
     epsilondust(idust)= epsilon_0/DBLE(Ndust)
  enddo
end subroutine init_dust_ratio

subroutine size_dust(sdust,dustratio)
  use amr_commons
  use hydro_commons
  implicit none
  real(dp), dimension(1:ndust):: sdust
  real(dp) :: dustratio
  real(dp) :: Anorm, Bnorm, epsilon_0 
  integer  :: idust

  epsilon_0 = dustratio/(1.0d0+dustratio)
  !We bin the distribution
  Anorm = epsilon_0 * (1.0d0-mrn_index)/(size_max**(1.0d0-mrn_index)-size_min**(1.0d0-mrn_index))
  Bnorm =(1.0-mrn_index)*epsilon_0/DBLE(ndust)/Anorm
  sdust(1)=size_min
  !We bin the distribution
  do idust =1,ndust-1
     sdust(idust+1) = (sdust(idust)**(1.0d0-mrn_index)+Bnorm)**(1.0/(1.0-mrn_index))
  enddo
 
end subroutine size_dust


