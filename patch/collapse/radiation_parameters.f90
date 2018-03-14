module radiation_parameters
  use amr_parameters
  use hydro_parameters,only:ngrp,nvar_bicg,nvar_trad

  ! DT adaptatif
  real(dp),allocatable,dimension(:,:)::rad_flux ! Flux entrant dans une cellule
  real(dp),allocatable,dimension(:,:)::urad     ! Old values of Erg in NR iterations
  real(dp),allocatable,dimension(:,:)::frad     ! Radiative force  
  real(dp)::Tr_floor=10.0 ! Background radiation field temperature - WARNING: it affects the pressure_fix in set_uold.
  integer::ntp,nfr

  real(dp):: alpha_imp = 1.0d0	!0.0:explicite 0.5:CN 1.0:implicite
  real(dp):: robin = 1.0d0	!0.0:Von Neumann 1.0:Dirichlet

  ! Multigroup
  integer,parameter::Nomega=100     ! Number of points in the omega data to compute Q moment term

  real(dp),parameter:: aR=7.56591469318689378e-015_dp
  real(dp),parameter::Tray_min=0.5d0 ! Minimum temperature in the radiative energy
  real(dp),parameter:: eray_min=(aR)*Tray_min**4 ! minimum rad energy inside frequency group
  real(dp),parameter:: deray_min=(4.0d0*aR)*Tray_min**3 ! minimum rad energy derivative inside frequency group
  real(dp):: small_er=1.0d-30       ! minimum rad energy inside frequency group in code units
  
  real(dp) :: numin=1.0d5,numax=1.0d19 ! Overall frequency boudaries
  real(dp) :: frequency_upperlimit=1.0d35 ! High end frequency if 'extra_end_group = .true.

  integer::Ninv_art4=1000                               ! Number of points in tabulated arT4 function
  real(dp),dimension(:    ),allocatable::dEr_inv_art4   ! Radiative energy increment
  real(dp),dimension(:,:  ),allocatable::inverse_art4_T ! array for tabulated arT4 function dT regular
  real(dp),dimension(:,:,:),allocatable::inverse_art4_E ! array for tabulated arT4 function dE regular
  
  real(dp), dimension(:), allocatable :: nu_min_hz ! minimum freq of given group in Hz
  real(dp), dimension(:), allocatable :: nu_max_hz ! maximum freq of given group in Hz
  real(dp), dimension(:), allocatable :: nu_min_ev ! minimum freq of given group in eV
  real(dp), dimension(:), allocatable :: nu_max_ev ! maximum freq of given group in eV
  
  real(dp),dimension(0:Nomega):: f_array,w1_array,dw1_array,w2_array,dw2_array ! Arrays of omega terms for Q computation

  logical :: freqs_in_Hz=.true.      ! Frequency units in Hz if true; if not eV
  logical :: read_groups=.false.     ! Read group boundaries from file if true
  logical :: split_groups_log=.true. ! Automatic splitting of group in log if true; if not use regular splitting
  logical :: extra_end_group=.false. ! The last group holds frequencies numax -> frequency_upperlimit if true
  logical :: grey_rad_transfer=.true.! Default: grey radiation transfer
  logical :: external_radiation_field=.false. ! Default: No external radiation background (@ Tr_floor)
  logical :: stellar_photon=.false.  ! Stellar photons are treated as a separate group (igrp=1). No emission for this group (radiation_source=0)

  ! Opacities
  character(len=12) :: opacity_type = 'grey'  ! 'grey' or 'multigroup'

  ! Radiation solver parameters
  real(dp)::epsilon_diff=1d-6                        ! CG iteration break criteria
  character(LEN=10)::fld_limiter='nolim'             ! Flux limiter (nolim, levermore or minerbo)
  integer::i_fld_limiter
  integer,parameter::i_fld_limiter_nolim=0
  integer,parameter::i_fld_limiter_minerbo=1
  integer,parameter::i_fld_limiter_levermore=2
  integer :: niter=0                                 ! Total number of iteration
  real(dp),dimension(1:10)::dtdiff_params=1d10       ! Conduction time step behaviour
  real(dp),dimension(1:10)::rosseland_params=1.0     ! Rosseland opacity coefficient's parameters
  real(dp),dimension(1:10)::planck_params=1.0        ! Planck opacity coefficient's parameters
  real(dp)::mu_gas=2.31d0                            ! For Saumon Chabrier EOS
  real(dp)::min_optical_depth=1.d-6        ! set the minimum optical depth in the cell (it may accelerate convergence in optically thin regions)

  ! Variables needed for BICG scheme
  real(dp),dimension(:,:,:,:),allocatable :: coeff_glob_left,coeff_glob_right
  real(dp),dimension(:,:,:  ),allocatable :: var_bicg,precond_bicg
  real(dp),dimension(:,:,:  ),allocatable :: mat_residual_glob
  real(dp),dimension(:,:    ),allocatable :: residual_glob
  real(dp),dimension(:,:    ),allocatable :: kappaR_bicg
  logical::block_diagonal_precond_bicg ! if .false. only diagonal, if .true. block diagonal
  integer :: i_rho,i_beta,i_y,i_pAp,i_s
  integer , dimension(1:nvar_bicg) :: ind_bicg
  real(dp), dimension(1:nvar_bicg) :: norm_bicg
  integer , dimension(1:nvar_trad) :: ind_trad
  real(dp), dimension(1:nvar_trad) :: norm_trad
  logical , dimension(1:nvar_trad) :: is_radiative_energy

  integer                                   :: irad_trans_model        !< Integer designating radiative transfer model: 0 = P1, 1 = M1
  integer, parameter                        :: irad_trans_model_p1 = 0 !< P1 radiative transfer model identifier
  integer, parameter                        :: irad_trans_model_m1 = 1 !< M1 radiative transfer model identifier
  integer                                   :: n_points                !< Number of points in the tabulated eigenvalues curve
  real(dp), dimension(:,:,:,:), allocatable :: valp                    !< Array to hold the tabulated eigenvalues as a function of \f$\theta\f$ and \f$\epsilon\f$
  real(dp)                                  :: valp_min=0.0_dp

  logical::store_matrix=.true.

#if USE_FLD==1 && NGRP == 1
  logical, parameter :: bicg_to_cg = .true.
#else
  logical, parameter :: bicg_to_cg = .false.
#endif

end module radiation_parameters

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Module MOD_OPACITIES:
!
!>  Module for opacities
!<
module mod_opacities
  use amr_parameters, only : dp
  implicit none
  real(dp), dimension(:,:,:,:), allocatable :: kappa_opmesh_p      !< Regular mesh of Planck opacities
  real(dp), dimension(:,:,:,:), allocatable :: kappa_opmesh_r      !< Regular mesh of Rosseland opacities
  real(dp), dimension(:,:,:  ), allocatable :: kappa_dustgas_p     !< Rosseland opacities from Dust + Gas
  real(dp), dimension(:,:,:  ), allocatable :: kappa_dustgas_r     !< Rosseland opacities from Dust + Gas
  real(dp), dimension(:      ), allocatable :: logt_dustgas        !< Temperature of points for Dust+Gas data
  real(dp), dimension(:      ), allocatable :: logd_dustgas        !< Density of points for Dust+Gas data
  real(dp), dimension(:      ), allocatable :: x_opmesh            !< X (density) coordinates of mesh points
  real(dp), dimension(:      ), allocatable :: y_opmesh            !< Y (temperature) coordinates of mesh points
  real(dp), dimension(:      ), allocatable :: z_opmesh            !< Z (rad temperature) coordinates of mesh points
  real(dp), dimension(:      ), allocatable :: numin_dustgas       !< Minimum freq. in Dust+Gas opacities
  real(dp), dimension(:      ), allocatable :: numax_dustgas       !< Maximum freq. in Dust+Gas opacities
  real(dp), dimension(:      ), allocatable :: extrapol_min        !< Extrapolation power for low frequencies
  real(dp), dimension(:      ), allocatable :: extrapol_max        !< Extrapolation power for high frequencies
  integer , dimension(:      ), allocatable :: nfreq_dustgas       !< Number of freq. in Dust+Gas opacities
  integer                                   :: nx_opmesh           !< Number of X points in opacity mesh
  integer                                   :: ny_opmesh           !< Number of Y points in opacity mesh
  integer                                   :: nz_opmesh           !< Number of Z points in opacity mesh
  integer                                   :: npoints             !< Number of points in freq. for Dust+Gas opacities
  real(dp)                                  :: tmin_op             !< 
  real(dp)                                  :: dx_opmesh           !< dx in opacity mesh
  real(dp)                                  :: dy_opmesh           !< dy in opacity mesh
  real(dp)                                  :: dz_opmesh           !< dy in opacity mesh
  real(dp)                                  :: opacity_dtemp
  real(dp), parameter                       :: t_grain_sublim = 1.5e+03_dp !< Grain sublimation temperature
end module mod_opacities

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Subroutine INIT_OPACITIES_DUST_AND_GAS:
!
!> Read opacities for dust from Semenov+Draine and gas from Franck.
!!    - The arbitrary points in Rho and T are read.
!!    - The Planck and Rosseland means are computed for each point for each group.
!!    - A Delaunay triangulation is computed from the points.
!!    - Each triangle defines a plane in 3D.
!!    - Then a finer but regular mesh in (Rho,T) is overlayed onto the triangulation.
!!    - This gives a regular mesh of Planck and Rosseland means to perform bicubic interpolations.
!<
subroutine init_opacities

  use radiation_parameters
  use mod_opacities
  use amr_commons, only : myid,nrestart

  implicit none

  integer                                 :: i,igroup,ipoint,j,it,nnumax,k,l,kw,kernel_size,istep
  integer                                 :: inu_min1,inu_max1,inu_min2,inu_max2,nparts,iprog,itrad
  integer                                 :: datasize,dataread,percentage,nneighbours,in,ii,jj,kk,np
  integer                                 :: iter,itermax,ntot,nmin,maxneighbours,nn,n,npasses,ntrad
  logical                                 :: enough_points_found,hole_filled
  real(dp), dimension(:    ), allocatable :: opnu,opknu
  real(dp), dimension(:,:,:), allocatable :: bin_count
  real(dp)                                :: integral1,integral2,integral3,integral4,integral5,integral6
  real(dp)                                :: dist,distmin,temp,x,y,xa,ya,xx,yy,dtrad
  real(dp)                                :: dmax_opmesh,dmin_opmesh,tmax_opmesh,tmin_opmesh,trmax_opmesh,trmin_opmesh
  real(dp)                                :: nu1,nu2,op1,op2,m,slope,dnu,kappa_min,grad1,grad2,minmod
  character (len=200)                     :: opfilename,fname
  logical , dimension(:,:,:), allocatable :: i_am_a_hole

  if(opacity_type == 'multigroup')then

     if(nrestart .eq. 0)then
        opfilename = 'VaytetEtAl2013AandA557A90_opacities.bin'

        if(myid==1)then
           write(*,*)
           write(*,*) '############# MULTIGROUP DUST AND GAS OPACITIES ##############'
           write(*,*) 'Reading opacity table: '//trim(opfilename)
           write(*,*) '=============================================================='
           write(*,*) 'Opacities: Vaytet et al. 2013, A&A, 557, A90'
           write(*,*) 'COMPUTING Planck and Rosseland mean opacities for:'
           write(*,*) ' - DUST (Semenov et al. 2003, Draine 2003) : 5K < T < 1500K'
           write(*,*) ' - MOLECULAR GAS (Ferguson et al. 2005): 1500K < T < 3200K'
           write(*,*) ' - ATOMIC GAS (Badnell et al. 2005): 3200K < T < 1.0e8K'
           write(*,*) '=============================================================='
        endif

        !####### Set variables #######
        nx_opmesh      = 100
        ny_opmesh      = 100
        nz_opmesh      = 100
        dmin_opmesh    = -24.00_dp
        dmax_opmesh    =   6.00_dp
        tmin_opmesh    =   0.00_dp
        tmax_opmesh    =   8.00_dp
        trmin_opmesh   =   0.00_dp
        trmax_opmesh   =   8.00_dp
        nparts         = 10
        ntrad          = 150
        kappa_min      = 1.0e-50_dp
        datasize       = 31637055
        !#############################

        dtrad = (trmax_opmesh-trmin_opmesh)/real(ntrad-1,dp)

        open(78,file=opfilename,form='unformatted')
        read(78) npoints,nnumax

        allocate(opnu(nnumax),opknu(nnumax))
        allocate(logt_dustgas(npoints),logd_dustgas(npoints))
        allocate(numin_dustgas(npoints),numax_dustgas(npoints),nfreq_dustgas(npoints))
        allocate(extrapol_min(npoints),extrapol_max(npoints))
        allocate(kappa_dustgas_p(ngrp,npoints,ntrad),kappa_dustgas_r(ngrp,npoints,ntrad))

        istep = 1
        iprog = 0
        dataread = 0

        do ipoint = 1,npoints

           percentage = nint(real(dataread)*100.0/real(datasize))

           if((myid==1) .and. (percentage .ge. iprog*istep))then
              write(*,'(i3,a)') percentage,'% complete'
              iprog = iprog + 1
           endif

           read(78) logd_dustgas(ipoint),logt_dustgas(ipoint)
           read(78) nfreq_dustgas(ipoint),numin_dustgas(ipoint),numax_dustgas(ipoint)
           read(78) extrapol_min(ipoint),extrapol_max(ipoint)
           read(78) opnu (1:nfreq_dustgas(ipoint))
           read(78) opknu(1:nfreq_dustgas(ipoint))

           dataread = dataread + nfreq_dustgas(ipoint)

           ! Begin loop over Trad
           do itrad = 1,ntrad

              temp = 10.0_dp**(real(itrad-1,dp)*dtrad + trmin_opmesh)

              ! find group boundaries
              inu_min1 = 0 ; inu_max1 = 0
              inu_min2 = 0 ; inu_max2 = 0

              do igroup = 1,ngrp

                 if(nu_min_hz(igroup) .lt. opnu(1))then
                    inu_min1 = 1 ; inu_min2 = 1
                 elseif(nu_min_hz(igroup) .ge. opnu(nfreq_dustgas(ipoint)))then
                    inu_min1 = nfreq_dustgas(ipoint) ; inu_min2 = nfreq_dustgas(ipoint)
                 else
                    do i = 1,nfreq_dustgas(ipoint)-1
                       if((nu_min_hz(igroup) .ge. opnu(i)).and.(nu_min_hz(igroup) .lt. opnu(i+1)))then
                          inu_min1 = i ; inu_min2 = i+1
                          exit
                       endif
                    enddo
                 endif

                 if(nu_max_hz(igroup) .lt. opnu(1))then
                    inu_max1 = 1 ; inu_max2 = 1
                 elseif(nu_max_hz(igroup) .ge. opnu(nfreq_dustgas(ipoint)))then
                    inu_max1 = nfreq_dustgas(ipoint) ; inu_max2 = nfreq_dustgas(ipoint)
                 else
                    do i = 1,nfreq_dustgas(ipoint)-1
                       if((nu_max_hz(igroup) .ge. opnu(i)).and.(nu_max_hz(igroup) .lt. opnu(i+1)))then
                          inu_max1 = i ; inu_max2 = i+1
                          exit
                       endif
                    enddo
                 endif

                 ! compute Planck and Rosseland mean opacities

                 integral1 = 0.0_dp
                 integral2 = 0.0_dp
                 integral3 = 0.0_dp
                 integral4 = 0.0_dp
                 integral5 = 0.0_dp
                 integral6 = 0.0_dp

                 ! If first frequency is outside opacity frequency range
                 if(inu_min1 .eq. inu_min2)then
                    ! Case where both frequencies are outside table on the same side
                    if(inu_min1 .eq. inu_max1)then
                       ! Select which slope to use
                       if(inu_min1 .eq. 1)then
                          slope = extrapol_min(ipoint)
                       else
                          slope = extrapol_max(ipoint)
                       endif
                       ! Split interval into 10 parts
                       dnu = (log10(nu_max_hz(igroup)) - log10(nu_min_hz(igroup)))/real(nparts,dp)
                       do i = 1,nparts
                          nu1 = 10.0_dp**(real(i-1,dp)*dnu + log10(nu_min_hz(igroup)))
                          nu2 = 10.0_dp**(real(i  ,dp)*dnu + log10(nu_min_hz(igroup)))
                          op1 = 10.0_dp**(slope * (log10(nu1)-log10(opnu(inu_min1))) + log10(opknu(inu_min1)))
                          op2 = 10.0_dp**(slope * (log10(nu2)-log10(opnu(inu_min1))) + log10(opknu(inu_min1)))
                          op1 = max(op1,kappa_min)
                          op2 = max(op2,kappa_min)
                          call compute_integral(integral1,integral2,integral3,integral4, &
                               integral5,integral6,nu1,nu2,op1,op2,temp)
                       enddo
                    else
                       !Case where only numin is outside of range
                       ! Select which slope to use
                       if(inu_min1 .eq. 1)then
                          slope = extrapol_min(ipoint)
                       else
                          slope = extrapol_max(ipoint)
                       endif
                       ! Split interval into 10 parts
                       dnu = (log10(opnu(inu_min1))-log10(nu_min_hz(igroup)))/real(nparts,dp)
                       do i = 1,nparts
                          nu1 = 10.0_dp**(real(i-1,dp)*dnu + log10(nu_min_hz(igroup)))
                          nu2 = 10.0_dp**(real(i  ,dp)*dnu + log10(nu_min_hz(igroup)))
                          op1 = 10.0_dp**(slope * (log10(nu1)-log10(opnu(inu_min1))) + log10(opknu(inu_min1)))
                          op2 = 10.0_dp**(slope * (log10(nu2)-log10(opnu(inu_min1))) + log10(opknu(inu_min1)))
                          op1 = max(op1,kappa_min)
                          op2 = max(op2,kappa_min)
                          call compute_integral(integral1,integral2,integral3,integral4, &
                               integral5,integral6,nu1,nu2,op1,op2,temp)
                       enddo
                    endif
                 endif

                 ! If last frequency is outside opacity frequency range
                 if((inu_max1 .eq. nfreq_dustgas(ipoint)) .and. (inu_max1 .eq. inu_max1))then
                    !Case where numax is outside of range
                    slope = extrapol_max(ipoint)
                    ! Split interval into 10 parts
                    dnu = (log10(nu_max_hz(igroup))-log10(opnu(inu_max1)))/real(nparts,dp)
                    do i = 1,nparts
                       nu1 = 10.0_dp**(real(i-1,dp)*dnu + log10(opnu(inu_max1)))
                       nu2 = 10.0_dp**(real(i  ,dp)*dnu + log10(opnu(inu_max1)))
                       op1 = 10.0_dp**(slope * (log10(nu1)-log10(opnu(inu_max1))) + log10(opknu(inu_max1)))
                       op2 = 10.0_dp**(slope * (log10(nu2)-log10(opnu(inu_max1))) + log10(opknu(inu_max1)))
                       op1 = max(op1,kappa_min)
                       op2 = max(op2,kappa_min)
                       call compute_integral(integral1,integral2,integral3,integral4, &
                            integral5,integral6,nu1,nu2,op1,op2,temp)
                    enddo
                 endif

                 ! First part of the curve between numin and inumin
                 if(inu_min1 .ne. inu_min2)then
                    nu1 = opnu (inu_min1)
                    nu2 = opnu (inu_min2)
                    m = ( opknu(inu_min2) - opknu(inu_min1) ) / (nu2 - nu1)
                    op1 = m * (nu_min_hz(igroup)-nu1) + opknu(inu_min1)
                    op2 = opknu(inu_min2)
                    nu1 = nu_min_hz(igroup)
                    call compute_integral(integral1,integral2,integral3,integral4, &
                         integral5,integral6,nu1,nu2,op1,op2,temp)
                 endif

                 ! Last part of the curve between inumax and numax
                 if(inu_max1 .ne. inu_max2)then
                    nu1 = opnu (inu_max1)
                    nu2 = opnu (inu_max2)
                    m = ( opknu(inu_max2) - opknu(inu_max1) ) / (nu2 - nu1)
                    op1 = opknu(inu_max1)
                    op2 = m * (nu_max_hz(igroup)-nu1) + opknu(inu_max1)
                    nu2 = nu_max_hz(igroup)
                    call compute_integral(integral1,integral2,integral3,integral4, &
                         integral5,integral6,nu1,nu2,op1,op2,temp)
                 endif

                 ! Middle part between inumin and inumax
                 do i = inu_min2,inu_max1-1
                    nu1 = opnu (i  )
                    nu2 = opnu (i+1)
                    op1 = opknu(i  )
                    op2 = opknu(i+1)
                    call compute_integral(integral1,integral2,integral3,integral4, &
                         integral5,integral6,nu1,nu2,op1,op2,temp)
                 enddo

                 if(abs(integral1) < 1.0e-50_dp)then
                    kappa_dustgas_p(igroup,ipoint,itrad) = integral5 / integral6
                 else
                    kappa_dustgas_p(igroup,ipoint,itrad) = integral1 / integral2
                 endif
                 if(abs(integral3) < 1.0e-50_dp)then
                    kappa_dustgas_r(igroup,ipoint,itrad) = integral5 / integral6
                 else
                    kappa_dustgas_r(igroup,ipoint,itrad) = integral3 / integral4
                 endif

                 !         if(logt_dustgas(ipoint) > 5.0)then
                 !            write(*,*) ipoint,igroup,integral1,integral2,integral3,integral4,&
                 !                       kappa_dustgas_p(igroup,ipoint),kappa_dustgas_r(igroup,ipoint)
                 !         endif

              enddo ! end do igroup = 1,ngrp

           enddo ! end do itrad = 1,ntrad

        enddo ! end do ipoint = 1,npoints

        close(78)

        if(myid==1)then
           write(*,*) 'Number of points:',npoints
           write(*,*) 'Rhomin,Rhomax:',minval(logd_dustgas),maxval(logd_dustgas)
           write(*,*) 'Tmin,Tmax:',minval(logt_dustgas),maxval(logt_dustgas)
        endif

        ! Create a mesh of (rho,T) points:
        if(myid==1) write(*,*) 'Computing regular mesh of opacities'

        allocate(bin_count(nx_opmesh,ny_opmesh,nz_opmesh),i_am_a_hole(nx_opmesh,ny_opmesh,nz_opmesh))
        allocate(kappa_opmesh_p(ngrp,nx_opmesh,ny_opmesh,nz_opmesh),kappa_opmesh_r(ngrp,nx_opmesh,ny_opmesh,nz_opmesh))
        allocate(x_opmesh(nx_opmesh),y_opmesh(ny_opmesh),z_opmesh(nz_opmesh))

        dx_opmesh = (dmax_opmesh-dmin_opmesh)/real(nx_opmesh-1,dp)
        dy_opmesh = (tmax_opmesh-tmin_opmesh)/real(ny_opmesh-1,dp)
        dz_opmesh = (trmax_opmesh-trmin_opmesh)/real(nz_opmesh-1,dp)     

        x_opmesh(1) = dmin_opmesh
        y_opmesh(1) = tmin_opmesh
        z_opmesh(1) = trmin_opmesh

        do i = 2,nx_opmesh
           x_opmesh(i) = x_opmesh(i-1) + dx_opmesh
        enddo
        do j = 2,ny_opmesh
           y_opmesh(j) = y_opmesh(j-1) + dy_opmesh
        enddo
        do k = 2,nz_opmesh
           z_opmesh(k) = z_opmesh(k-1) + dz_opmesh
        enddo

        bin_count = 0.0_dp
        kappa_opmesh_p = 0.0_dp
        kappa_opmesh_r = 0.0_dp
        i_am_a_hole = .true.

        do ipoint = 1,npoints

           i = floor((logd_dustgas(ipoint)-dmin_opmesh)/dx_opmesh - 0.5_dp) + 2
           j = floor((logt_dustgas(ipoint)-tmin_opmesh)/dy_opmesh - 0.5_dp) + 2

           ! now search along Trad direction
           do itrad = 1,ntrad
              temp = real(itrad-1,dp)*dtrad + trmin_opmesh

              k = floor((temp-tmin_opmesh)/dz_opmesh - 0.5_dp) + 2

              do igroup = 1,ngrp
                 kappa_opmesh_p(igroup,i,j,k) = kappa_opmesh_p(igroup,i,j,k) + kappa_dustgas_p(igroup,ipoint,itrad)
                 kappa_opmesh_r(igroup,i,j,k) = kappa_opmesh_r(igroup,i,j,k) + kappa_dustgas_r(igroup,ipoint,itrad)
              enddo

              bin_count(i,j,k) = bin_count(i,j,k) + 1.0_dp
              i_am_a_hole(i,j,k) = .false.

           enddo

        enddo

        do k = 1,nz_opmesh
           do j = 1,ny_opmesh
              do i = 1,nx_opmesh
                 if(i_am_a_hole(i,j,k))then
                    do igroup = 1,ngrp
                       kappa_opmesh_p(igroup,i,j,k) = 1.0e-30_dp
                       kappa_opmesh_r(igroup,i,j,k) = 1.0e-30_dp
                    enddo
                 else
                    do igroup = 1,ngrp
                       kappa_opmesh_p(igroup,i,j,k) = kappa_opmesh_p(igroup,i,j,k) / bin_count(i,j,k)
                       kappa_opmesh_r(igroup,i,j,k) = kappa_opmesh_r(igroup,i,j,k) / bin_count(i,j,k)
                    enddo
                 endif
              enddo
           enddo
        enddo

        ! Convert to log
        kappa_opmesh_p = log10(kappa_opmesh_p)
        kappa_opmesh_r = log10(kappa_opmesh_r)

        ! Now fill holes in table
        if(myid==1)write(*,*) 'Filling holes in table...'

        ! First find holes which are surrounded by many points, and then go to bigger ang bigger holes

        npasses = 10 ! number of times to perform filling process

        do n = 1,npasses

           nneighbours = 26
           maxneighbours = 0
           itermax = 8
           do in = nneighbours,0,-1

              do k = 1,nz_opmesh
                 do j = 1,ny_opmesh
                    do i = 1,nx_opmesh
                       if(i_am_a_hole(i,j,k))then

                          ! Find number of neighbours
                          nn = 0
                          do ii = max(i-1,1),min(i+1,nx_opmesh)
                             do jj = max(j-1,1),min(j+1,ny_opmesh)
                                do kk = max(k-1,1),min(k+1,nz_opmesh)
                                   if(.not.i_am_a_hole(ii,jj,kk))then
                                      nn = nn + 1
                                   endif
                                enddo
                             enddo
                          enddo

                          if(nn == in)then

                             hole_filled = .false.

                             do igroup = 1,ngrp

                                do iter = 1,itermax

                                   np = 0
                                   xx = 0.0_dp
                                   yy = 0.0_dp
                                   ntot = (1 + 2*iter)**3 - 1
                                   nmin = int(2.0_dp*sqrt(real(ntot,dp)))   ! ntot / 2 + 1
                                   enough_points_found = .false.

                                   do ii = max(i-iter,1),min(i+iter,nx_opmesh)
                                      do jj = max(j-iter,1),min(j+iter,ny_opmesh)
                                         do kk = max(k-iter,1),min(k+iter,nz_opmesh)

                                            if(.not.i_am_a_hole(ii,jj,kk))then
                                               xx = xx + kappa_opmesh_p(igroup,ii,jj,kk)
                                               yy = yy + kappa_opmesh_r(igroup,ii,jj,kk)
                                               np = np + 1
                                            endif

                                         enddo
                                      enddo
                                   enddo

                                   if(np .ge. nmin)then
                                      enough_points_found = .true.
                                      exit
                                   endif

                                enddo

                                if(enough_points_found)then
                                   kappa_opmesh_p(igroup,i,j,k) = xx / real(np,dp)
                                   kappa_opmesh_r(igroup,i,j,k) = yy / real(np,dp)
                                   hole_filled = .true.
                                endif

                             enddo

                             if(hole_filled) i_am_a_hole(i,j,k) = .false.

                          endif

                       endif

                    enddo ! end do i = 1,nx_opmesh
                 enddo ! end do j = 1,ny_opmesh
              enddo ! end do k = 1,nz_opmesh

           enddo ! end do in = nneighbours,0,-1

        enddo ! end do n = 1,npasses

        if(myid==1)write(*,*) 'Done'

        ! Free memory
        deallocate(opnu,opknu)
        deallocate(logt_dustgas,logd_dustgas)
        deallocate(numin_dustgas,numax_dustgas,nfreq_dustgas)
        deallocate(extrapol_min,extrapol_max)
        deallocate(kappa_dustgas_p,kappa_dustgas_r)
        deallocate(bin_count,i_am_a_hole)

        if(myid==1)then
           open (79,file='multigroup_opacity.bin',form='unformatted')
           write (79) nx_opmesh,ny_opmesh,nz_opmesh,dx_opmesh,dy_opmesh,dz_opmesh,dmin_opmesh,dmax_opmesh,tmin_opmesh,tmax_opmesh,trmin_opmesh,trmax_opmesh

           write (79) x_opmesh
           write (79) y_opmesh
           write (79) z_opmesh
           write (79) kappa_opmesh_p(1:ngrp,1:nx_opmesh,1:ny_opmesh,1:nz_opmesh)
           write (79) kappa_opmesh_r(1:ngrp,1:nx_opmesh,1:ny_opmesh,1:nz_opmesh)
           close(79)
        end if
     else
        if(myid==1)write(*,*) 'Reading opacity table from previous run'

        open (80,file='multigroup_opacity.bin',form='unformatted')
        read (80) nx_opmesh,ny_opmesh,nz_opmesh,dx_opmesh,dy_opmesh,dz_opmesh,dmin_opmesh,dmax_opmesh,tmin_opmesh,tmax_opmesh,trmin_opmesh,trmax_opmesh

        allocate(x_opmesh(nx_opmesh),y_opmesh(ny_opmesh),z_opmesh(nz_opmesh))
        allocate(kappa_opmesh_p(ngrp,nx_opmesh,ny_opmesh,nz_opmesh),kappa_opmesh_r(ngrp,nx_opmesh,ny_opmesh,nz_opmesh))

        read (80) x_opmesh
        read (80) y_opmesh
        read (80) z_opmesh
        read (80) kappa_opmesh_p(1:ngrp,1:nx_opmesh,1:ny_opmesh,1:nz_opmesh)
        read (80) kappa_opmesh_r(1:ngrp,1:nx_opmesh,1:ny_opmesh,1:nz_opmesh)
        close(80)

     end if

  else

     if(ngrp .gt. 1 .and. myid==1)then
        write(*,*) 
        write(*,*) '=============================================================='
        write(*,*) 'WARNING: using grey opacities with NGROUP>1!'
        write(*,*) '=============================================================='
     end if

     opfilename = 'vaytet_grey_opacities3D.bin'

     if(myid==1)then
        write(*,*)
        write(*,*) '################ GREY DUST AND GAS OPACITIES #################'
        write(*,*) 'Reading opacity table: '//trim(opfilename)
        write(*,*) '=============================================================='
        write(*,*) 'Opacities: Vaytet et al. 2013, A&A, 557, A90'
        write(*,*) 'READING Planck and Rosseland mean opacities for:'
        write(*,*) ' - DUST (Semenov et al. 2003, Draine 2003) : 5K < T < 1500K'
        write(*,*) ' - MOLECULAR GAS (Ferguson et al. 2005): 1500K < T < 3200K'
        write(*,*) ' - ATOMIC GAS (Badnell et al. 2005): 3200K < T < 1.0e8K'
        write(*,*) '=============================================================='
     endif

     open (78,file=trim(opfilename),form='unformatted')
     read (78) nx_opmesh,ny_opmesh,nz_opmesh,dx_opmesh,dy_opmesh,dz_opmesh,dmin_opmesh,dmax_opmesh,tmin_opmesh,tmax_opmesh,trmin_opmesh,trmax_opmesh

     allocate(x_opmesh(nx_opmesh),y_opmesh(ny_opmesh),z_opmesh(nz_opmesh))
     allocate(kappa_opmesh_p(ngrp,nx_opmesh,ny_opmesh,nz_opmesh),kappa_opmesh_r(ngrp,nx_opmesh,ny_opmesh,nz_opmesh))

     read (78) x_opmesh
     read (78) y_opmesh
     read (78) z_opmesh
     read (78) kappa_opmesh_p(1,1:nx_opmesh,1:ny_opmesh,1:nz_opmesh)
     read (78) kappa_opmesh_r(1,1:nx_opmesh,1:ny_opmesh,1:nz_opmesh)
     close(78)

     if(ngrp .gt. 1)then
        do k = 1,nz_opmesh
           do j = 1,ny_opmesh
              do i = 1,nx_opmesh
                 do igroup=2,ngrp
                    kappa_opmesh_p(igroup,i,j,k) = kappa_opmesh_p(1,i,j,k)
                    kappa_opmesh_r(igroup,i,j,k) = kappa_opmesh_r(1,i,j,k)
                 end do
              end do
           end do
        end do
     end if
  endif

  if(myid==1)then
     write(*,*) 'INIT_OPACITIES complete'
     write(*,*) '##############################################################'
  endif

  return

end subroutine init_opacities

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Function PLANCK_ANA:
!
!> Compute Planck average opacity.
!<
function planck_ana(dens,Tp,Tr,igroup)

  use radiation_parameters
  use mod_opacities
  use pm_commons,only:Teff_sink

  implicit none

  integer, intent(in) :: igroup
  real(dp),intent(in) :: dens,Tp,Tr
  
  integer             :: ival,jval,kval
  real(dp)            :: x,y,z,dx,dy,dz
  real(dp)            :: x0,x1,y0,y1,z0,z1
  real(dp)            :: c00,c01,c10,c11,c0,c1
  real(dp)            :: planck_ana
  
  ! compute dust and gas opacities
  x = log10(dens)
  y = log10(Tp)
  z = log10(Tr)
  if(stellar_photon)then
     if(igroup==1 .and. maxval(Teff_sink).gt.0)z = log10(maxval(Teff_sink)) 
  end if
    
  ival = floor((x - x_opmesh(1)) / dx_opmesh) + 1
  jval = floor((y - y_opmesh(1)) / dy_opmesh) + 1
  kval = floor((z - z_opmesh(1)) / dz_opmesh) + 1
  
  ! enforce to be in the table
  ival = min(max(ival,1),nx_opmesh)
  jval = min(max(jval,1),ny_opmesh)
  kval = min(max(kval,1),nz_opmesh)

  ! Perform tri-linear interpolation
  
  ! Compute coordinate deltas
  x0 = x_opmesh(ival  )
  x1 = x_opmesh(ival+1)
  y0 = y_opmesh(jval  )
  y1 = y_opmesh(jval+1)
  z0 = z_opmesh(kval  )
  z1 = z_opmesh(kval+1)
  
  dx = (x-x0)/(x1-x0)
  dy = (y-y0)/(y1-y0)
  dz = (z-z0)/(z1-z0)
  
  ! First linear interpolation along x
  c00 = kappa_opmesh_p(igroup,ival  ,jval  ,kval  )*(1.0_dp-dx) + kappa_opmesh_p(igroup,ival+1,jval  ,kval  )*dx
  c10 = kappa_opmesh_p(igroup,ival  ,jval+1,kval  )*(1.0_dp-dx) + kappa_opmesh_p(igroup,ival+1,jval+1,kval  )*dx
  c01 = kappa_opmesh_p(igroup,ival  ,jval  ,kval+1)*(1.0_dp-dx) + kappa_opmesh_p(igroup,ival+1,jval  ,kval+1)*dx
  c11 = kappa_opmesh_p(igroup,ival  ,jval+1,kval+1)*(1.0_dp-dx) + kappa_opmesh_p(igroup,ival+1,jval+1,kval+1)*dx
  
  ! Second linear interpolation along y
  c0 = c00*(1.0_dp-dy) + c10*dy
  c1 = c01*(1.0_dp-dy) + c11*dy
  
  ! Third linear interpolation along z
  planck_ana = c0*(1.0_dp-dz) + c1*dz
  
  planck_ana = dens*10.0_dp**(planck_ana)

end function planck_ana

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Function ROSSELAND_ANA:
!
!> Compute Rosseland mean opacity.
!<
function rosseland_ana(dens,Tp,Tr,igroup)

  use radiation_parameters
  use mod_opacities
  use pm_commons,only:Teff_sink

  implicit none

  integer, intent(in) :: igroup
  real(dp),intent(in) :: dens,Tp,Tr
  
  integer             :: ival,jval,kval
  real(dp)            :: x,y,z,dx,dy,dz
  real(dp)            :: x0,x1,y0,y1,z0,z1
  real(dp)            :: c00,c01,c10,c11,c0,c1
  real(dp)            :: rosseland_ana
  
  ! compute dust and gas opacities
  x = log10(dens)
  y = log10(Tp)
  z = log10(Tr)
  if(stellar_photon)then
     if(igroup==1 .and. maxval(Teff_sink).gt.0)z = log10(maxval(Teff_sink)) 
  end if

  ival = floor((x - x_opmesh(1)) / dx_opmesh) + 1
  jval = floor((y - y_opmesh(1)) / dy_opmesh) + 1
  kval = floor((z - z_opmesh(1)) / dz_opmesh) + 1
  
  ! enforce to be in the table
  ival = min(max(ival,1),nx_opmesh)
  jval = min(max(jval,1),ny_opmesh)
  kval = min(max(kval,1),nz_opmesh)

  ! Perform tri-linear interpolation
  
  ! Compute coordinate deltas
  x0 = x_opmesh(ival  )
  x1 = x_opmesh(ival+1)
  y0 = y_opmesh(jval  )
  y1 = y_opmesh(jval+1)
  z0 = z_opmesh(kval  )
  z1 = z_opmesh(kval+1)
  
  dx = (x-x0)/(x1-x0)
  dy = (y-y0)/(y1-y0)
  dz = (z-z0)/(z1-z0)
  
  ! First linear interpolation along x
  c00 = kappa_opmesh_r(igroup,ival  ,jval  ,kval  )*(1.0_dp-dx) + kappa_opmesh_r(igroup,ival+1,jval  ,kval  )*dx
  c10 = kappa_opmesh_r(igroup,ival  ,jval+1,kval  )*(1.0_dp-dx) + kappa_opmesh_r(igroup,ival+1,jval+1,kval  )*dx
  c01 = kappa_opmesh_r(igroup,ival  ,jval  ,kval+1)*(1.0_dp-dx) + kappa_opmesh_r(igroup,ival+1,jval  ,kval+1)*dx
  c11 = kappa_opmesh_r(igroup,ival  ,jval+1,kval+1)*(1.0_dp-dx) + kappa_opmesh_r(igroup,ival+1,jval+1,kval+1)*dx
  
  ! Second linear interpolation along y
  c0 = c00*(1.0_dp-dy) + c10*dy
  c1 = c01*(1.0_dp-dy) + c11*dy
  
  ! Third linear interpolation along z
  rosseland_ana = c0*(1.0_dp-dz) + c1*dz
  
  rosseland_ana = dens*10.0_dp**(rosseland_ana)

end function rosseland_ana

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Function SCATTERING_ANA:
!
!> This routine computes the scattering opacity kappa_s*rho
!! as a function of density and temperature.
!! Units are supposed to be in cgs here (as in units.f90)
!<
function scattering_ana(dens,Tp,Tr,igroup)

  use amr_commons
  use const

  implicit none

  integer ,intent(in)    :: igroup
  real(dp),intent(in)    :: dens,Tp,Tr
  real(dp)               :: scattering_ana
  
  scattering_ana = zero

end function scattering_ana

!##################################################################################################
!##################################################################################################
!##################################################################################################
!##################################################################################################

!  Subroutine COMPUTE_INTEGRAL:
!
!> Computes integrals for Planck and Rosseland means between nu1 and nu2.
!<
subroutine compute_integral(int1,int2,int3,int4,int5,int6,nu1,nu2,op1,op2,temp)

  use amr_parameters, only : dp

  implicit none
  
  real(dp), intent(inout) :: int1,int2,int3,int4,int5,int6
  real(dp), intent(in   ) :: nu1,nu2,op1,op2,temp
  real(dp)                :: BPlanck,Div_BPlanck

  int1 = int1 + 0.5_dp * (nu2-nu1) * (     BPlanck(nu1,temp) * op1 + &
                                           BPlanck(nu2,temp) * op2 )
  int2 = int2 + 0.5_dp * (nu2-nu1) * (     BPlanck(nu1,temp)       + &
                                           BPlanck(nu2,temp)       )
  int3 = int3 + 0.5_dp * (nu2-nu1) * ( Div_BPlanck(nu1,temp)       + &
                                       Div_BPlanck(nu2,temp)       )
  int4 = int4 + 0.5_dp * (nu2-nu1) * ( Div_BPlanck(nu1,temp) / op1 + &
                                       Div_BPlanck(nu2,temp) / op2 )
  int5 = int5 + 0.5_dp * (nu2-nu1) * ( op1 + op2 )
  int6 = int6 +          (nu2-nu1)

  return

end subroutine compute_integral
