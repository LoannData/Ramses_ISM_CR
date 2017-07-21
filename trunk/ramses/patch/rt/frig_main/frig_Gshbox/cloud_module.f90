module cloud_module
  use amr_parameters
  use hydro_parameters,only:Msun
  use rt_parameters

  ! TODO - CLEAN THIS OUT

  !initial temperature used for the isothermal run
  real(dp)::temper
  real(dp)::temper_iso

  !feedback from jet
  logical:: jet = .false., rad_jet=.false. 
  real(dp)::Ucoef=1.
  real(dp):: mass_jet_sink=0. !mass above which a jets is included

  !Initial conditions parameter for the dense core
  real(dp)::bx_bound=0.
  real(dp)::by_bound=0.
  real(dp)::bz_bound=0.
  real(dp)::turb=0.
  real(dp)::dens0=0.
  real(dp)::V0=0.
  real(dp)::Height0=0.


  real(dp)::bl_fac=1.   !multiply calculated boxlen by this factor


  !Initial conditions parameters for the dense core
  logical ::bb_test=.false. ! Activate Boss & Bodenheimer inital conditions instead of 1/R^2 density profile
  logical ::uniform_bmag=.false. ! Activate uniform magnetic field initial conditions for BE-like initial density profile
  real(dp)::mass_c=1.         !cloud mass in solar mass
  real(dp)::contrast=100.d0   !density contrast (used when bb_test=.true.)
  real(dp)::cont=1.           !density contrast (used when bb_test=.false.)
  real(dp)::rap=1.            !axis ratio
  real(dp)::ff_sct=1.         !freefall time / sound crossing time
  real(dp)::ff_rt=1.          !freefall time / rotation time
  real(dp)::ff_act=1.         !freefall time / Alfven crossing time
  real(dp)::ff_vct=1.         !freefall time / Vrms crossing time
  real(dp)::theta_mag=0.      !angle between magnetic field and rotation axis
  real(dp)::thet_mag=0.      !angle between magnetic field and rotation axis

  real(dp):: C2_vis=0.0d0 !Von Neumann & Richtmeyer artificial viscosity coefficient 3 en principe
  real(dp):: alpha_dense_core=0.5d0
  real(dp):: beta_dense_core=0.0d0
  real(dp):: crit=0.0d0
  real(dp):: delta_rho=0.0d0
  real(dp):: Mach=0.0d0


  real(dp)::Vshear=0.    !value of the shear in km/s
  logical:: shear=.false.        !add Corriolis and centrifuge forces



  ! PMS evolution related stuff
  logical :: rt_feedback=.false.       ! take into account RT feedback
  logical :: PMS_evol=.false.          ! Take into account PMS evolution subgrid model
  logical :: Hosokawa_track=.false.    ! Take into account PMS evolution subgrid model
  real(dp):: dt_lsink_update=50        ! frequency of the sink luminosity update with PMS evolution (in yr)
  real(dp):: epsilonlib=0.0            ! Fraction of energy absorbed by the prostostars at the accretion shock
  real(dp):: mprotostar=0.0009546*Msun ! initial mass of the protostar (1 Mjup)
  real(dp):: rstar_init=2.5            ! Initial radius of the protostar in Rsun
  integer :: modell=0
  integer :: modrestart=0              ! name of model you want to restart from, this is an input
  real(dp):: facc_star_lum=0.75d0      ! fraction of the accretion luminosity radiated by the sinks
  real(dp):: facc_star=0.5d0           ! fraction of the sink accreted mass actually accreted by the star
  integer::nmdot_PMS,nm_PMS,ndata_PMS
  integer ,allocatable,dimension(:)::nb_ligne_PMS
  real(dp),allocatable,dimension(:,:,:)::data_PMS

  ! Scale time by a factor? (e.g. t_ff)
  real(dp)::scale_tout=1d0


end module cloud_module

!================================================================
!================================================================
!================================================================
!================================================================

subroutine calc_dmin(d_c)
  use amr_commons
  use hydro_commons
  use cloud_module
  implicit none

  ! NOTE!! - IS THIS REALLY NECESSARY? - SAM GEEN OCTOBER 2015

  real(dp):: d_c, cont_ic, dmin

  cont_ic = 10.
  dmin = d_c / cont / cont_ic

  if (myid == 1) then
    write(*,*) "dmin = ", dmin
  endif
end subroutine calc_dmin
!================================================================
!================================================================
!================================================================
!================================================================
subroutine calc_boxlen
  use amr_commons
  use amr_parameters
  use hydro_commons
  use poisson_parameters
  use cloud_module
!  use const
  implicit none
  !================================================================
  !this routine calculate boxlen
  !================================================================
  integer :: i
  real(dp):: pi
  real(dp):: d_c,zeta
  real(dp):: res_int,r_0,C_s
  integer::  np
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),save:: first
  real(dp):: mu=1.4d0 ! NOTE - MUST BE THE SAME AS IN units.f90!!
!  real(dp)::myid

!   myid=1

    if (first .eq. 0.) then

    pi=acos(-1.0d0)

    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_T2 = scale_T2 * mu

    !calculate the mass in code units (Msolar / Mparticle / pc^3
    mass_c = mass_c * (2.d33 / (scale_d * scale_l**3) )

    !calculate the sound speed
    C_s = sqrt( T2_star / scale_T2)

    !calculate  zeta=r_ext/r_0
    zeta = sqrt(cont - 1.)

    !calculate an integral used to compute the cloud radius
    np=1000
    res_int=0.
    do i=1,np
     res_int = res_int + log(1.+(zeta/np*i)**2) * zeta/np
    enddo
    res_int = zeta*log(1.+zeta**2) - res_int

    !now we determine the central density and the external cloud radius
    !we have mass = 2 pi rho_c r_0^2 z_0 * res_int
    !which results from the integration of rho = dc/(1.+(x^2+y^2)/r_O^2+z^2/z_0^2)
    !for (x^2+y^2)/r_O^2+z^2/z_0^2 < zeta
    !we also have ff_sct = sqrt(3. pi / 32 / G / d_c) C_s / (r_0)
    !which just state the ratio of freefall time over sound crossing time
    !from these 2 formula, rho_c and r_0 are found to be:



    r_0 = mass_c / (2.*pi*rap*res_int) * (ff_sct)**2 / (3.*pi/32.) / C_s**2

    d_c = mass_c / (2.*pi*rap*res_int) / r_0**3

    !it is equal to twice the length of the major axis
    boxlen = r_0 * zeta * max(rap,1.) * 4.

    ! Multiply boxlen by an extra factor
    boxlen = bl_fac * boxlen

    if (myid == 1) then
    write(*,*) '** Cloud parameters estimated in calc-boxlen **'
    write(*,*) 'inner radius (pc) ', r_0
    write(*,*) 'peak density (cc) ', d_c
    write(*,*) 'total box length (pc) ', boxlen
    write(*,*) 'cloud mass (code units) ', mass_c
    write(*,*) 'boxlen (code units) ',boxlen
    write(*,*)
    endif



    first=1.
    endif

    call calc_dmin(d_c)

end subroutine calc_boxlen

subroutine read_cloud_params(nml_ok)

  use amr_parameters
  use feedback_module
  use clfind_commons
  use cloud_module

  implicit none
  logical::nml_ok
  real(dp)::cellsize
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp),parameter::pcincm=3.086d18

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------
  namelist/cloud_params/mass_c,rap,cont,ff_sct,ff_rt,ff_act,ff_vct,thet_mag &
       & ,bl_fac,switch_solv,turb,Height0,dens0,bx_bound,by_bound,bz_bound,Vshear,shear

  ! Read namelist file
  rewind(1)
  read(1,NML=cloud_params,END=101)
101 continue                                   ! No harm if no namelist

  ! Get some units out there
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

!  ! Calculate boxlen
!  if (mass_c .gt. 0) then
!     call calc_boxlen
!  end if


  write(*,*) 'shear ',shear


  ! normalise the shear in code units 
  Vshear = Vshear*1.e5/scale_v

  !since boxlen is not known initialy we must multiply the
  !refining parameters by boxlen here
  x_refine = x_refine*boxlen
  y_refine = y_refine*boxlen
  z_refine = z_refine*boxlen
  r_refine = r_refine*boxlen

  ! Set the sink formation threshold based on the Jeans criterion
!  cellsize = boxlen * 0.5**nlevelmax * pcincm / scale_l
!  n_sink = 881.0 / cellsize**2 ! Scaled to give 1e6 for 30pc/1024
!  n_clfind = 0.1 * n_sink
!  if(myid==1) write(*,*) "SETTING n_sink, n_clfind TO", n_sink, n_clfind

  ! Feedback parameters
  call read_feedback_params(nml_ok)


end subroutine read_cloud_params




!#########################################################
!#########################################################
!#########################################################
subroutine boundary_frig(ilevel)
  Use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
!  use hydro_parameters !, ONLY: nvar,boundary_var,gamma,bx_bound,by_bound,bz_bound,turb,dens0,V0
  use hydro_commons
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine set up open boundary conditions which deals properly with div B 
  ! it uses the 2 last cells of the domain
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz,j
  integer::info,ibound,nx_loc,idim,neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp),dimension(1:nvector,1:3),save::vv
  real(dp),dimension(1:nvector,1:nvar+3)::q   ! Primitive variables
  real(dp)::pi,time
  integer ::ivar,jgrid,ind_cell_vois
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,Cwnm
  real(dp)::dx_min, fact, Emag,Emag0

! STG HACK - ignore if not MHD
! TODO: Take boundary cleaner and use for non-MHD solver
#ifndef SOLVERmhd
  return
#endif 


return


  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  Cwnm = sqrt(8000./scale_T2)

  pi=ACOS(-1.0d0)

  time = t * Cwnm / boxlen



  if(numbtot(1,ilevel)==0)return

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  
  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=dble(nx_loc)/boxlen
  dx_loc=dx/scale

  dx_min = (0.5D0**levelmin)/scale

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do
  
  !-------------------------------------
  ! Compute analytical velocity field
  !-------------------------------------
  ncache=active(ilevel)%ngrid
  
  ! Loop over grids by vector sweeps
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim
        
        ! Gather cell indices
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=iskip+ind_grid(i)
        end do
        
        ! Gather cell centre positions
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              xx(i,idim)=(xx(i,idim)-skip_loc(idim))/scale
           end do
        end do
        

       do i=1,ngrid



        ! STG HACK CHECK FOR BORKED
        if (uold(ind_cell(i),1) .lt. 0d0) then
           write(*,*) "DENSITY < 0 BEFORE VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        if (uold(ind_cell(i),5) .lt. 0d0) then
           write(*,*) "TOTAL CELL ENERGY < 0 BEFORE VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        do j=5,8
           if (isnan(uold(ind_cell(i),j))) then
              write(*,*) "VARIABLE IS NAN BEFORE VELOCITY_FINE, OH NO", ind_cell(i),j,uold(ind_cell(i),1)
              call clean_stop
           end if
        end do


        !impose vanishing gradient conditions at the y  faces
        if(  xx(i,2) .lt. 2.*dx_min ) then 


             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),4))


           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 3 (1) is 4 (2)
           if(ind .eq. 3 .or. ind .eq. 4 .or. ind .eq. 7 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - 2*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 3 .or. ind .eq. 4 .or. ind .eq. 7 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+2) = uold(ind_cell_vois,7)
 
              uold(ind_cell(i),7)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6) - uold(ind_cell(i),8) 
           else
              !should be equal to uold(ind_cell(i),7) of the preceeding case 
              uold(ind_cell(i),nvar+2) =  uold(ind_cell(i),nvar+1 ) + uold(ind_cell_vois,7) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6)  - uold(ind_cell(i),8) 

              !ensure div B
              uold(ind_cell(i),7) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),6) - uold(ind_cell(i),8) 
           endif


           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif

        if(  xx(i,2) .gt. boxlen-2.*dx_min ) then 
             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),3))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 3 (4) is 1 (2)
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 5 .or. ind .eq. 6) then 
             ind_cell_vois = ind_cell_vois + 2*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 5 .or. ind .eq. 6) then 
              uold(ind_cell(i),7) = uold(ind_cell_vois,nvar+2)
 
              uold(ind_cell(i),nvar+2)  = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+3) 
           else
              !should be equal to uold(ind_cell(i),10) of the preceeding case 
              uold(ind_cell(i),7) =  uold(ind_cell(i),6 ) + uold(ind_cell_vois,nvar+2) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1)  - uold(ind_cell(i),nvar+3) 

              !ensure div B
              uold(ind_cell(i),nvar+2) =  uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8)  -uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+3) 
           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 


        endif



        if(  xx(i,3) .lt. 2.*dx_min ) then 


             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),6))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 5 (6) is 1 (2)
           if(ind .eq. 5 .or. ind .eq. 6 .or. ind .eq. 7 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - 4*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 5 .or. ind .eq. 6 .or. ind .eq. 7 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+3) = uold(ind_cell_vois,8)
 
              uold(ind_cell(i),8)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),6) - uold(ind_cell(i),7) 
           else
              !should be equal to uold(ind_cell(i),8) of the preceeding case 
              uold(ind_cell(i),nvar+3) =  uold(ind_cell(i), nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell_vois,8) - uold(ind_cell(i),6)  - uold(ind_cell(i),7) 

              !ensure div B
              uold(ind_cell(i),8) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),6) - uold(ind_cell(i),7) 

           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))


           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 


        endif


        if(  xx(i,3) .gt. boxlen-2.*dx_min ) then 


             !look for the grid neigbour of the bottom father
             jgrid = son(nbor(ind_grid(i),5))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
             !we must add 2*ngridmax because the neighbour of 1 (2) is 5 (6)
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 3 .or. ind .eq. 4) then 
             ind_cell_vois = ind_cell_vois + 4*ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag


           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 2 .or. ind .eq. 3 .or. ind .eq. 4) then 
              uold(ind_cell(i),8) = uold(ind_cell_vois,nvar+3)
 
              uold(ind_cell(i),nvar+3)  = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+2) 
           else
              !should be equal to uold(ind_cell(i),nvar+3) of the preceeding case 
              uold(ind_cell(i),8) =  uold(ind_cell(i), 6) + uold(ind_cell(i),7) + uold(ind_cell_vois,nvar+3) - uold(ind_cell(i),nvar+1)  - uold(ind_cell(i),nvar+2) 

              !ensure div B
              uold(ind_cell(i),nvar+3) =  uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8)  -uold(ind_cell(i),nvar+1) - uold(ind_cell(i),nvar+2) 

           endif

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif

        uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

        ! STG HACK CHECK FOR BORKED
        if (uold(ind_cell(i),5) .lt. 0d0) then
           write(*,*) "TOTAL ENERGY < 0 AFTER VELOCITY_FINE, OH NO", ind_cell(i)
           call clean_stop
        end if
        do j=5,8
           if (isnan(uold(ind_cell(i),j))) then
              write(*,*) "VARIABLE IS NAN AFTER VELOCITY_FINE, OH NO", ind_cell(i),j
              call clean_stop
           end if
        end do


       enddo



       
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine boundary_frig
!#########################################################
!#########################################################
!#########################################################
!#########################################################

! CC 03/17
! Get the typical height scale z0 of the gas (mass_all does not work with amr)
subroutine get_height_scale(ilevel,z0,mass_all)
   use amr_commons
   use hydro_commons
   implicit none
   include 'mpif.h'
   integer::ilevel
   integer::i,ind,ncache,igrid,iskip
   integer::info,nleaf,ngrid,nx_loc
   real(dp),dimension(1:3)::skip_loc
   integer::rang,nb_procs,rang_z
   integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf
   real(dp)::dx,vol,scale,z0
   real(kind=8)::mass_loc,mass_max,mass_diff,mass_all
   real(kind=8),allocatable::mass_tab(:)

  nx_loc=icoarse_max-icoarse_min+1
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim

  mass_loc=0.0d0

  call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,info)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rang,info)

  allocate(mass_tab(0:nb_procs-1))

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     
     ! Loop over cells
     do ind=1,twotondim        
        iskip=ncoarse+(ind-1)*ngridmax
        do i=1,ngrid
           ind_cell(i)=ind_grid(i)+iskip
        end do
        
        ! Gather leaf cells
        nleaf=0
        do i=1,ngrid
           if(son(ind_cell(i))==0)then
              nleaf=nleaf+1
              ind_leaf(nleaf)=ind_cell(i)
           end if
        end do

        do i=1,nleaf
           mass_loc=mass_loc+uold(ind_leaf(i),1)*vol
        end do
     
      end do
 
   end do

   !mass_all=0d0
   !z0=0d0
   !mass_tab(:)=0d0
   !mass_max=0d0
   !nb_procs=0
   !rang=0

   call MPI_ALLGATHER(mass_loc,1,MPI_DOUBLE_PRECISION,mass_tab,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(mass_loc,mass_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,info)
   call MPI_ALLREDUCE(mass_loc,mass_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
        
        ! The height z0 is such that M(z0)=Mmax/2**1.5
   mass_tab(:)=abs(mass_tab(:)-mass_max/(2**1.5)) 
 
   mass_diff=mass_tab(0)
   rang_z=0
   do i=1,nb_procs-1
      if(mass_tab(i)<mass_diff) then  
        mass_diff=mass_tab(i)
        rang_z=i
      end if
   end do

   if(rang==rang_z) then
      z0=abs((xg(ind_grid(1),3)-skip_loc(3)-0.5))*scale
   end if
   call MPI_BCAST(z0,1,MPI_DOUBLE_PRECISION,rang_z,MPI_COMM_WORLD,info)

end subroutine get_height_scale
        
           

subroutine get_cell_index4(cell_index,cell_levl,xpart,xtcell,ilevel,np)
  use amr_commons
  implicit none
  integer                                :: np,ilevel
  integer,dimension(1:nvector)           :: cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)      :: xpart
  ! This function returns the index of the cell, at maximum level
  ! ilevel, in which the input particle sits
  real(dp)                               :: xx,yy,zz
  integer                                :: i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0,igrid_old

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim)::xtcell

  integer::idim
  integer::ix,iy,iz
  real(dp)::dx
  real(dp),dimension(1:3)::skip_loc



  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)

  
  !if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  !else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  !else
     !write(*,*)"nx=ny=nz != 1,3 is not supported."
     !call clean_stop
  !end if
  
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,np
     xx = xpart(i,1)
     yy = xpart(i,2)
     zz = xpart(i,3)
     if( ((xx .le. 0) .or. (xx .ge. 1.)) .or. ((yy .le. 0) .or. (yy .ge. 1.)) .or. ((zz .le. 0) .or. (zz .ge. 1.)) ) then 
        cell_index(i)=-1.
     else 
        xx = xx + (nx-1.)/2.
        yy = yy + (ny-1.)/2.
        zz = zz + (nz-1.)/2.
        igrid=igrid0
        do j=1,ilevel
           ii=1; jj=1; kk=1
           if(xx<xg(igrid,1))ii=0
           if(yy<xg(igrid,2))jj=0
           if(zz<xg(igrid,3))kk=0
           ind=1+ii+2*jj+4*kk
           iskip=ncoarse+(ind-1)*ngridmax
           ind_cell=iskip+igrid
           igrid_old = igrid
           igrid=son(ind_cell)
           if(igrid==0.or.j==ilevel)  exit
        end do

        do idim=1,ndim
           xtcell(i,idim)=xg(igrid_old,idim)+xc(ind,idim)
           xtcell(i,idim)=(xtcell(i,idim)-skip_loc(idim))
        end do

        cell_index(i)=ind_cell
        cell_levl(i)=j
     endif
  end do

  
end subroutine get_cell_index4

