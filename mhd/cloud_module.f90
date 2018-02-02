module cloud_module
  use amr_parameters
  use hydro_parameters,only:Msun
  
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

  real(dp):: C2_vis=0.0d0 !Von Neumann & Richtmeyer artificial viscosity coefficient 3 en principe
  real(dp):: alpha_dense_core=0.5d0
  real(dp):: beta_dense_core=0.0d0
  real(dp):: crit=0.0d0
  real(dp):: delta_rho=0.0d0
  real(dp):: Mach=0.0d0
  real(dp):: r0_box=4.0d0

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
  real(dp):: facc_star_mom=1.0d0       ! fraction of the angular momentum accreted by the sinks
  integer::nmdot_PMS,nm_PMS,ndata_PMS
  integer ,allocatable,dimension(:)::nb_ligne_PMS
  real(dp),allocatable,dimension(:,:,:)::data_PMS

end module cloud_module

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
  namelist/cloud_params/bl_fac
!  namelist/cloud_params/alpha_dense_core,beta_dense_core,crit,delta_rho &
!       & ,mass_c,rap,cont,ff_sct,ff_rt,ff_act,ff_vct,theta_mag,bb_test &
!       & ,contrast,Mach,uniform_bmag

  ! Read namelist file
  rewind(1)
  read(1,NML=cloud_params,END=101)
101 continue                                   ! No harm if no namelist

  ! Get some units out there
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)



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



  !if you want to activate zero gradient boundary conditions with proper div B 
  ! remove the return
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


        !impose vanishing gradient conditions at the x  faces
        if(  xx(i,1) .lt. 2.*dx_min ) then 

             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),2))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
           if(ind .eq. 2 .or. ind .eq. 4 .or. ind .eq. 6 .or. ind .eq. 8) then 
             ind_cell_vois = ind_cell_vois - ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! Subtract non-thermal pressure terms
!#if NENER>0
!           do irad=1,nener
!              uold(ind_cell(i),5) = uold(ind_cell(i),5) - uold(ind_cell(i),8+irad)
!           end do
!#endif

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 2 .or. ind .eq. 4 .or. ind .eq. 6 .or. ind .eq. 8) then 
              uold(ind_cell(i),nvar+1) = uold(ind_cell_vois,6)
 

              uold(ind_cell(i),6)  = uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),7) - uold(ind_cell(i),8) 
           else
              !should be equal to uold(ind_cell(i),7) of the preceeding case 
              uold(ind_cell(i),nvar+1) =  uold(ind_cell_vois,6) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3) - uold(ind_cell(i),7)  - uold(ind_cell(i),8) 

              !ensure div B
              uold(ind_cell(i),6) =  uold(ind_cell(i),nvar+1) + uold(ind_cell(i),nvar+2) + uold(ind_cell(i),nvar+3)  -uold(ind_cell(i),7) - uold(ind_cell(i),8) 
           endif




           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

           ! Add back the non-thermal pressure
!#if NENER>0
!           do irad=1,nener
!              uold(ind_cell(i),5) = uold(ind_cell(i),5) + uold(ind_cell(i),8+irad)
!           end do
!#endif

        endif



        !impose vanishing gradient conditions at the x  faces
        if(  xx(i,1) .gt. boxlen-2.*dx_min ) then 

             !look for the grid neigbour of the top father
             jgrid = son(nbor(ind_grid(i),1))

           ind_cell_vois = iskip + jgrid 
             !remember iskip is calculated above
           if(ind .eq. 1 .or. ind .eq. 3 .or. ind .eq. 5 .or. ind .eq. 7) then 
             ind_cell_vois = ind_cell_vois + ngridmax
           endif

           uold(ind_cell(i),1:nvar+3) =  uold(ind_cell_vois,1:nvar+3) 
           uold(ind_cell(i),1) = MAX(uold(ind_cell(i),1),smallr)

           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) = uold(ind_cell(i),5)  - Emag

           ! we have to modify the 2 normal components of the magnetic field
           if(ind .eq. 1 .or. ind .eq. 3 .or. ind .eq. 5 .or. ind .eq. 7) then 
              uold(ind_cell(i),6) = uold(ind_cell_vois,nvar+1)
 
              uold(ind_cell(i),nvar+1) = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 
           else
              !should be equal to uold(ind_cell(i),9) of the preceeding case 
              uold(ind_cell(i),6) =  uold(ind_cell(i),7) + uold(ind_cell(i),8)  + uold(ind_cell_vois,nvar+1) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 

              !ensure div B
              uold(ind_cell(i),nvar+1) = uold(ind_cell(i),6) + uold(ind_cell(i),7) + uold(ind_cell(i),8) - uold(ind_cell(i),nvar+2) - uold(ind_cell(i),nvar+3) 
           endif


           A=0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
           B=0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
           C=0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))

           Emag = 0.5*(A**2+B**2+C**2)

           uold(ind_cell(i),5) =  uold(ind_cell(i),5) + Emag 

        endif





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
