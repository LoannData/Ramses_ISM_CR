subroutine courant_fine(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use radiation_parameters,only:frad,dtdiff_params
#if USE_TURB==1
  use turb_commons
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
  real(kind=8),dimension(4)::comm_buffin,comm_buffout
#endif
  integer::ilevel
  !----------------------------------------------------------------------
  ! Using the Courant-Friedrich-Levy stability condition,               !
  ! this routine computes the maximum allowed time-step.                !
  !----------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::nleaf,ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf

  real(dp)::dt_lev,dx,vol,scale
  real(kind=8)::mass_loc,ekin_loc,eint_loc,emag_loc,dt_loc
  real(kind=8)::mass_all,ekin_all,eint_all,emag_all,dt_all
  real(dp),dimension(1:nvector,1:nvar+3+ndust),save::uu
  real(dp),dimension(1:nvector,1:ndim),save::gg
#if NDUST>0
  real(dp)::dt_dust,vdust
  integer::idust
  !real(dp),dimension(1:nvector)::vdu
#endif
#if NIMHD==1
  ! modif nimhd
  real(dp)::dtwad_loc,dtwad_lev,dtwad_all
  real(dp)::dtambdiff_loc,dtambdiff_lev,dtambdiff_all
  real(dp)::dtmagdiff_loc,dtmagdiff_lev,dtmagdiff_all
  real(dp)::dthall_loc,dthall_lev,dthall_all
  real(dp)::tmag1,tmag2
  ! fin modif nimhd
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  mass_all=0.0d0; mass_loc=0.0d0
  ekin_all=0.0d0; ekin_loc=0.0d0
  emag_all=0.0d0; emag_loc=0.0d0
  eint_all=0.0d0; eint_loc=0.0d0
  
  dt_all=dtnew(ilevel); dt_loc=dt_all
#if NIMHD==1
  ! modif nimhd
  dtambdiff_all=dtambdiff(ilevel); dtambdiff_loc=dtambdiff_all
  dtmagdiff_all=dtmagdiff(ilevel); dtmagdiff_loc=dtmagdiff_all
  dtwad_all=dtwad(ilevel); dtwad_loc=dtwad_all
  dthall_all=dthall(ilevel); dthall_loc=dthall_all
  ! fin modif nimhd
#endif

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  vol=dx**ndim
  call velocity_fine(ilevel)

  if (ischeme .eq. 1) then
     CALL velocity_fine(ilevel)
     do ivar=1,nvar+3
        call make_virtual_fine_dp(uold(1,ivar),ilevel)
     end do
     if(simple_boundary)call make_boundary_hydro(ilevel)
  endif

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

        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nleaf
              uu(i,ivar)=uold(ind_leaf(i),ivar)
           end do
        end do
#if NDUST>0
        do i=1,nleaf
        !do idim=1,ndim
           do idust=1,ndust
#if NDIM==1              
              uu(i,nvar+3+idust)=  abs(v_dust(ind_leaf(i),idust,1))
#endif
#if NDIM==2              
              uu(i,nvar+3+idust)= abs(v_dust(ind_leaf(i),idust,1))+abs(v_dust(ind_leaf(i),idust,2))
#endif
#if NDIM==3              
              uu(i,nvar+3+idust)= abs(v_dust(ind_leaf(i),idust,1))+abs(v_dust(ind_leaf(i),idust,2))+abs(v_dust(ind_leaf(i),idust,3))
#endif
              if(.not.dust_diffusion) uu(i,nvar+3+idust)=0.0d0

           end do
        !end do
        end do
#endif
        ! Gather gravitational acceleration
        gg=0.0d0
        if(poisson)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=f(ind_leaf(i),idim)
              end do
           end do
        end if
        ! Gather radiative force
#if USE_FLD==1 || USE_M_1==1
        if(fld)then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=gg(i,idim)+frad(ind_leaf(i),idim)
              end do
           end do
        end if
#endif
        ! Gather turbulent force
#if USE_TURB==1
        if (turb .AND. turb_type/=3) then
           do idim=1,ndim
              do i=1,nleaf
                 gg(i,idim)=gg(i,idim)+fturb(ind_leaf(i),idim)
              end do
           end do
        end if
#endif

        ! Compute total mass
        do i=1,nleaf
           mass_loc=mass_loc+uu(i,1)*vol
        end do

        ! Compute total energy
        do i=1,nleaf
           ekin_loc=ekin_loc+uu(i,5)*vol
        end do

        ! Compute total magnetic energy
        do ivar=1,3
           do i=1,nleaf
              emag_loc=emag_loc+0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do

        ! Compute total internal energy
        do i=1,nleaf
           eint_loc=eint_loc+uu(i,5)*vol
        end do
        do ivar=1,3
           do i=1,nleaf
              eint_loc=eint_loc-0.5d0*uu(i,1+ivar)**2/uu(i,1)*vol &
                   & -0.125d0*(uu(i,5+ivar)+uu(i,nvar+ivar))**2*vol
           end do
        end do
#if NENER>0
        do ivar=1,nener
           do i=1,nleaf
              eint_loc=eint_loc-uu(i,8+ivar)*vol
           end do
        end do
#endif

        ! Compute CFL time-step
        if(nleaf>0)then

#if NIMHD==1
           ! modif nimhd

           call cmpdt(uu,gg,dx,dt_lev,nleaf,dtambdiff_lev,dtmagdiff_lev,dthall_lev)
           dt_loc=min(dt_loc,dt_lev)
           dtambdiff_loc=min(dtambdiff_loc,dtambdiff_lev)
           dtmagdiff_loc=min(dtmagdiff_loc,dtmagdiff_lev)
           dtwad_loc=min(dtwad_loc,dt_lev)
           dthall_loc=min(dthall_loc,dthall_lev)
           ! fin modif nimhd
#else
           call cmpdt(uu,gg,dx,dt_lev,nleaf)
           dt_loc=min(dt_loc,dt_lev)
#endif
        end if

     end do
     ! End loop over cells

  end do
  ! End loop over grids

  ! Compute global quantities
#ifndef WITHOUTMPI
  comm_buffin(1)=mass_loc
  comm_buffin(2)=ekin_loc
  comm_buffin(3)=eint_loc
  comm_buffin(4)=emag_loc
  call MPI_ALLREDUCE(comm_buffin,comm_buffout,4,MPI_DOUBLE_PRECISION,MPI_SUM,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dt_loc     ,dt_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
#if NIMHD==1
  ! modif nimhd
  call MPI_ALLREDUCE(dtambdiff_loc     ,dtambdiff_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dtmagdiff_loc     ,dtmagdiff_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dtwad_loc     ,dtwad_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dthall_loc     ,dthall_all      ,1,MPI_DOUBLE_PRECISION,MPI_MIN,&
       &MPI_COMM_WORLD,info)
  ! fin modif nimhd
#endif
  mass_all=comm_buffout(1)
  ekin_all=comm_buffout(2)
  eint_all=comm_buffout(3)
  emag_all=comm_buffout(4)
#endif
#ifdef WITHOUTMPI
  mass_all=mass_loc
  ekin_all=ekin_loc
  eint_all=eint_loc
  emag_all=emag_loc
  dt_all=dt_loc
#if NIMHD==1
  ! modif nimhd
  dtambdiff_all=dtambdiff_loc
  dtmagdiff_all=dtmagdiff_loc
  dtwad_all=dtwad_loc
  dthall_all=dthall_loc
  ! fin modif nimhd
#endif
#endif

  mass_tot=mass_tot+mass_all
  ekin_tot=ekin_tot+ekin_all
  eint_tot=eint_tot+eint_all
  emag_tot=emag_tot+emag_all
  dtnew(ilevel)=MIN(dtnew(ilevel),dt_all)

#if NIMHD==1
  ! modif nimhd
  dtambdiff(ilevel)=MIN(dtambdiff(ilevel), dtambdiff_all)
  dtmagdiff(ilevel)=MIN(dtmagdiff(ilevel), dtmagdiff_all)
  dtwad(ilevel)=MIN(dtwad(ilevel), dtwad_all) 
  dthall(ilevel)=MIN(dthall(ilevel), dthall_all)
  
!   if(myid==1) write(*,*) 'dtnew(ilevel),dtwad(ilevel)',dtnew(ilevel),dtwad(ilevel)

  tmag1 = 1.0e+30 ; tmag2 = 1.0e+30
  
  if(nambipolar.eq.1) then
     ! ambipolar diffusion
     ! WARNING this should not be done for tests
     if (nminitimestep.eq.1) then
        ! alfven time alone maybe not correct
        ! comparison with global time step
        tmag1=max(dtambdiff(ilevel),dtwad(ilevel)*coefalfven)
!         if(myid==1) write(*,*) 'tmag,dtambdiff(ilevel),dtnew(ilevel)*coefalfven,dtnew(ilevel)',tmag1,dtambdiff(ilevel),dtnew(ilevel)*coefalfven,dtnew(ilevel)
     else
        tmag1=dtambdiff(ilevel)
     endif
!neil      dtnew(ilevel)=MIN(dtnew(ilevel),tmag)
     ! Barenblatt ambipolar diffusion
     !dtnew(ilevel)=dtambdiff(ilevel)
  endif

  if  (nmagdiffu == 1) then
     if (nminitimestep.eq.1) then
        ! alfven time alone maybe not correct
        ! comparison with global time step
        tmag2=max(dtmagdiff(ilevel),dtwad(ilevel)*coefdtohm)
!         if(myid==1) write(*,*) 'tmag,dtmagdiff(ilevel),dtnew(ilevel)*coefdtohm,dtnew(ilevel)',tmag2,dtmagdiff(ilevel),dtnew(ilevel)*coefdtohm,dtnew(ilevel)
     else
        tmag2=dtmagdiff(ilevel)
     endif
!neil      dtnew(ilevel)=MIN(dtnew(ilevel),tmag,dthall(ilevel))
  end if
  
!   if  (nmagdiffu2 == 0) then
! !neil     dtnew(ilevel)=MIN(dtnew(ilevel),dtmagdiff(ilevel),dthall(ilevel))
!      dtnew(ilevel)=MIN(dtnew(ilevel),tmag1,tmag2,dthall(ilevel))
!   else
!      dtnew(ilevel)=MIN(dtnew(ilevel),tmag1,dthall(ilevel))     ! subcycling, on ne prend pas en compte le temps ohmique
!   end if
  
  if  (nambipolar2 == 0) then ! no subcycling, on prend en compte le temps ambipolaire
     dtnew(ilevel)=MIN(dtnew(ilevel),tmag1)
  endif
  if  (nmagdiffu2 == 0) then ! no subcycling, on prend en compte le temps ohmique
     dtnew(ilevel)=MIN(dtnew(ilevel),tmag2)
  endif
  
  ! and finally Hall dt
  dtnew(ilevel)=MIN(dtnew(ilevel),dthall(ilevel))

  ! Magnetic diffusion alone
  !dtnew(ilevel)=dtmagdiff(ilevel)
  
  ! fin modif nimhd
#endif

  if(dt_control)dtnew(ilevel)=dtdiff_params(1)*dtdiff_params(2)**nstep_coarse

111 format('   Entering courant_fine for level ',I2)

end subroutine courant_fine

!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine velocity_fine(ilevel)
  use amr_commons      !, ONLY: dp,ndim,nvector,boxlen,t
!  use hydro_parameters !, ONLY: nvar,boundary_var,gamma,bx_bound,by_bound,bz_bound,turb,dens0,V0
  use hydro_commons
  use poisson_parameters
  implicit none
  integer::ilevel
  !----------------------------------------------------------
  ! This routine computes the gravitational acceleration,
  ! the maximum density rho_max, and the potential energy
  !----------------------------------------------------------
  integer::igrid,ngrid,ncache,i,ind,iskip,ix,iy,iz
  integer::info,ibound,nx_loc,idim,neul=5
  real(dp)::dx,dx_loc,scale,d,u,v,w,A,B,C
  real(kind=8)::rho_max_loc,rho_max_all,epot_loc,epot_all
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:3)::skip_loc

  integer ,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp),dimension(1:nvector,1:ndim),save::x
  real(dp),dimension(1:nvector,1:3),save::vv
  real(dp),dimension(1:nvector,1:nvar+3)::q   ! Primitive variables
  real(dp)::pi,time
  integer ::ivar,jgrid,ind_cell_vois
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,Cwnm
  real(dp)::dx_min, fact, Emag,Emag0,Cs,Cs0,rc,r0,r_acc
  real(dp)::omega0,x0,y0,z0,M00,emass,xx,yy,zz,dens1,d00

!  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

!  Cwnm = sqrt(8000./scale_T2)

  Cs0=sqrt(temper_iso)

  pi=ACOS(-1.0d0)

!  time = t * Cwnm / boxlen

 

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
              x(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
           end do
        end do
        ! Rescale position from code units to user units
        do idim=1,ndim
           do i=1,ngrid
              x(i,idim)=(x(i,idim)-skip_loc(idim))/scale
           end do
        end do
        
  M00=gravity_params(1)
  emass=gravity_params(2)

  x0=gravity_params(3)
  y0=gravity_params(4)
  z0=gravity_params(5)
  omega0 = sqrt(M00)
  r0 = 5.5!boxlen / 8.
  d00=1e-11/scale_d

       do i=1,ngrid


     xx=x(i,1)-x0
     yy=x(i,2)-y0
     zz=x(i,3)-z0
     rc=sqrt(xx**2+yy**2)
     dens1=1.d-22/scale_d

     Cs=max(Cs0*( (rc**5+(10.*emass)**5)**(1./5.)/r0)**(-temper_expo/2.),Cs0/2.)

     Emag=0.

        !confine the gas outside the disk close to  the equatorial plane 
        !this is to force it to be accreted from the equatorial plane
        if(rc .ge. 0.5*r0 .and. confine) then 
           if(zz .ge. 0.) uold(ind_cell(i),4) = min(uold(ind_cell(i),4) , 0.)
           if(zz .le. 0.) uold(ind_cell(i),4) = max(uold(ind_cell(i),4) , 0.)
        endif

        if(rc .ge. 1.*r0 .and. confine .and. uold(ind_cell(i),1) .gt. 1.d-5) then 
!           if(zz .ge. 0.) uold(ind_cell(i),4) = min(uold(ind_cell(i),4) , -5.*Cs*min(abs(zz)/r0,1.)*uold(ind_cell(i),1) )
!           if(zz .le. 0.) uold(ind_cell(i),4) = max(uold(ind_cell(i),4) ,  5.*Cs*min(abs(zz)/r0,1.)*uold(ind_cell(i),1) )
           if(zz .ge. 0.) uold(ind_cell(i),4) =  -5.*Cs*min(abs(zz)/r0,1.)*uold(ind_cell(i),1) 
           if(zz .le. 0.) uold(ind_cell(i),4) =   5.*Cs*min(abs(zz)/r0,1.)*uold(ind_cell(i),1) 
        endif

        if(rc .ge. 3.5*r0) then 
      
           dens1 =  max(d00 * r0**2 / sqrt(rc**2 + emass**2)**(3-temper_expo/2.) * exp( - M00 / Cs**2 / 2.  * zz**2 / (rc**2+emass**2)**(1.5) ) /1000. ,1.d-22/scale_d)
           uold(ind_cell(i),1) = dens1

          
           uold(ind_cell(i),2) = dens1*   omega0 * yy / (rc**2+emass**2)**(0.75)
           uold(ind_cell(i),3) = dens1*( -omega0 * xx / (rc**2+emass**2)**(0.75))
           uold(ind_cell(i),4) = 0.

           uold(ind_cell(i),6)  = 0.!sqrt(dens0*2.)*Cs*bx_bound
           uold(ind_cell(i),9)  = 0.!sqrt(dens0*2.)*Cs*bx_bound

           uold(ind_cell(i),7)  = 0.!sqrt(dens0*2.)*Cs*by_bound
           uold(ind_cell(i),10) = 0.!sqrt(dens0*2.)*Cs*by_bound

           uold(ind_cell(i),8)  = 0.!sqrt(dens0*2.)*Cs*bz_bound
           uold(ind_cell(i),11) = 0.!sqrt(dens0*2.)*Cs*bz_bound



!          if(accret_circ) then 

!           uold(ind_cell(i),1) = dens0 + uold(ind_cell(i),1)


             !add some mass flux
!             uold(ind_cell(i),2) = dens0*U0*Cs + uold(ind_cell(i),2)


             !gives enough momentum for the gas to reach the outer part of the disk 
!             uold(ind_cell(i),3) = dens0*V0 * sqrt(gravity_params(1)/r0) *(r0/rc) + uold(ind_cell(i),3)


!          else 

          r_acc = sqrt(yy**2+(zz-shifting*r0)**2)

          if( r_acc .lt. r0/accret_rad .and. xx .lt. 0. .and. accret_nonsym .eqv. .true.) then 

           uold(ind_cell(i),1) = dens0 + uold(ind_cell(i),1)


             !add some mass flux
             uold(ind_cell(i),2) = dens0*U0*Cs + uold(ind_cell(i),2)


             !gives enough momentum for the gas to reach the outer part of the disk 
             uold(ind_cell(i),3) = dens0*V0 * sqrt(gravity_params(1)/r0) *(r0/rc) + uold(ind_cell(i),3)

!           endif

           else if (accret_nonsym .eqv. .false. .and. abs(zz) .lt. r0 ) then

             !the total mass flux is the same for the 2 accretion schemes
             !pi r0^2 / (2.*pi*3.5*r0*r0) = > 1/(2*3.5)
             uold(ind_cell(i),1) = dens0 / (2.*3.5) + uold(ind_cell(i),1)

             !add some mass flux
             uold(ind_cell(i),2) = -dens0 / (2.*3.5)*U0*Cs*xx/rc + uold(ind_cell(i),2)

             uold(ind_cell(i),3) = -dens0 / (2.*3.5)*U0*Cs*yy/rc + uold(ind_cell(i),3)

             !gives enough momentum for the gas to reach the outer part of the disk 
             uold(ind_cell(i),2) = dens0 / (2.*3.5)*V0 * (yy/rc) * sqrt(gravity_params(1)/r0) *(r0/rc) + uold(ind_cell(i),2)

             uold(ind_cell(i),3) = dens0 / (2.*3.5)*V0 * (-xx/rc) * sqrt(gravity_params(1)/r0) *(r0/rc) + uold(ind_cell(i),3)


           endif


           uold(ind_cell(i),neul)= uold(ind_cell(i),1)*Cs**2/(gamma-1.) + 0.5/uold(ind_cell(i),1)*(uold(ind_cell(i),2)**2+uold(ind_cell(i),3)**2+uold(ind_cell(i),4)**2) + Emag


          endif




       enddo



       
     end do
     ! End loop over cells

  end do
  ! End loop over grids

end subroutine velocity_fine
!#########################################################
!#########################################################
!#########################################################
!#########################################################
