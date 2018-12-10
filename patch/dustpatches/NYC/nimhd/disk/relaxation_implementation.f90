subroutine relaxation_implementation(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none

  integer::ilevel
  !--------------------------------------------------------------------------
  ! Enforcing relaxation damping to velocity components at early times,     !
  ! of in geometric regions. This now happens after set_uold, so do to uold !
  !--------------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::ix,iy,iz
  integer::ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell

  real(dp)::dx,scale,v2,rad2,relax_f
  real(dp)::Rin,vcyl,vrad,pres,dens,x_rad,y_rad,z_rad,rad,dx_max,ek,eint,emag
  real(dp),dimension(ndim)::vel

  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel


  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)

  ! Mesh size at level ilevel in coarse cell units
  dx=0.5D0**ilevel
  dx_max = 0.5D0**nlevelmax
  Rin = gravity_params(2)

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Local constants
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)

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
      ! Gather cell centre positions
      do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
         end do
      end do
      ! Rescale position from code units to user units
      do idim=1,ndim
         do i=1,ngrid
            xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
         end do
      end do

      ! To do, allow relax_factor to decay to 1 after some time.
      if ( (nstep .lt. abs(relax_ndecay*relaxation_step)) .or. (t .lt. abs(relax_ndecay*relaxation_time)) ) then
        if ( (relaxation_step .ge. 0.d0) .and. (relaxation_time .ge. 0.d0) )then
          if ( (nstep .lt. abs(relaxation_step)) .or. (t .lt. abs(relaxation_time)) ) then
            relax_f = relax_factor
          else
            if (relaxation_step .ne. 0) then
              relax_f = relax_factor**( (relax_ndecay - nstep/abs(relaxation_step))/(relax_ndecay - 1.d0) )  
            elseif (relaxation_time .ne. 0.d0) then
              relax_f = relax_factor**( (relax_ndecay - t/abs(relaxation_time))/(relax_ndecay - 1.d0) )  
            endif
          endif
          do i=1,ngrid
            uold(ind_cell(i),1+ndim)=uold(ind_cell(i),1+ndim)/relax_f
            x_rad = xx(i,1)-scale*(0.5d0 + dx_max*disk_xc)
            y_rad = xx(i,2)-scale*(0.5d0 + dx_max*disk_yc)
            rad = (x_rad**2 + y_rad**2)**0.5
            ! note, vrad has units dens * vel
            if (rad .gt. 0.d0) then
              vrad = (uold(ind_cell(i),1+1)*x_rad + uold(ind_cell(i),1+2)*y_rad)/rad
              uold(ind_cell(i),1+1)=uold(ind_cell(i),1+1) - (1.d0 - 1.d0/relax_f)*vrad*x_rad/rad
              uold(ind_cell(i),1+2)=uold(ind_cell(i),1+2) - (1.d0 - 1.d0/relax_f)*vrad*y_rad/rad
            else
              vrad = 0.d0
            endif
            ! KE = KE - (dens*Old_Vz)**2/dens/2 + (dens*New_Vz)**2/dens/2
            ! old_Vz = New_Vz*relax_factor
            ! KE = KE - (dens*New_Vz)**2/dens/2*(relax_factor**2 - 1)
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim)- &
              (relax_f**2 - 1.d0)*(vrad**2/relax_f**2 + uold(ind_cell(i),1+3)**2)/ &
              (2.d0*max(uold(ind_cell(i),1),smallr))
          enddo
        else ! forcing vertical and radial vel terms to be zero
          do i=1,ngrid
            uold(ind_cell(i),2+ndim) = uold(ind_cell(i),2+ndim) - uold(ind_cell(i),1+3)**2/ &
                    (2.d0*max(uold(ind_cell(i),1),smallr))
            uold(ind_cell(i),1+ndim)=0.d0
            x_rad = xx(i,1)-scale*(0.5d0 + dx_max*disk_xc)
            y_rad = xx(i,2)-scale*(0.5d0 + dx_max*disk_yc)
            rad = (x_rad**2 + y_rad**2)**0.5
            ! note, vrad has units dens * vel
            if (rad .gt. 0.d0) then
              vrad = (uold(ind_cell(i),1+1)*x_rad + uold(ind_cell(i),1+2)*y_rad)/rad
              uold(ind_cell(i),1+1)=uold(ind_cell(i),1+1) - vrad*x_rad/rad
              uold(ind_cell(i),1+2)=uold(ind_cell(i),1+2) - vrad*y_rad/rad
            else
              vrad = 0.d0
            endif
              uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim)- vrad**2/ &
              (2.d0*max(uold(ind_cell(i),1),smallr))
          enddo
        endif
      endif

      ! Here is where you can enforce certain hydro states; e.g., inner or outer boundaries
      ! Inner/Outer radia > 0 are assumed to be in spherical coordinates,
      !                   < 0 are assumed to be in cylindrical coordinates
      ! Assumption is this is rerunning condinit in these regions.
      if ( relaxation_InnerRadius .gt. 0.d0) then
        do i=1,ngrid
          x_rad = xx(i,1)-scale*(0.5d0 + dx_max*disk_xc)
          y_rad = xx(i,2)-scale*(0.5d0 + dx_max*disk_yc)
          z_rad = xx(i,3)-scale*(0.5d0 + dx_max*disk_zc)
          rad2 = (x_rad**2 + y_rad**2 + z_rad**2)

          if (rad2 .lt. relaxation_InnerRadius**2) then

            ! Get original internal energy
            ek=0.0d0
            do idim=1,ndim
              ek = ek + 0.5d0*uold(ind_cell(i),1)*uold(ind_cell(i),1+idim)**2
            enddo
            eint = uold(ind_cell(i),2+ndim) - ek
            ! Get and save magnetic energy
            emag=0.d0
            do idim=1,ndim
              emag=emag+0.125d0*(uold(ind_cell(i),5+idim)+uold(ind_cell(i),nvar+idim))**2
            enddo

            ! Update variables. relax_revert = 1.0 means you revert to the ICs.
            ! To do, allow relax_revert to decay to 0 after some time.
            uold(ind_cell(i),1)=(1.d0 - relax_revert)*uold(ind_cell(i),1) + relax_revert*dens
            ! velocity -> momentum
            uold(ind_cell(i),2:1+ndim)=(1.d0 - relax_revert)*uold(ind_cell(i),2:1+ndim) + relax_revert*dens*vel(1:ndim)
            ! kinetic energy
            uold(ind_cell(i),2+ndim)=0.d0
            do idim=1,ndim
              uold(ind_cell(i),2+ndim) = uold(ind_cell(i),2+ndim) + 0.5d0*uold(ind_cell(i),1)*uold(ind_cell(i),1+idim)**2
            enddo
            ! pressure -> total fluid energy
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim) + (1.d0 - relax_revert)*eint + relax_revert*pres/(gamma-1.0d0)
            ! magnetic
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim) + emag
          endif
        enddo
      elseif ( relaxation_InnerRadius .lt. 0.d0) then
        do i=1,ngrid
          x_rad = xx(i,1)-scale*(0.5d0 + dx_max*disk_xc)
          y_rad = xx(i,2)-scale*(0.5d0 + dx_max*disk_yc)
          z_rad = xx(i,3)-scale*(0.5d0 + dx_max*disk_zc)
          rad2 = (x_rad**2 + y_rad**2)

          if (rad2 .lt. relaxation_InnerRadius**2) then
            call relax_condinit(x_rad,y_rad,z_rad,dx,dens,vel,pres)

            ! Get original internal energy
            ek=0.0d0
            do idim=1,ndim
              ek = ek + 0.5d0*uold(ind_cell(i),1)*uold(ind_cell(i),1+idim)**2
            enddo
            eint = uold(ind_cell(i),2+ndim) - ek
            ! Get and save magnetic energy
            emag=0.d0
            do idim=1,ndim
              emag=emag+0.125d0*(uold(ind_cell(i),5+idim)+uold(ind_cell(i),nvar+idim))**2
            enddo

            ! Update variables. relax_revert = 1.0 means you revert to the ICs.
            ! To do, allow relax_revert to decay to 0 after some time.
            uold(ind_cell(i),1)=(1.d0 - relax_revert)*uold(ind_cell(i),1) + relax_revert*dens
            
            ! velocity -> momentum
            uold(ind_cell(i),2:1+ndim)=(1.d0 - relax_revert)*uold(ind_cell(i),2:1+ndim) + relax_revert*dens*vel(1:ndim)
            ! kinetic energy
            uold(ind_cell(i),2+ndim)=0.d0
            do idim=1,ndim
              uold(ind_cell(i),2+ndim) = uold(ind_cell(i),2+ndim) + 0.5d0*uold(ind_cell(i),1)*uold(ind_cell(i),1+idim)**2
            enddo
            ! pressure -> total fluid energy
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim) + (1.d0 - relax_revert)*eint + relax_revert*pres/(gamma-1.0d0)
            ! magnetic
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim) + emag
          endif
        enddo
      endif

      if ( relaxation_OuterRadius .gt. 0.d0) then
        do i=1,ngrid
          x_rad = xx(i,1)-scale*(0.5d0 + dx_max*disk_xc)
          y_rad = xx(i,2)-scale*(0.5d0 + dx_max*disk_yc)
          z_rad = xx(i,3)-scale*(0.5d0 + dx_max*disk_zc)
          rad2 = (x_rad**2 + y_rad**2 + z_rad**2)

          if (rad2 .gt. relaxation_OuterRadius**2) then
            call relax_condinit(x_rad,y_rad,z_rad,dx,dens,vel,pres)

            ! Get original internal energy
            ek=0.0d0
            do idim=1,ndim
              ek = ek + 0.5d0*uold(ind_cell(i),1)*uold(ind_cell(i),1+idim)**2
            enddo
            eint = uold(ind_cell(i),2+ndim) - ek
            ! Get and save magnetic energy
            emag=0.d0
            do idim=1,ndim
              emag=emag+0.125d0*(uold(ind_cell(i),5+idim)+uold(ind_cell(i),nvar+idim))**2
            enddo

            ! Update variables. relax_revert = 1.0 means you revert to the ICs.
            ! To do, allow relax_revert to decay to 0 after some time.
            uold(ind_cell(i),1)=(1.d0 - relax_revert)*uold(ind_cell(i),1) + relax_revert*dens
            ! velocity -> momentum
            uold(ind_cell(i),2:1+ndim)=(1.d0 - relax_revert)*uold(ind_cell(i),2:1+ndim) + relax_revert*dens*vel(1:ndim)
            ! kinetic energy
            uold(ind_cell(i),2+ndim)=0.d0
            do idim=1,ndim
              uold(ind_cell(i),2+ndim) = uold(ind_cell(i),2+ndim) + 0.5d0*uold(ind_cell(i),1)*uold(ind_cell(i),1+idim)**2
            enddo
            ! pressure -> total fluid energy
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim) + (1.d0 - relax_revert)*eint + relax_revert*pres/(gamma-1.0d0)
            ! magnetic
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim) + emag
          endif
        enddo
      elseif ( relaxation_OuterRadius .lt. 0.d0) then
        do i=1,ngrid
          x_rad = xx(i,1)-scale*(0.5d0 + dx_max*disk_xc)
          y_rad = xx(i,2)-scale*(0.5d0 + dx_max*disk_yc)
          z_rad = xx(i,3)-scale*(0.5d0 + dx_max*disk_zc)
          rad2 = (x_rad**2 + y_rad**2)

          if (rad2 .gt. relaxation_OuterRadius**2) then
            call relax_condinit(x_rad,y_rad,z_rad,dx,dens,vel,pres)

            ! Get original internal energy
            ek=0.0d0
            do idim=1,ndim
              ek = ek + 0.5d0*uold(ind_cell(i),1)*uold(ind_cell(i),1+idim)**2
            enddo
            eint = uold(ind_cell(i),2+ndim) - ek
            ! Get and save magnetic energy
            emag=0.d0
            do idim=1,ndim
              emag=emag+0.125d0*(uold(ind_cell(i),5+idim)+uold(ind_cell(i),nvar+idim))**2
            enddo

            ! Update variables. relax_revert = 1.0 means you revert to the ICs.
            ! To do, allow relax_revert to decay to 0 after some time.
            uold(ind_cell(i),1)=(1.d0 - relax_revert)*uold(ind_cell(i),1) + relax_revert*dens
            ! velocity -> momentum
            uold(ind_cell(i),2:1+ndim)=(1.d0 - relax_revert)*uold(ind_cell(i),2:1+ndim) + relax_revert*dens*vel(1:ndim)
            ! kinetic energy
            uold(ind_cell(i),2+ndim)=0.d0
            do idim=1,ndim
              uold(ind_cell(i),2+ndim) = uold(ind_cell(i),2+ndim) + 0.5d0*uold(ind_cell(i),1)*uold(ind_cell(i),1+idim)**2
            enddo
            ! pressure -> total fluid energy
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim) + (1.d0 - relax_revert)*eint + relax_revert*pres/(gamma-1.0d0)
            ! magnetic
            uold(ind_cell(i),2+ndim)=uold(ind_cell(i),2+ndim) + emag
          endif
        enddo
      endif

    enddo
  enddo

111 format('   Entering relaxation_implementation for level ',I2)

end subroutine relaxation_implementation
