!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell,v)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx,dx_max                     ! Cell size
  real(dp),dimension(1:nvector,1:3)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  real(dp),optional,intent(in),dimension(1:nvector,1:3)::v ! Cell center velocity
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) are cell center velocity in grid units (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  integer::idim,i,i_table,ir,ir1,ir2
  real(dp)::gmass,emass,cmass,FofC,FofR,xmass,ymass,zmass,rr,rc,rx,ry,rz,Rmin_grav,Rmax_grav,kpc=3.085678d21
  real(dp)::vcyl_rc,f_table,maxRmass,rmass,fr,dr_b,dummy,acc_bulge_x,acc_bulge_y,acc_bulge_z
  logical, save :: first=.true.
  integer, parameter :: nbulge=1000
  real(dp), save, dimension(0:nbulge) :: save_bulge

  real(dp), save ::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,pi,gm_rescale,zmin,zmax,dz

  if (first) then
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    ! units only used for gravity_type=3, where first needs to remain true until end of subroutine
  endif

  pi = acos(-1.d0)

  dx_max = boxlen/2.d0**nlevelmax

  ! Constant vector
  if(gravity_type==1)then
     do idim=1,ndim
        do i=1,ncell
           f(i,idim)=gravity_params(idim)
        end do
     end do
  end if

  ! Point mass
  if(gravity_type==2 .or. gravity_type==-2)then
     if (gravity_type==2) f = 0.d0
     gmass=gravity_params(1) ! GM in code units
     emass=dx
     emass=gravity_params(2) ! Softening length in code units
     xmass=boxlen/2.d0 + gravity_params(3)*dx_max ! Point mass coordinates relative to centre of box
     ymass=boxlen/2.d0 + gravity_params(4)*dx_max
     zmass=boxlen/2.d0 + gravity_params(5)*dx_max
     Rmin_grav=gravity_params(6) ! In code units
     Rmax_grav=gravity_params(7) ! In code units

     if (Rmax_grav .eq. 0.d0) Rmax_grav = 100.d0*boxlen
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        if (disk_testR .lt. 0.d0) then
          rx=x(i,1)-xmass
#if NDIM>1
          ry=x(i,2)-ymass
#endif
        else
          rx = disk_testR
#if NDIM>1
          ry=0.d0
#endif
        endif
#if NDIM>2
        if (disk_testz .lt. 0.d0) then
          rz=x(i,3)-zmass
        else
          rz = disk_testz
        endif
#endif
        rr=sqrt(rx**2+ry**2+rz**2)
        rc=sqrt(rx**2+ry**2)

        if ( (Rmin_grav < rc) .and. (rc < Rmax_grav) ) then

          if (disk_setup .lt. 10001) then
            if (rc .gt. emass) then
              f(i,1)=f(i,1) - gmass*rx/rr**3
#if NDIM>1
              f(i,2)=f(i,2) - gmass*ry/rr**3
#endif
#if NDIM>2
              f(i,3)=f(i,3) - gmass*rz/rr**3
#endif
            else
              if (disk_setup .eq. 201) then
                f_table = abs(rz)/(boxlen/2.d0)*2**(disk_TableLevel-1) + 0.5d0
                i_table = max(min(int(f_table),2**(disk_TableLevel-1)-1),1)
                f_table = f_table - i_table
                if (rc < disk_vramp*emass) then
                  vcyl_rc = (1.d0/(disk_vramp*emass))*(vrot_table(i_table)*(1.d0-f_table) + &
                                                       vrot_table(i_table+1)*f_table)
                else
                  vcyl_rc = (1.d0/rc)*(vrot_table(i_table)*(1.d0-f_table) + &
                                       vrot_table(i_table+1)*f_table)
                endif
                f(i,1)=f(i,1) - vcyl_rc**2*rx
#if NDIM>1
                f(i,2)=f(i,2) - vcyl_rc**2*ry
#endif
#if NDIM>2
                f(i,3)=f(i,3) - gmass*rz/(sqrt(rz**2 + emass**2))**3
#endif
              else
                stop "Error: Not supported for gravity with new dPdR"

!                if (disk_setup .eq. 1111) then
!                  pres = C3*xscl**3 + C2*xscl**2 + C1*xscl + C0
!                  if (r_cyl .gt. Rin_2.5d0*dR) then
!                    dPdR = 3.d0*C3*xscl**2 + 2.0d0*C2*xscl + C1
!                  else
!                    dPdR = 0.d0
!                  endif
!                  vcyl_ = ?
!                elseif (disk_setup .eq. 11111) then
!                  temp =
!                  rho =
!                  pres =
!                  if (r_cyl .gt. Rin_2.5d0*dR) then
!                    dPdR = rho*(3.d0*C3*xscl**2 + 2.0d0*C2*xscl + C1) + &
!                           temp*rho*( )
!                  else
!                    dPdR = 0.d0
!                  endif
!                endif

              endif
            endif
          elseif (disk_setup .lt. 24001 .or. disk_setup .ge. 30001) then
            if (grav_angle .gt. 0.d0) then
              zmin = grav_angle*(1.d0 - 0.5d0*grav_width)*max(rc,sqrt(emass**2 - rz**2))
              zmax = grav_angle*(1.d0 + 0.5d0*grav_width)*max(rc,sqrt(emass**2 - rz**2))
              dz = zmax - zmin
              gm_rescale = gmass * (cos( min(1.d0,max(0.d0,(abs(rz)-zmin)/dz))*pi ) + 1.d0)/2.d0
            else
              gm_rescale = gmass
            endif
            if (disk_testR .lt. 0.d0) then
              f(i,1) = f(i,1) - gm_rescale*rx/max(rr,emass)**3
#if NDIM>1
              f(i,2) = f(i,2) - gm_rescale*ry/max(rr,emass)**3
#endif
            endif
#if NDIM>2
            if (disk_testz .lt. 0.d0) then
              f(i,3) = f(i,3) - gm_rescale*rz/max(rr,emass)**3
            endif
#endif
          else
            if (grav_angle .gt. 0.d0) then
              zmin = grav_angle*(1.d0 - 0.5d0*grav_width)*max(rc,sqrt(emass**2 - rz**2))
              zmax = grav_angle*(1.d0 + 0.5d0*grav_width)*max(rc,sqrt(emass**2 - rz**2))
              dz = zmax - zmin
              gm_rescale = gmass * (cos( min(1.d0,max(0.d0,(abs(rz)-zmin)/dz))*pi ) + 1.d0)/2.d0
            else
              gm_rescale = gmass
            endif
            if (rr .lt. emass) then
              gm_rescale = gm_rescale * (rr/emass) * (99.d0/8.d0*(rr/emass)**3-77.d0/4.d0*(rr/emass)**4+63.d0/8.d0*(rr/emass)**5)**2
            endif
            if (disk_testR .lt. 0.d0) then
              if (rr .gt. 0.d0) then
                f(i,1) = f(i,1) - gm_rescale*rx/rr**3
#if NDIM>1
                f(i,2) = f(i,2) - gm_rescale*ry/rr**3
#endif
              endif

            endif
#if NDIM>2
            if (disk_testz .lt. 0.d0) then
              if (rr .gt. 0.d0) then
                f(i,3) = f(i,3) - gm_rescale*rz/rr**3
              endif
            endif
#endif
          endif
        else
          f(i,1:ndim) = 0.d0
        endif
     end do

  else if (abs(gravity_type)==3)then
#if NDIM<3
     stop "NFW + Bulge only in 3D"
#endif
     gmass=gravity_params(1)*6.67d-8/scale_l**3*scale_t**2 ! MDM * G in code units, params should have M200_DM in Msun
     emass=gravity_params(2) ! Softening length
     xmass=gravity_params(3) ! Center mass coordinates
     ymass=gravity_params(4)
     zmass=gravity_params(5)
     cmass=gravity_params(6) ! concentration
     rmass=gravity_params(7) ! virial radius
     maxRmass=gravity_params(8) ! maximum radius beyond which grav isn't included

     ! File units are in kpc , (km/s)**2/kpc
     if (first) then ! read data
       first = .false.
       ! Read Bulge info: first read r into dr_b, second into fr. Bulge is spherical
       open(42,file='Bulge_acc_1d.dat',status='old')
       read(42,*) dr_b, save_bulge(0)
       read(42,*) fr, save_bulge(1)
       do ir=2,nbulge
         read(42,*) dummy, save_bulge(ir)
       enddo
       ! Units are (km/s)**2 / kpc
       save_bulge = save_bulge * 1d10 / kpc
       ! Units are cm/s**2
       save_bulge = save_bulge/scale_l*scale_t**2
       ! Units are in code_length now

       ! set dr / dz
       dr_b = fr - dr_b
       close(42)
     endif


     FofC = log(1.d0 + cmass) - cmass/(1.d0 + cmass)
     do i=1,ncell
        rx=0.0d0; ry=0.0d0; rz=0.0d0
        rx=x(i,1)-boxlen/2.d0-xmass
        ry=x(i,2)-boxlen/2.d0-ymass
        rz=x(i,3)-boxlen/2.d0-zmass

        rr=max(sqrt(rx**2+ry**2+rz**2),emass)

        if (rr < maxRmass) then
          FofR = log(1.d0 + rr/(rmass/cmass)) - rr/(rmass/cmass)/(1.d0 + rr/(rmass/cmass))

          fr = rr/dr_b
          ir1 = int(fr)
          fr = fr - ir1
          ir2 = min(ir1+1,nbulge)
          acc_bulge_x = (fr*save_bulge(ir2) + (1.d0-fr)*save_bulge(ir1))*rx/rr
          acc_bulge_y = (fr*save_bulge(ir2) + (1.d0-fr)*save_bulge(ir1))*ry/rr
          acc_bulge_z = (fr*save_bulge(ir2) + (1.d0-fr)*save_bulge(ir1))*rz/rr

          f(i,1)=f(i,1)-Gmass*rx/rr**3*FofR/FofC + acc_bulge_x
          f(i,2)=f(i,2)-Gmass*ry/rr**3*FofR/FofC + acc_bulge_y
          f(i,3)=f(i,3)-Gmass*rz/rr**3*FofR/FofC + acc_bulge_z
        else
          f(i,1:3) = 0.d0
        endif
     end do

  end if


end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=4.D0*ACOS(-1.0D0)

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2d0*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2d0*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana
