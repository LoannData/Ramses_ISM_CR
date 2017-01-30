subroutine backup_sink(filename)
  use amr_commons
  use pm_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'  
#endif 

  character(LEN=80)::filename

  integer::ilun,idim,i
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::ii
  logical,allocatable,dimension(:)::nb
  integer,parameter::tag=1135
  integer::dummy_io,info2

  if(.not. sink) return

  if(verbose)write(*,*)'Entering backup_sink'

  ilun=4*ncpu+myid+10

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  open(unit=ilun,file=TRIM(fileloc),form='unformatted')
  rewind(ilun)

  write(ilun)nsink
  write(ilun)nindsink
  if(nsink>0)then
     allocate(xdp(1:nsink))
     do i=1,nsink
        xdp(i)=msink(i)
     end do
     write(ilun)xdp ! Write sink mass
     do i=1,nsink
        xdp(i)=tsink(i)
     end do


     !to be reintroduced
     !!ADDED by PH 09/2013
!     do i=1,nsink
!        xdp(i)=dmfsink(i)
!     end do
!     write(ilun)xdp ! Write sink mass



     write(ilun)xdp ! Write sink birth epoch
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=xsink(i,idim)
        end do
        write(ilun)xdp ! Write sink position
     enddo
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=vsink(i,idim)
        end do
        write(ilun)xdp ! Write sink velocity
     enddo
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=lsink(i,idim)
        end do
        write(ilun)xdp ! Write sink angular momentum
     enddo
     do i=1,nsink
        xdp(i)=delta_mass(i)
     end do
     write(ilun)xdp ! Write sink accumulated rest mass energy
     do i=1,nsink
        xdp(i)=acc_rate(i)
     end do
     write(ilun)xdp ! Write sink accretion rate


     !PH to be reintroduced
!     do i=1,nsink
!        xdp(i)=Eioni(i)
!     end do
!     write(ilun)xdp ! Write sink accumulated ionising photons

     do i=1,nsink
        xdp(i)=Teff_sink(i)
     end do
     write(ilun)xdp ! Write sink stellar effective temperature
     do i=1,nsink
        xdp(i)=rsink_star(i)
     end do
     write(ilun)xdp ! Write sink stellar radius
     deallocate(xdp)
     allocate(ii(1:nsink))
     do i=1,nsink
        ii(i)=idsink(i)
     end do
     write(ilun)ii ! Write sink index
     !        do i=1,nsink
     !           ii(i)=level_sink(i)
     !        end do
     !        write(ilun)ii ! Write sink level
     deallocate(ii)
     !        write(ilun)ncloud_sink   ! Write ncloud
     allocate(nb(1:nsink))
     do i=1,nsink
        nb(i)=new_born(i)
     end do
     write(ilun)nb ! Write level at which sinks where integrated
     deallocate(nb)
     write(ilun)sinkint_level ! Write level at which sinks where integrated
  endif
  close(ilun)
 
  ! Send the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
        dummy_io=1
        call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
             & MPI_COMM_WORLD,info2)
     end if
  endif
#endif
  
end subroutine backup_sink



subroutine output_sink(filename)
  use amr_commons
  use hydro_commons
  use pm_commons
  use units_commons
  use cloud_module
  implicit none
  character(LEN=80)::filename

  integer::i,idim,ipart,isink
  integer::nx_loc,ny_loc,nz_loc,ilun,icpu,idom
  real(dp)::scale,l_abs,rot_period,dx_min
!  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::star_mass,pi
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  
  pi=acos(-1.0d0)

  if(verbose)write(*,*)'Entering output_sink'

  ilun=myid+10

  ! Conversion factor from user units to cgs units                                                                   
!  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp
  
  if(verbose)write(*,*)'Entering output_sink'
  
  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace')
  !======================
  ! Write sink properties
  !======================
  write(ilun,*)'Number of sink = ',nsink

  write(ilun,'(" ======================================================================================================================================================================================= ")')
  write(ilun,'("  Id     M[Msol]    x           y           z           vx       vy       vz     rot_period[y] lx/|l|  ly/|l|  lz/|l| acc_rate[Msol/y] acc_lum[Lsol]  age[y]  int_lum[Lsol]     Teff [K] ")')
  write(ilun,'(" ======================================================================================================================================================================================= ")')  
  
  do isink=1,nsink
     star_mass=facc_star*msink(isink)
     l_abs=max((lsink(isink,1)**2+lsink(isink,2)**2+lsink(isink,3)**2)**0.5,1.d-50)
     rot_period=32*pi*star_mass*(dx_min)**2/(5*l_abs+tiny(0.d0))
     write(ilun,'(I5,2X,F9.5,3(2X,F10.7),3(2X,F7.4),2X,E13.3,3(2X,F6.3),5(2X,E11.3))')&
          idsink(isink),msink(isink)*scale_m/Msun, &
          xsink(isink,1:ndim),vsink(isink,1:ndim),&
          rot_period*scale_t/year,lsink(isink,1)/l_abs,lsink(isink,2)/l_abs,lsink(isink,3)/l_abs,&
          acc_rate(isink)*scale_m/Msun/(scale_t)*year,acc_lum(isink)/scale_t**2*scale_l**3*scale_d*scale_l**2/scale_t/Lsun,&
          (t-tsink(isink))*scale_t/year,&
          int_lum(isink)*scale_d*scale_l**3*scale_v**2/scale_t/Lsun,Teff_sink(isink)
  end do
  write(ilun,'(" ====================================================================================================================================================================================== ")')
  close(ilun)

end subroutine output_sink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_sink_csv(filename)
  use amr_commons
  use pm_commons
  use hydro_commons
  use units_commons
  use cloud_module
  implicit none
  character(LEN=80)::filename,fileloc

  integer::ilun,icpu,isink
  integer::i,idim,ipart
  integer::nx_loc,ny_loc,nz_loc
  real(dp)::scale,l_abs,star_mass,rot_period,dx_min,pi
  character(LEN=5)::nchar

  if(verbose)write(*,*)'Entering output_sink_csv'

  pi=acos(-1.0d0)

  ilun=2*ncpu+myid+10

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**nlevelmax/aexp

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace', recl=500)
  !======================
  ! Write sink properties
  !======================
  do isink=1,nsink
     star_mass=facc_star*msink(isink)
     l_abs=max((lsink(isink,1)**2+lsink(isink,2)**2+lsink(isink,3)**2)**0.5,1.d-50)
     rot_period=32*pi*star_mass*(dx_min)**2/(5*l_abs+tiny(0.d0))
     write(ilun,'(I10,16(A1,ES20.10))')&
          idsink(isink),',', &
          msink(isink)*scale_m/Msun,',',&
          xsink(isink,1),',',xsink(isink,2),',',xsink(isink,3),',',&
          vsink(isink,1),',',vsink(isink,2),',',vsink(isink,3),',',&
          rot_period*scale_t/year,',',&
          lsink(isink,1)/l_abs,',',lsink(isink,2)/l_abs,',',lsink(isink,3)/l_abs,',',&
          acc_rate(isink)*scale_m/Msun/(scale_t)*year,',',&
          acc_lum(isink)/scale_t**2*scale_l**3*scale_d*scale_l**2/scale_t/Lsun,',',&
          (t-tsink(isink))*scale_t/year,',',&
          int_lum(isink)*scale_d*scale_l**3*scale_v**2/scale_t/Lsun,',',&
          Teff_sink(isink)
  end do

  close(ilun)

end subroutine output_sink_csv
