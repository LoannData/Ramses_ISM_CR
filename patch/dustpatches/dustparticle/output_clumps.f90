! Output clumps every main step
! Use to construct clump evolution / merger trees
! Sam Geen, November 2016

subroutine write_clumps_each_step
  use amr_commons
  use clfind_commons
  use pm_commons
  implicit none
  integer::info,stat
  logical::folder_exists
  character(LEN=5)::soutput,scpu,nchar
  character(LEN=80)::clumpfilename,sinkfilename,filedir,readmefile,filecmd
  if(npeaks.gt.0)then
     if(verbose)write(*,*)'Entering write_clump_each_step'
     ! Check whether we need to reset the current file
     call title(ifout-1,nchar)
     ! Timestep number
     call title(myid,scpu)
     ! Shouldn't need to make the folder, but just in case
     filedir='output_'//TRIM(nchar)//'/'
     inquire(file=readmefile,exist=folder_exists)
     ! Make folder?
     if(.not.folder_exists) then
        filecmd='mkdir -p '//TRIM(filedir)
#ifdef NOSYSTEM
        call PXFMKDIR(TRIM(filedir),LEN(TRIM(filedir)),O'755',info)
#else
        call system(filecmd)
#endif
     endif
     ! Make file names
     clumpfilename=TRIM(filedir)//'clump_'//TRIM(nchar)//'.dat'//TRIM(scpu)
     sinkfilename=TRIM(filedir)//'sink_'//TRIM(nchar)//'.dat'//TRIM(scpu)
     ! Check for old file to wipe in case of restart
     if (clump_file_id.ne.ifout-1) then
        clump_file_id = ifout-1
        ! Delete old files from a previous run before writing to them
        open(unit=9099135, iostat=stat, file=clumpfilename, status='old')
        if (stat == 0) close(9099135, status='delete')
        open(unit=9099135, iostat=stat, file=sinkfilename, status='old')
        if (stat == 0) close(9099135, status='delete')
     endif
     ! Make readme file
     readmefile='output_clumps.txt'
     ! Make folder
     ! Write readme file (to allow checking for folder existing)
     open(unit=13337,file=TRIM(readmefile),form='formatted',status='replace', recl=500)
     write(13337,*) "These files store clump and sink information for each timestep"
     write(13337,*) "These are used for making clump merger trees"
     write(13337,*) "clump_?????.dat format:"
     write(13337,*)"index t lev parent ncell peak_x peak_y peak_z v_x v_y v_z "//&
          "rho- rho+ rho_av mass_cl relevance"
     write(13337,*) "sink_?????.dat format:"
     write(13337,*) "idsink,t,msink,x,y,z,vx,vy,vz,lx,ly,lz,age,accretionrate,eioni"
     close(13337)
     ! Write sink CSV
     call output_sink_properties(sinkfilename)
     ! Write clump CSV
     call output_clump_properties(clumpfilename)
  endif

end subroutine write_clumps_each_step
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_clump_properties(fileloc)
  use amr_commons
  use pm_commons,ONLY:mp
  use hydro_commons,ONLY:mass_sph
  use clfind_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer,parameter::tag=1101
  integer::dummy_io,info2
  !---------------------------------------------------------------------------
  ! this routine writes the clump properties to screen and to file
  !---------------------------------------------------------------------------

  logical::activepeaks
  integer::i,j,jj,ilun,n_rel,n_rel_tot,info,nx_loc
  real(dp)::rel_mass,rel_mass_tot,scale,particle_mass,particle_mass_tot
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu
  real(dp),dimension(1:npeaks)::peakd
  integer,dimension(1:npeaks)::ind_sort

  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  if(ivar_clump==0)then
     particle_mass=MINVAL(mp, MASK=(mp.GT.0.))
#ifndef WITHOUTMPI  
     call MPI_ALLREDUCE(particle_mass,particle_mass_tot,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,info)
     particle_mass=particle_mass_tot  
#endif
  else
     if(hydro)then
        particle_mass=mass_sph
     endif
  endif

  ! sort clumps by peak density in ascending order
  do i=1,npeaks
     peakd(i)=max_dens(i)
     ind_sort(i)=i
  end do
  call quick_sort_dp(peakd,ind_sort,npeaks) 

  ! Check for peaks in this file
  activepeaks = .false.
  do j=npeaks,1,-1
     jj=ind_sort(j)
     if (relevance(jj) > relevance_threshold .and. halo_mass(jj) > mass_threshold*particle_mass)then           
        activepeaks = .true.
        exit ! Peak found, all ok
     endif
  enddo
  if (.not.activepeaks) return

  ilun=20

  ! print results in descending order to screen/file
  rel_mass=0.
  n_rel=0
  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif
  call title(ifout-1,nchar)

  open(unit=ilun,file=fileloc,form='formatted',position="append")
  
  do j=npeaks,1,-1
     jj=ind_sort(j)
     if (relevance(jj) > relevance_threshold .and. halo_mass(jj) > mass_threshold*particle_mass)then           
        write(ilun,'(I8,X,1PE18.9E2,X,I2,X,I10,X,I10,14(X,1PE18.9E2))')&
             jj+ipeak_start(myid)&
             ,t&
             ,lev_peak(jj)&
             ,new_peak(jj)&
             ,n_cells(jj)&
             ,peak_pos(jj,1)&
             ,peak_pos(jj,2)&
             ,peak_pos(jj,3)&
             ,center_of_mass(jj,1)&
             ,center_of_mass(jj,2)&
             ,center_of_mass(jj,3)&
             ,clump_velocity(jj,1)&
             ,clump_velocity(jj,2)&
             ,clump_velocity(jj,3)&
             ,min_dens(jj)&
             ,max_dens(jj)&
             ,clump_mass(jj)/clump_vol(jj)&
             ,clump_mass(jj)&
             ,relevance(jj)
        rel_mass=rel_mass+clump_mass(jj)
        n_rel=n_rel+1
     end if
  end do

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

end subroutine output_clump_properties

subroutine output_sink_properties(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename,fileloc

  integer::ilun,icpu,isink


  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',position="append")
  !======================
  ! Write sink properties
  !======================
  do isink=1,nsink
     write(ilun,'(I10,14(A1,ES20.10))')idsink(isink),',',t,',',msink(isink),&
          ',',xsink(isink,1),',',xsink(isink,2),',',xsink(isink,3),&
          ',',vsink(isink,1),',',vsink(isink,2),',',vsink(isink,3),&
          ',',lsink(isink,1),',',lsink(isink,2),',',lsink(isink,3),&
          ',',t-tsink(isink),',',dMBHoverdt(isink), ',', Eioni(isink)
  end do

  close(ilun)

end subroutine output_sink_properties
