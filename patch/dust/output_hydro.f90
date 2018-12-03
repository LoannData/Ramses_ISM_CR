subroutine file_descriptor_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif

  character(LEN=80)::filename
  character(LEN=80)::fileloc
  integer::ivar,ilun,idust
  
  if(verbose)write(*,*)'Entering file_descriptor_hydro'

  ilun=11

  ! Open file
  fileloc=TRIM(filename)
  open(unit=ilun,file=fileloc,form='formatted')

  ! Write run parameters
  write(ilun,'("nvar        =",I11)')nvar+4+ndust*ndim
  ivar=1
  write(ilun,'("variable #",I2,": density")')ivar
  if(write_conservative) then
     ivar=2
     write(ilun,'("variable #",I2,": momentum_x")')ivar
     ivar=3
     write(ilun,'("variable #",I2,": momentum_y")')ivar
     ivar=4
     write(ilun,'("variable #",I2,": momentum_z")')ivar
  else
     ivar=2
     write(ilun,'("variable #",I2,": velocity_x")')ivar
     ivar=3
     write(ilun,'("variable #",I2,": velocity_y")')ivar
     ivar=4
     write(ilun,'("variable #",I2,": velocity_z")')ivar
  endif
  ivar=5
  write(ilun,'("variable #",I2,": B_left_x")')ivar
  ivar=6
  write(ilun,'("variable #",I2,": B_left_y")')ivar
  ivar=7
  write(ilun,'("variable #",I2,": B_left_z")')ivar
  ivar=8
  write(ilun,'("variable #",I2,": B_right_x")')ivar
  ivar=9
  write(ilun,'("variable #",I2,": B_right_y")')ivar
  ivar=10
  write(ilun,'("variable #",I2,": B_right_z")')ivar
#if NENER>NGRP
  if(write_conservative) then
#if NCR>0  
     ! CR energies
     do ivar=1,ncr
        write(ilun,'("variable #",I2,": cosmic_rays_energy_",I1)')10+ivar,ivar
     end do
#endif
     ! Non-thermal energies
     do ivar=1+ncr,nent
        write(ilun,'("variable #",I2,": non_thermal_energy_",I1)')10+ivar,ivar
     end do
  else
#if NCR>0
     ! CR pressures
     do ivar=1,ncr
        write(ilun,'("variable #",I2,": cosmic_rays_pressure_",I1)')10+ivar,ivar
     end do
#endif
     ! Non-thermal pressures
     do ivar=1+ncr,nent
        write(ilun,'("variable #",I2,": non_thermal_pressure_",I1)')10+ivar,ivar
     end do
  end if
#endif
  if(write_conservative) then
     ivar=11+nent
     write(ilun,'("variable #",I2,": total_energy")')ivar
  else
     ivar=11+nent
     write(ilun,'("variable #",I2,": thermal_pressure")')ivar
  endif
#if NGRP>0
  ! Radiative energies
  do ivar=1,ngrp
     write(ilun,'("variable #",I2,": radiative_energy_",I1)')firstindex_er+3+ivar,ivar
  end do
#endif
#if USE_M_1==1
  ! Radiative fluxes
  do ivar=1,ngrp
     write(ilun,'("variable #",I2,": radiative_flux_x",I1)')firstindex_fr+3       +ivar,ivar
  end do
if(ndim>1) then
  do ivar=1,ngrp
     write(ilun,'("variable #",I2,": radiative_flux_y",I1)')firstindex_fr+3+  ngrp+ivar,ivar
  end do
endif
if(ndim>2) then
  do ivar=1,ngrp
     write(ilun,'("variable #",I2,": radiative_flux_z",I1)')firstindex_fr+3+2*ngrp+ivar,ivar
  end do
endif
#endif
#if NEXTINCT>0
  ! Extinction
  do ivar=1,nextinct
     write(ilun,'("variable #",I2,": extinction",I1)')firstindex_extinct+3+ivar,ivar
  end do
#endif
#if NPSCAL>0
  if(write_conservative) then
     ! Passive scalars
     do ivar=1,npscal
        write(ilun,'("variable #",I2,": passive_scalar_cons_",I2)')firstindex_pscal+3+ivar,ivar
     end do
  else
#if NDUST>0
     ! dust ratio is separated from pscal to avoid confusion while reading the files
     do ivar=firstindex_pscal+1,firstindex_ndust
        write(ilun,'("variable #",I2,": passive_scalar_",I1)')3+ivar,ivar-firstindex_pscal
     end do
     do ivar=firstindex_ndust+1,firstindex_ndust+ndust
        if(ivar-firstindex_ndust<10) write(ilun,'("variable #",I2,": epsilon_",I1)')3+ivar,ivar-firstindex_ndust
        if(ivar-firstindex_ndust.ge.10) write(ilun,'("variable #",I2,": epsilon_",I2)')3+ivar,ivar-firstindex_ndust

     end do
     do ivar=firstindex_ndust+ndust+1,firstindex_pscal+npscal
        write(ilun,'("variable #",I2,": passive_scalar_",I1)')3+ivar,ivar-firstindex_ndust-ndust
     end do
#else
  ! Passive scalars
     do ivar=1,npscal
        write(ilun,'("variable #",I2,": passive_scalar_",I1)')firstindex_pscal+3+ivar,ivar
     end do
#endif
  endif

#endif
  ! Temperature
  ivar=firstindex_pscal+3+npscal+1
  write(ilun,'("variable #",I2,": temperature")')ivar
#if NDUST>0
  do idust=1,ndust
     if(idust<10)write(ilun,'("variable #",I2,":vdx_",I1)')ivar+idust,idust
     if(idust.ge.10)write(ilun,'("variable #",I2,":vdx_",I2)')ivar+idust,idust

#if NDIM>1
     if(idust<10)write(ilun,'("variable #",I2,":vdy_",I1)')ivar+idust+1,idust
     if(idust.ge.10)write(ilun,'("variable #",I2,":vdy_",I2)')ivar+idust+1,idust

#endif
#if NDIM>2
     if(idust<10) write(ilun,'("variable #",I2,":vdz_",I1)')ivar+idust+2,idust
     if(idust.ge.10) write(ilun,'("variable #",I2,":vdz_",I2)')ivar+idust+2,idust

#endif          
     ivar = ivar +(ndim-1)
  end do
#endif   
  close(ilun)

end subroutine file_descriptor_hydro

subroutine backup_hydro(filename)
  use amr_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::dummy_io,info2
#endif
  character(LEN=80)::filename
  integer::i,ivar,ncache,ind,ilevel,igrid,iskip,ilun,istart,ibound,ht,idim
  real(dp)::d,u,v,w,A,B,C,e
  real(dp):: sum_dust
#if NDUST>0
  integer:: idust
#endif  
  integer,allocatable,dimension(:)::ind_grid
  real(dp)::cmp_temp,p
  real(dp),allocatable,dimension(:)::xdp
  character(LEN=5)::nchar
  character(LEN=80)::fileloc
  integer,parameter::tag=1121
#if NENER>0
  integer::irad
#endif

  if(verbose)write(*,*)'Entering backup_hydro'

  ilun=ncpu+myid+10

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
  open(unit=ilun,file=fileloc,form='unformatted')
  write(ilun)ncpu
!   if(eos) then 
!      write(ilun)nvar+3+1
!   else
!      write(ilun)nvar+3
!   endif
  write(ilun)nvar+4
  write(ilun)ndim
  write(ilun)nlevelmax
  write(ilun)nboundary
  write(ilun)gamma
  do ilevel=1,nlevelmax
     do ibound=1,nboundary+ncpu
        if(ibound<=ncpu)then
           ncache=numbl(ibound,ilevel)
           istart=headl(ibound,ilevel)
        else
           ncache=numbb(ibound-ncpu,ilevel)
           istart=headb(ibound-ncpu,ilevel)
        end if
        write(ilun)ilevel
        write(ilun)ncache
        if(ncache>0)then
           allocate(ind_grid(1:ncache),xdp(1:ncache))
           ! Loop over level grids
           igrid=istart
           do i=1,ncache
              ind_grid(i)=igrid
              igrid=next(igrid)
           end do
           ! Loop over cells
           do ind=1,twotondim
              iskip=ncoarse+(ind-1)*ngridmax
              do ivar=1,4
                 if(ivar==1)then ! Write density
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,1)
                    end do
                 else ! Write velocity field
                    if(write_conservative) then
                       do i=1,ncache
                          xdp(i)=uold(ind_grid(i)+iskip,ivar)
                       end do
                    else
                       do i=1,ncache
                          xdp(i)=uold(ind_grid(i)+iskip,ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                       end do
                    endif
                 endif
                 write(ilun)xdp
              end do
              do ivar=6,8 ! Write left B field
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
              do ivar=nvar+1,nvar+3 ! Write right B field
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,ivar)
                 end do
                 write(ilun)xdp
              end do
#if NENER>NGRP
              ! Write non-thermal pressures
              if(write_conservative) then
                 do ivar=1,nent
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,8+ivar)
                    end do
                    write(ilun)xdp
                 end do
              else
                 do ivar=1,nent
                    do i=1,ncache
                       xdp(i)=(gamma_rad(ivar)-1d0)*uold(ind_grid(i)+iskip,8+ivar)
                    end do
                    write(ilun)xdp
                 end do
              endif
#endif
              if(write_conservative) then
                 do i=1,ncache ! Write total energy
                    xdp(i)=uold(ind_grid(i)+iskip,5)
                 enddo
                 write(ilun)xdp
              else
                 do i=1,ncache ! Write thermal pressure
                    d=max(uold(ind_grid(i)+iskip,1),smallr)
                    u=uold(ind_grid(i)+iskip,2)/d
                    v=uold(ind_grid(i)+iskip,3)/d
                    w=uold(ind_grid(i)+iskip,4)/d
                    A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                    B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                    C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                    e=uold(ind_grid(i)+iskip,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
                    do irad=1,nener
                       e=e-uold(ind_grid(i)+iskip,8+irad)
                    end do
#endif
                    sum_dust=0.0d0
#if NDUST>0
                    do idust=1,ndust
                       sum_dust=sum_dust+uold(ind_grid(i)+iskip,firstindex_ndust+idust)/d
                    end do
#endif                    
                    call pressure_eos((1.0d0-sum_dust)*d,e,p)
                    xdp(i)=p
                 end do
                 write(ilun)xdp
              endif

#if NGRP>0
              do ivar=1,ngrp ! Write radiative energy if any
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,firstindex_er+ivar)
                 end do
                 write(ilun)xdp
              end do
#if USE_M_1==1
              do ivar=1,nfr ! Write radiative flux if any
                 do i=1,ncache
                    xdp(i)=uold(ind_grid(i)+iskip,firstindex_fr+ivar)
                 end do
                 write(ilun)xdp
              end do
#endif
#endif
#if NEXTINCT>0
              ! Write extinction if activated
              do i=1,ncache
                 xdp(i)=uold(ind_grid(i)+iskip,firstindex_extinct+1)
              end do
              write(ilun)xdp
#endif

#if NPSCAL>0
              if(write_conservative) then
                 do ivar=1,npscal ! Write conservative passive scalars if any
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+ivar)
                    end do
                    write(ilun)xdp
                 end do
              else
                 do ivar=1,npscal ! Write passive scalars if any
                    do i=1,ncache
                       xdp(i)=uold(ind_grid(i)+iskip,firstindex_pscal+ivar)/max(uold(ind_grid(i)+iskip,1),smallr)
                    end do
                    write(ilun)xdp
                 end do
              endif
#endif
              
              ! Write temperature
              do i=1,ncache
                 d=max(uold(ind_grid(i)+iskip,1),smallr)
                 if(energy_fix) then
                    e=uold(ind_grid(i)+iskip,nvar)
                 else
                    u=uold(ind_grid(i)+iskip,2)/d
                    v=uold(ind_grid(i)+iskip,3)/d
                    w=uold(ind_grid(i)+iskip,4)/d
                    A=0.5*(uold(ind_grid(i)+iskip,6)+uold(ind_grid(i)+iskip,nvar+1))
                    B=0.5*(uold(ind_grid(i)+iskip,7)+uold(ind_grid(i)+iskip,nvar+2))
                    C=0.5*(uold(ind_grid(i)+iskip,8)+uold(ind_grid(i)+iskip,nvar+3))
                    e=uold(ind_grid(i)+iskip,5)-0.5*d*(u**2+v**2+w**2)-0.5*(A**2+B**2+C**2)
#if NENER>0
                    do irad=1,nener
                       e=e-uold(ind_grid(i)+iskip,8+irad)
                    end do
#endif
                 endif
                 sum_dust=0.0d0
#if NDUST>0
                 do idust=1,ndust
                    sum_dust=sum_dust+uold(ind_grid(i)+iskip,firstindex_ndust+idust)/d
                 end do
#endif                    
                 call temperature_eos((1.0d0-sum_dust)*d,e,cmp_temp,ht,sum_dust)
                 xdp(i)=cmp_temp
             
              end do
              write(ilun)xdp   
#if NDUST>0
           do idust=1,ndust
              do idim=1,ndim           
                 do i=1,ncache
                    xdp(i)=v_dust(ind_grid(i)+iskip,idust,idim)!/(1.0d0-sum_dust)
                 end do
                 write(ilun)xdp
              end do
           end do
#endif
        end do
        
        deallocate(ind_grid, xdp)
        end if
     end do
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

end subroutine backup_hydro





