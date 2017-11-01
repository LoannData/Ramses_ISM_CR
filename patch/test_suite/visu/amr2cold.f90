program amr2cube
  !--------------------------------------------------------------------------
  ! Ce programme calcule le cube cartesien pour les
  ! variables hydro d'une simulation RAMSES. 
  ! Version F90 par R. Teyssier le 01/04/01.
  ! Retourne les data au format vtk (vecteur pour vitesse et champ magnetique)
  ! Retourne les data au format radmc:
  !      - amr_grid.inp : contient la structure de la grille regulier 
  !        (coordonnees en cgs)
  !      - dust_density.inp : contient la densite (cgs) des poussieres   
  !        (typiquement 1/100 de la densite du gaz)
  !      - dust_temperature.inp : contient la temperature (K) des poussieres 
  !--------------------------------------------------------------------------
  !--------------------------------------------------------------------------
  !  En plus de l'original : format gnuplot de sortie. Unités changées : unités
  !  astronomique, densité en particules par cc.
  !--------------------------------------------------------------------------




  implicit none
  integer::ndim,n,i,j,k,l,twotondim,ncoarse,type=0,domax=0
  integer::ivar,nvar,ncpu,ncpuh,lmax=0,nboundary,ngrid_current
  integer::nx,ny,nz,ilevel,idim,jdim,kdim,icell
  integer::nlevelmax,ilevel1,ngrid1
  integer::nlevelmaxs,nlevel,iout
  integer::ind,ipos,ngrida,ngridh,ilevela,ilevelh
  integer::ngridmax,nstep_coarse,icpu,ncpu_read
  integer::nhx,nhy,ihx,ihy,ivar1,ivar2,igrida
  real::gamma,smallr,smallc,gammah
  real::boxlen,boxlen2
  real::t,aexp,hexp,t2,aexp2,hexp2
  real::omega_m,omega_l,omega_k,omega_b
  real::omega_m2,omega_l2,omega_k2,omega_b2,ud,ul,ut,utemp,mu_gas,maxtp

  integer::nx_sample=0,ny_sample=0,nz_sample=0
  integer::ngrid,imin,imax,jmin,jmax,kmin,kmax
  integer::ncpu2,npart2,ndim2,nlevelmax2,nstep_coarse2
  integer::nx2,ny2,nz2,ngridmax2,nvarh,ndimh,nlevelmaxh
  integer::nx_full,ny_full,nz_full,lmin,levelmin
  integer::ix,iy,iz,ixp1,iyp1,izp1,ndom,impi,bit_length,maxdom,ixx,iyy,izz
  integer,dimension(1:8)::idom,jdom,kdom,cpu_min,cpu_max
  real(KIND=8),dimension(1:8)::bounding_min,bounding_max
  real(KIND=8)::dkey,order_min,dmax,dummy
  real(KIND=8)::xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1
  real(KIND=8)::xxmin,xxmax,yymin,yymax,zzmin,zzmax
  real(KIND=8)::ddx,ddy,ddz,dex,dey,dez,xx,yy,zz,pos,dx,dy,dz,eps
  real(KIND=8),dimension(:,:),allocatable::x,xg,x_tilde
  real(KIND=8),dimension(:,:,:),allocatable::var
  real(KIND=4),dimension(:,:,:),allocatable::toto
  real(KIND=8),dimension(:)  ,allocatable::rho
  real(KIND=8),dimension(:,:)  ,allocatable::vec

  logical,dimension(:)  ,allocatable::ref
  integer,dimension(:)  ,allocatable::isp
  integer,dimension(:,:),allocatable::son,ngridfile,ngridlevel,ngridbound
  real(KIND=8),dimension(1:8,1:3)::xc
  real(KIND=8),dimension(1:3)::xbound=(/0d0,0d0,0d0/)
  character(LEN=5)::nchar,ncharcpu
  character(LEN=80)::ordering
  character(LEN=128)::nomfich,repository,outfich,filetype='bin'
  logical::ok,ok_part,ok_cell
  real(KIND=8),dimension(:),allocatable::bound_key
  logical,dimension(:),allocatable::cpu_read
  integer,dimension(:),allocatable::cpu_list
  character(LEN=1)::proj='z'
  real(KIND=4) :: norm
  ! nouvelles variables de normalisation
  real(KIND=8)   :: scale_t_kyrs,scale_ua,jnorm
  integer      :: old, direction
  real(KIND=8) :: tempmin,tempmax 
  real(KIND=8),dimension(3) :: center,J10,J100,J500,J1000,z_tilde,y_tilde
  real(KIND=8),dimension(:,:),allocatable :: rvectJ
  real(KIND=8),dimension(:,:,:),allocatable :: e_tilde
  integer :: corr
  character(len=20) :: str

  type level
     integer::ilevel
     integer::ngrid
     real(KIND=4),dimension(:,:,:),pointer::cube
     real(KIND=4),dimension(:,:,:,:),pointer::cubevec
     integer::imin
     integer::imax
     integer::jmin
     integer::jmax
     integer::kmin
     integer::kmax
  end type level

  type(level),dimension(1:100)::grid

  ! Temporary space for reading labels from the info file.
  character(LEN=128)::temp_label

  call read_params

  !-----------------------------------------------
  ! Normalisation
  !-----------------------------------------------

  J100 = (/0.0,0.0,1.0/)
  
!   open(unit=1234,file="angmom.dat",STATUS="OLD")
!   read(1234,*)J10, J100,J500,J1000
!   close(1234)
!   if ( abs(J100(2)) .ge. max(abs(J100(1)),abs(J100(3))) ) then
!      corr=2
!   else if ( abs(J100(2)) .ge. max(abs(J100(1)),abs(J100(3))) ) then
!      corr=1
!   else
!      corr=0
!   end if


  !-----------------------------------------------
  ! Lecture du fichier hydro au format RAMSES
  !-----------------------------------------------
  !print*, repository
  ipos=INDEX(repository,'output_')
  nchar=repository(ipos+7:ipos+13)
  nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     !print*, "la"
     stop
  endif
  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  inquire(file=nomfich, exist=ok) ! verify input file 
  if ( .not. ok ) then
     print *,TRIM(nomfich)//' not found.'
     stop
  endif

  nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
  open(unit=10,file=nomfich,status='old',form='unformatted')
  read(10)ncpu
  read(10)ndim
  read(10)nx,ny,nz
  read(10)nlevelmax
  read(10)ngridmax
  read(10)nboundary
  read(10)ngrid_current
  read(10)boxlen
  close(10)
  twotondim=2**ndim
  xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)

  allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
  allocate(ngridlevel(1:ncpu,1:nlevelmax))
  if(nboundary>0)allocate(ngridbound(1:nboundary,1:nlevelmax))

  if(ndim==2)then
     write(*,*)'Output file contains 2D data'
     write(*,*)'Aborting'
     stop
  endif

  nomfich=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  open(unit=10,file=nomfich,form='formatted',status='old')
  read(10,*)
  read(10,*)
  read(10,'(A13,I11)')temp_label,levelmin
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)

  read(10,'(A13,E23.15)')temp_label,boxlen
  read(10,'(A13,E23.15)')temp_label,t
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,'(A13,E23.15)')temp_label,ul
  read(10,'(A13,E23.15)')temp_label,ud
  read(10,'(A13,E23.15)')temp_label,ut
  read(10,'(A13,E23.15)')temp_label,mu_gas
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  read(10,*)
  !print*,'ud',ud,ul,ut
  scale_t_kyrs = ut/(3600.*24*365.25*1.e3) ! convert code units to kyrs
  scale_ua = ul*6.6846e-14 ! convert code units to cm (ul) and then to AU
  utemp=mu_gas*1.66d-24/1.38d-16
  utemp=utemp*ul*ul/ut/ut
  !print*,utemp
  read(10,'(A14,A80)')temp_label,ordering
  write(*,'(XA14,A20)')temp_label,TRIM(ordering)
  read(10,*)
  allocate(cpu_list(1:ncpu))
  if(TRIM(ordering).eq.'hilbert')then
     allocate(bound_key(0:ncpu))
     allocate(cpu_read(1:ncpu))
     cpu_read=.false.
     do impi=1,ncpu
        read(10,'(I8,1X,E23.15,1X,E23.15)')i,bound_key(impi-1),bound_key(impi)
     end do
  endif
  close(10)

  !-----------------------
  ! Map parameters
  !-----------------------
  if(lmax==0)then
     lmax=nlevelmax
  endif
  write(*,*)'time=',t
  write(*,*)'Working resolution =',2**lmax
!center(1)=0.49645996
!center(2)=0.4998779
!center(3)=0.5004883
  xmin=xmin+(center(1)-0.5)
  xmax=xmax+(center(1)-0.5)
  ymin=ymin+(center(2)-0.5)
  ymax=ymax+(center(2)-0.5)
  zmin=zmin+(center(3)-0.5)
  zmax=zmax+(center(3)-0.5)
!if  ( direction == 1 ) then     ! axe x
!    ! rien à faire
!else if ( direction == 2 ) then   ! axe y : xx sera y, yy sera z, zz sera x
!    tempmin = xmin
!    tempmax = xmax
!    xmin = ymin
!    xmax = ymax
!    ymin = zmin
!    ymax = zmax
!    zmin = tempmin
!    zmax = tempmax
!else if ( direction == 3 ) then  ! axe z : permutation circulaire 
!    tempmin = xmin
!    tempmax = xmax
!    xmin = zmin
!    xmax = zmax
!    zmin = ymin 
!    zmax = ymax
!    ymin = tempmin
!    ymax = tempmax
!end if
  xxmin=xmin ; xxmax=xmax
  yymin=ymin ; yymax=ymax
  zzmin=zmin ; zzmax=zmax

  if(TRIM(ordering).eq.'hilbert')then

     dmax=max(xxmax-xxmin,yymax-yymin,zzmax-zzmin)
     do ilevel=1,lmax
        dx=0.5d0**ilevel
        if(dx.lt.dmax)exit
     end do
     lmin=ilevel
     bit_length=lmin-1
     maxdom=2**bit_length
     imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0
     if(bit_length>0)then
        imin=int(xxmin*dble(maxdom))
        imax=imin+1
        jmin=int(yymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zzmin*dble(maxdom))
        kmax=kmin+1
     endif

     dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
     ndom=1
     if(bit_length>0)ndom=8
     idom(1)=imin; idom(2)=imax
     idom(3)=imin; idom(4)=imax
     idom(5)=imin; idom(6)=imax
     idom(7)=imin; idom(8)=imax
     jdom(1)=jmin; jdom(2)=jmin
     jdom(3)=jmax; jdom(4)=jmax
     jdom(5)=jmin; jdom(6)=jmin
     jdom(7)=jmax; jdom(8)=jmax
     kdom(1)=kmin; kdom(2)=kmin
     kdom(3)=kmin; kdom(4)=kmin
     kdom(5)=kmax; kdom(6)=kmax
     kdom(7)=kmax; kdom(8)=kmax
     
     do i=1,ndom
        if(bit_length>0)then
           call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
           order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
     end do
     
     cpu_min=0; cpu_max=0
     do impi=1,ncpu
        do i=1,ndom
           if (   bound_key(impi-1).le.bounding_min(i).and.&
                & bound_key(impi  ).gt.bounding_min(i))then
              cpu_min(i)=impi
           endif
           if (   bound_key(impi-1).lt.bounding_max(i).and.&
                & bound_key(impi  ).ge.bounding_max(i))then
              cpu_max(i)=impi
           endif
        end do
     end do
     
     ncpu_read=0
     do i=1,ndom
        do j=cpu_min(i),cpu_max(i)
           if(.not. cpu_read(j))then
              ncpu_read=ncpu_read+1
              cpu_list(ncpu_read)=j
              cpu_read(j)=.true.
           endif
        enddo
     enddo
  else
     ncpu_read=ncpu
     do j=1,ncpu
        cpu_list(j)=j
     end do
  end  if

  !-----------------------------
  ! Compute hierarchy
  !-----------------------------
  do ilevel=1,lmax
     nx_full=2**ilevel
     ny_full=2**ilevel
     nz_full=2**ilevel
     imin=int(xxmin*dble(nx_full))+1
     imax=int(xxmax*dble(nx_full))+1
     jmin=int(yymin*dble(ny_full))+1
     jmax=int(yymax*dble(ny_full))+1
     kmin=int(zzmin*dble(nz_full))+1
     kmax=int(zzmax*dble(nz_full))+1
if (direction ==1) then
     allocate(grid(ilevel)%cube(1,jmin:jmax,kmin:kmax))
     allocate(grid(ilevel)%cubevec(1,jmin:jmax,kmin:kmax,1:3))
else if (direction ==2) then
     allocate(grid(ilevel)%cube(1,imin:imax,kmin:kmax))
     allocate(grid(ilevel)%cubevec(1,imin:imax,kmin:kmax,1:3))
else if (direction ==3) then
     allocate(grid(ilevel)%cube(1,jmin:jmax,imin:imax))
     allocate(grid(ilevel)%cubevec(1,jmin:jmax,imin:imax,1:3))
end if
     grid(ilevel)%cube=0.0
     grid(ilevel)%cubevec=0.0
     grid(ilevel)%imin=imin
     grid(ilevel)%imax=imax
     grid(ilevel)%jmin=jmin
     grid(ilevel)%jmax=jmax
     grid(ilevel)%kmin=kmin
     grid(ilevel)%kmax=kmax
  end do

  !-----------------------------------------------
  ! Compute projected variables
  !----------------------------------------------

  ! Loop over processor files
  do k=1,ncpu_read
     icpu=cpu_list(k)
     call title(icpu,ncharcpu)

     ! Open AMR file and skip header
     nomfich=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=10,file=nomfich,status='old',form='unformatted')
     write(*,*)'Processing file '//TRIM(nomfich)
     do i=1,21
        read(10)
     end do
     ! Read grid numbers
     read(10)ngridlevel
     ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
     read(10)
     if(nboundary>0)then
        do i=1,2
           read(10)
        end do
        read(10)ngridbound
        ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
     endif
     read(10)
! ROM: comment the single follwing line for old stuff
     read(10)
     if(TRIM(ordering).eq.'bisection')then
        do i=1,5
           read(10)
        end do
     else
        read(10)
     endif
     read(10)
     read(10)
     read(10)

     ! Open HYDRO file and skip header
     nomfich=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
     open(unit=11,file=nomfich,status='old',form='unformatted')
     read(11)
     read(11)nvarh
     read(11)
     read(11)
     read(11)
     read(11)

     ! Loop over levels
     do ilevel=1,lmax

        ! Geometry
        dx=0.5**ilevel
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel
        do ind=1,twotondim
           iz=(ind-1)/4
           iy=(ind-1-4*iz)/2
           ix=(ind-1-2*iy-4*iz)
           xc(ind,1)=(dble(ix)-0.5D0)*dx
           xc(ind,2)=(dble(iy)-0.5D0)*dx
           xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do

        ! Allocate work arrays
        ngrida=ngridfile(icpu,ilevel)
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
           allocate(xg(1:ngrida,1:ndim))
           allocate(son(1:ngrida,1:twotondim))
           allocate(var(1:ngrida,1:twotondim,1:nvarh))
           allocate(x  (1:ngrida,1:ndim))
           allocate(rho(1:ngrida))
           allocate(ref(1:ngrida))
           allocate(x_tilde  (1:ngrida,1:ndim))
           allocate(vec(1:ngrida,1:12))
           allocate(rvectJ(1:ngrida,1:3))
           allocate(e_tilde(1:ngrida,1:3,1:3))
        endif

        ! Loop over domains
        do j=1,nboundary+ncpu

           ! Read AMR data
           if(ngridfile(j,ilevel)>0)then
              read(10) ! Skip grid index
              read(10) ! Skip next index
              read(10) ! Skip prev index
              ! Read grid center
              do idim=1,ndim
                 if(j.eq.icpu)then
                    read(10)xg(:,idim)
                 else
                    read(10)
                 endif
              end do
              read(10) ! Skip father index
              do ind=1,2*ndim
                 read(10) ! Skip nbor index
              end do
              ! Read son index
              do ind=1,twotondim
                 if(j.eq.icpu)then
                    read(10)son(:,ind)
                 else
                    read(10)
                 end if
              end do
              ! Skip cpu map
              do ind=1,twotondim
                 read(10)
              end do
              ! Skip refinement map
              do ind=1,twotondim
                 read(10)
              end do
           endif

           ! Read HYDRO data
           read(11)
           read(11)
           if(ngridfile(j,ilevel)>0)then
              ! Read hydro variables
              do ind=1,twotondim
                 do ivar=1,nvarh
                    if(j.eq.icpu)then
                       read(11)var(:,ind,ivar)
                    else
                       read(11)
                    end if
                 end do
              end do
           end if
        end do

        ! Compute map
        if(ngrida>0)then

           ! Loop over cells
           do ind=1,twotondim
               rho=0.

              ! Compute cell center
              do i=1,ngrida
                 x(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
                 x(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
                 x(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
              end do
              ! Check if cell is refined
              do i=1,ngrida
                 ref(i)=son(i,ind)>0.and.ilevel<lmax
              end do
              ! Extract variable
              select case (type)
              case (-1)
                 rho = icpu
              case (0)
                 rho = ilevel
              case (1)

                 vec(:,5)=(dble((x(:,1)-center(1))**2+(x(:,2)-center(2))**2&
                 &+(x(:,3)-center(3))**2)**0.5)*scale_ua*boxlen
                 vec(:,4)=sqrt(0.25*((var(:,ind,5)+var(:,ind,8))**2 &
                 &         +(var(:,ind,6)+var(:,ind,9))**2 &
                 &         +(var(:,ind,7)+var(:,ind,10))**2) &
                 &         *4.*acos(-1.)*ud*(ul/ut)**2)
                 !                     rho(:) = var(:,ind,1)*ud
                 rho(:) = var(:,ind,1)*ud*1.d6 !rho en g/mcube
                 !!!!!!!!!!!!! tests pour des angles !
                 do  i=1,3
                     e_tilde(:,3,i) = J100(i)      ! ez tilde, deja normalise
                 end do
                 rvectJ(:,1) = (x(:,2)-center(2))*J100(3) - (x(:,3)-center(3))*J100(2)
                 rvectJ(:,2) = (x(:,3)-center(3))*J100(1) - (x(:,1)-center(1))*J100(3)
                 rvectJ(:,3) = (x(:,1)-center(1))*J100(2) - (x(:,2)-center(2))*J100(1)
                 do  i=1,3
                     e_tilde(:,2,i) = -rvectJ(:,i)/&
                     &(sqrt(rvectJ(:,1)**2+rvectJ(:,2)**2+rvectJ(:,3)**2)) ! -rvectJ/norm(rvectJ)
                 end do
                 e_tilde(:,1,1) = e_tilde(:,2,2)*e_tilde(:,3,3)&
                 &- e_tilde(:,2,3)*e_tilde(:,3,2)
                 e_tilde(:,1,2) = e_tilde(:,2,3)*e_tilde(:,3,1)&
                 &- e_tilde(:,2,1)*e_tilde(:,3,3)
                 e_tilde(:,1,3) = e_tilde(:,2,1)*e_tilde(:,3,2)&
                 &- e_tilde(:,2,2)*e_tilde(:,3,1)
                 ! trouver v_r, v_theta, v_z dans la nouvelle base
                 !v_z : vec(3). J.v=J v_z
                 Jnorm=sqrt(J100(1)**2+J100(2)**2+J100(3)**2)
                 vec(:,3) = (J100(1)*var(:,ind,2))/(Jnorm)+&
                 &(J100(2)*var(:,ind,3))/(Jnorm)+&
                 &(J100(3)*var(:,ind,4))/(Jnorm)
                 !v_theta : vec(2).
                 vec(:,2)=var(:,ind,2)*e_tilde(:,2,1)&
                 &+var(:,ind,3)*e_tilde(:,2,2)&
                 &+var(:,ind,4)*e_tilde(:,2,3)
                 !v_r : vec(1)
                 vec(:,1)=var(:,ind,2)*e_tilde(:,1,1)&
                 &+var(:,ind,3)*e_tilde(:,1,2)&
                 &+var(:,ind,4)*e_tilde(:,1,3)

                 !trouver b_theta et b_r b_phi
                 !b_z = b.e_ztilde
!                 vec(:,6) = 0.5*((var(:,ind,5)+var(:,ind,8))*e_tilde(:,3,1) &
!                 & + (var(:,ind,6)+var(:,ind,9))*e_tilde(:,3,2) &
!                 & + (var(:,ind,7)+var(:,ind,10))*e_tilde(:,3,3)) &
!                 & * sqrt(4.*acos(-1.)*ud*(ul/ut)**2)
!                 print*, vec(:,6)

                 !b_r = b.e_xtilde
!                 vec(:,7) = 0.5*((var(:,ind,5)+var(:,ind,8))*e_tilde(:,1,1) &
!                 & + (var(:,ind,6)+var(:,ind,9))*e_tilde(:,1,2) &
!                 & + (var(:,ind,7)+var(:,ind,10))*e_tilde(:,1,3)) &
!                 & * sqrt(4.*acos(-1.)*ud*(ul/ut)**2)

                 !b_theta = b.e_ytilde
!                 vec(:,8) = 0.5*((var(:,ind,5)+var(:,ind,8))*e_tilde(:,2,1) &
!                 & + (var(:,ind,6)+var(:,ind,9))*e_tilde(:,2,2) &
!                 & + (var(:,ind,7)+var(:,ind,10))*e_tilde(:,2,3)) &
!                 & * sqrt(4.*acos(-1.)*ud*(ul/ut)**2)

!J100(1)=1.
!J100(2)=0.
!J100(3)=0.
                 !trouver position dans les nouvelles coordonnees
                 z_tilde(3)=1.d0/sqrt(1.d0+(J100(3)/J100(1))**2)
                 z_tilde(2)=0.d0
                 z_tilde(1)=sqrt(1.d0-z_tilde(3)**2)
                 y_tilde(1)=z_tilde(2)*J100(3)-z_tilde(3)*J100(2)
                 y_tilde(2)=z_tilde(3)*J100(1)-z_tilde(1)*J100(3)
                 y_tilde(3)=z_tilde(1)*J100(2)-z_tilde(2)*J100(1)
                 x_tilde(:,1)=abs((x(:,1))*J100(1)+(x(:,2))*J100(2)&
                 & +(x(:,3))*J100(3))
                 x_tilde(:,2)=abs((x(:,1))*y_tilde(1)+(x(:,2))*y_tilde(2)&
                 & +(x(:,3))*y_tilde(3))
                 x_tilde(:,3)=abs((x(:,1))*z_tilde(1)+(x(:,2))*z_tilde(2)&
                 & +(x(:,3))*z_tilde(3))
                 x_tilde(:,1)=((x(:,1)-center(1))*J100(1)+(x(:,2)-center(2))*J100(2)&
                 & +(x(:,3)-center(3))*J100(3))
                 x_tilde(:,2)=((x(:,1)-center(1))*y_tilde(1)+(x(:,2)-center(2))*y_tilde(2)&
                 & +(x(:,3)-center(3))*y_tilde(3))
                 x_tilde(:,3)=((x(:,1)-center(1))*z_tilde(1)+(x(:,2)-center(2))*z_tilde(2)&
                 & +(x(:,3)-center(3))*z_tilde(3))


                 rho = var(:,ind,1)*ud*dx*boxlen*ul !rho en g/cc
              case (2)
                 rho = var(:,ind,2)*ul/ut !u en cm/s
              case (3)
                 rho = var(:,ind,3)*ul/ut ! v
              case (4)
                 rho = var(:,ind,4)*ul/ut !w
              case (5)
                 rho = var(:,ind,5) !bxl
              case (6)
                 rho = var(:,ind,6) ! byl
              case (7)
                 rho = var(:,ind,7) ! bzl
              case (8)
                 rho = var(:,ind,8) ! bxr
              case (9)
                 rho = var(:,ind,9) ! byr
              case (10)
                 rho = var(:,ind,10) ! bzr
              case (11)
                 rho = var(:,ind,11) ! P 
              case (12)
                 rho = var(:,ind,12) ! Er
              case (13)
                 rho = var(:,ind,13) ! Ei
              case (14)
                 rho = 0.5*var(:,ind,1)*(var(:,ind,2)**2+var(:,ind,4)**2+var(:,ind,3)**2) !rho*u^2
              case (15)
                 rho = sqrt(0.25*((var(:,ind,5)+var(:,ind,8))**2 & 
                      &    +(var(:,ind,6)+var(:,ind,9))**2 &
                      !&    +(var(:,ind,7)+var(:,ind,10))**2)*(4.*acos(-1.)*ud*(ul/ut)**2) &  ! B^2 en gauss
                      &    +(var(:,ind,7)+var(:,ind,10))**2)*(ud*(ul/ut)**2) &  ! B^2 en gauss
                      &    /(var(:,ind,1)*ud))/1.d5 !    
              case(16)
                 rho = sqrt(0.25*((var(:,ind,5)+var(:,ind,8))**2 & 
                     &         +(var(:,ind,6)+var(:,ind,9))**2 &
                     &         +(var(:,ind,7)+var(:,ind,10))**2) &
                     &         *4.*acos(-1.)*ud*(ul/ut)**2)            !  abs(B) en Gauss
                 !rho = (-var(:,ind,5)+var(:,ind,8))+&
                 !&(-var(:,ind,6)+var(:,ind,9))+&
                 !&(-var(:,ind,7)+var(:,ind,10))
                 !rho = sqrt(((var(:,ind,8))**2 & 
                     !&         +(var(:,ind,9))**2 &
                     !&         +(var(:,ind,10))**2) &
                     !&         *4.*acos(-1.)*ud*(ul/ut)**2)            !  abs(B) en Gauss
                 vec(:,1) = 0.5*(var(:,ind,5)+var(:,ind,8))*sqrt(4.*acos(-1.)*ud*(ul/ut)**2)  !bx
                 vec(:,2) = 0.5*(var(:,ind,6)+var(:,ind,9))*sqrt(4.*acos(-1.)*ud*(ul/ut)**2)  !by
                 vec(:,3) = 0.5*(var(:,ind,7)+var(:,ind,10))*sqrt(4.*acos(-1.)*ud*(ul/ut)**2) !bz
              case(17)
                 rho = var(:,ind,1)*ud !rho en g/cc
                 vec(:,1) = var(:,ind,2)*ul/ut      !u
                 vec(:,2) = var(:,ind,3)*ul/ut      !v
                 vec(:,3) = var(:,ind,4)*ul/ut      !w
              case(18)
!                 rho=var(:,ind,13)*utemp*0.66667 ! tp
                 rho=var(:,ind,11)/var(:,ind,1)*utemp
                 vec(:,1) = var(:,ind,1)*ud      
              case (19)
                 rho = ((var(:,ind,2)**2+var(:,ind,3)**2+var(:,ind,4)**2)**0.5)*ul/ut !u_r
                 do igrida=1,ngrida
                    if(var(igrida,ind,4).lt.0.0) rho(igrida)=-rho(igrida)
                 end do
              end select
              ! Store data cube
              do i=1,ngrida
                 ok_cell= .not.ref(i)
                 if(ok_cell)then
                    ix=int(x(i,1)*dble(nx_full))+1
                    iy=int(x(i,2)*dble(ny_full))+1
                    iz=int(x(i,3)*dble(nz_full))+1
!                    ixx=int((x_tilde(i,1)+center(1))*dble(nx_full))+1
!                    iyy=int((x_tilde(i,2)+center(2))*dble(ny_full))+1
!                    izz=int((x_tilde(i,3)+center(3))*dble(nz_full))+1
                    !if(    ix>=grid(ilevel)%imin.and.&
if (direction ==1) then
!   print*, ix
                    if(    iy>=grid(ilevel)%jmin.and.&
                         & iz>=grid(ilevel)%kmin.and.&
                         & ix>=grid(ilevel)%imin.and.&
                         & ix<=grid(ilevel)%imax.and.&
!                         & ix<=grid(ilevel)%imax.and.&
!                         & ix>=grid(ilevel)%imin.and.&
                         & iy<=grid(ilevel)%jmax.and.&
                         & iz<=grid(ilevel)%kmax)then
!                         if (x(i,2) .ne. x_tilde(i,2)) then
!                         print*, x(i,1),x(i,2),x(i,3)
!                         print*, x_tilde(i,1),x_tilde(i,2),x_tilde(i,3)
!                         read(*,*)
!                      end if
!   print*, ix
                       
                       select case (type)
                       case(-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19)
                          grid(ilevel)%cube(1,iy,iz)=rho(i)+grid(ilevel)%cube(1,iy,iz)
                       case(16,17)   
                          grid(ilevel)%cube(ix,iy,iz)=rho(i)
                          grid(ilevel)%cubevec(ix,iy,iz,1)=vec(i,1)
                          grid(ilevel)%cubevec(ix,iy,iz,2)=vec(i,2)
                          grid(ilevel)%cubevec(ix,iy,iz,3)=vec(i,3)
                       end select
                    endif
else if (direction ==2) then
                    if(    ix>=grid(ilevel)%imin.and.&
                         & iz>=grid(ilevel)%kmin.and.&
                         & iy>=grid(ilevel)%jmin.and.&
                         & iy<=grid(ilevel)%jmax.and.&
                         & ix<=grid(ilevel)%imax.and.&
                         & iz<=grid(ilevel)%kmax)then
                       select case (type)
                       case(-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19)
                          grid(ilevel)%cube(1,ix,iz)=rho(i)+grid(ilevel)%cube(1,ix,iz)
                       case(16,17)   
                          grid(ilevel)%cube(ix,iy,iz)=rho(i)
                          grid(ilevel)%cubevec(ix,iy,iz,1)=vec(i,1)
                          grid(ilevel)%cubevec(ix,iy,iz,2)=vec(i,2)
                          grid(ilevel)%cubevec(ix,iy,iz,3)=vec(i,3)
                       end select
                    endif
else if (direction ==3) then
                    if(    iy>=grid(ilevel)%jmin.and.&
                         & ix>=grid(ilevel)%imin.and.&
                         & iz>=grid(ilevel)%kmin.and.&
                         & iz<=grid(ilevel)%kmax.and.&
                         !& ix<=grid(ilevel)%imax.and.&
                         & iy<=grid(ilevel)%jmax.and.&
                         & ix<=grid(ilevel)%imax)then
                       
                       select case (type)
                       case(-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19)
                          grid(ilevel)%cube(1,iy,ix)=rho(i)+grid(ilevel)%cube(1,iy,ix)
                       case(16,17)   
                          grid(ilevel)%cube(ix,iy,iz)=rho(i)
                          grid(ilevel)%cubevec(ix,iy,iz,1)=vec(i,1)
                          grid(ilevel)%cubevec(ix,iy,iz,2)=vec(i,2)
                          grid(ilevel)%cubevec(ix,iy,iz,3)=vec(i,3)
                       end select
                    endif
end if
                 end if
              end do

           end do
           ! End loop over cell

           deallocate(xg,son,var,ref,rho,x,vec,rvectJ,e_tilde,x_tilde)
        endif

     end do
     ! End loop over levels

     close(10)
     close(11)

  end do
  ! End loop over cpus

if  ( direction == 1 ) then     ! axe x
    ! rien à faire
else if ( direction == 2 ) then   ! axe y : xx sera y, yy sera z, zz sera x
    tempmin = xmin
    tempmax = xmax
    xmin = ymin
    xmax = ymax
    ymin = tempmin
    ymax = tempmax
!    ymin = zmin
!    ymax = zmax
!    zmin = tempmin
!    zmax = tempmax
else if ( direction == 3 ) then  ! axe z : permutation circulaire 
    tempmin = xmin
    tempmax = xmax
    xmin = zmin
    xmax = zmax
    zmin = tempmin 
    zmax = tempmax
!    ymin = tempmin
!    ymax = tempmax
end if

  nx_full=2**lmax
  ny_full=2**lmax
  nz_full=2**lmax
  imin=int(xmin*dble(nx_full))+1
  imax=int(xmax*dble(nx_full))
  jmin=int(ymin*dble(ny_full))+1
  jmax=int(ymax*dble(ny_full))
  kmin=int(zmin*dble(nz_full))+1
  kmax=int(zmax*dble(nz_full))

  write(*,'(" Zoom in",3(1x,I5,"-",I5))')imin,imax,jmin,jmax,kmin,kmax
  write(*,*)'Please wait...'
  
  do ix=imin,imin!imax
     xmin=((ix-0.5)/2**lmax)
     do iy=jmin,jmax
        ymin=((iy-0.5)/2**lmax)
        do iz=kmin,kmax
           zmin=((iz-0.5)/2**lmax)
           do ilevel=max(levelmin,lmin),lmax-1
              ndom=2**ilevel
              i=int(xmin*ndom)+1
              j=int(ymin*ndom)+1
              k=int(zmin*ndom)+1
              select case (type)
              case(-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19)
!                 grid(lmax)%cube(1,ix,iz)=grid(lmax)%cube(1,ix,iz)+ &
!                      & grid(ilevel)%cube(1,i,k)
                 grid(lmax)%cube(1,iy,iz)=grid(lmax)%cube(1,iy,iz)+ &
                      & grid(ilevel)%cube(1,j,k)
              case(16,17)
                 grid(lmax)%cube(ix,iy,iz)=max(grid(lmax)%cube(ix,iy,iz), &
                      & grid(ilevel)%cube(i,j,k))
                 grid(lmax)%cubevec(ix,iy,iz,1)=grid(lmax)%cubevec(ix,iy,iz,1)+ &
                      & grid(ilevel)%cubevec(i,j,k,1)!*0.5**(-ilevel+lmax)
                 grid(lmax)%cubevec(ix,iy,iz,2)=grid(lmax)%cubevec(ix,iy,iz,2)+ &
                      & grid(ilevel)%cubevec(i,j,k,2)!*0.5**(-ilevel+lmax)
                 grid(lmax)%cubevec(ix,iy,iz,3)=grid(lmax)%cubevec(ix,iy,iz,3)+ &
                      & grid(ilevel)%cubevec(i,j,k,3)!*0.5**(-ilevel+lmax)
              end select
           end do
        end do
     end do
  end do
  select case (type)
  case(-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19)
     
     write(*,*)'Min value:',minval(grid(lmax)%cube(1,jmin:jmax,kmin:kmax))
     write(*,*)'Max value:',maxval(grid(lmax)%cube(1,jmin:jmax,kmin:kmax))
  case(16)
     write(*,*)'Min value:',minval(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))
     write(*,*)'Max value:',maxval(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))
     write(*,*)'Norm:     ',sum   (dble(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))) &
          & /(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1)
     write(*,*)'Min value:',minval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,1)) &
          &, minval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,2)) &
          &, minval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,3))
     write(*,*)'Max value:',maxval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,1)) &
          &, maxval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,2)) &
          &, maxval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,3)) 
     write(*,*)'Norm:     ',sum   (dble(abs(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,:)))) &
          & /(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1)
  !norm = 25000. !max(25000.,0.05*sum(dble(abs(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,:)))) &
          !& /(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1))

  case(17)
     write(*,*)'Min value:',minval(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))
     write(*,*)'Max value:',maxval(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))

     write(*,*)'Norm:     ',sum   (dble(grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax))) &
          & /(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1)
     write(*,*)'Min value:',minval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,1)) &
          &, minval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,2)) &
          &, minval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,3))
     write(*,*)'Max value:',maxval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,1)) &
          &, maxval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,2)) &
          &, maxval(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,3)) 
     write(*,*)'Norm:     ',sum   (dble(abs(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,:)))) &
          & /(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1)
  !norm = 25000. !max(25000.,0.05*sum(dble(abs(grid(lmax)%cubevec(imin:imax,jmin:jmax,kmin:kmax,:)))) &
          !& /(imax-imin+1)/(jmax-jmin+1)/(kmax-kmin+1))

  end select

  ! Output file
  if(TRIM(filetype).eq.'bin')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='unformatted')
     if(nx_sample==0)then
        write(20)imax-imin+1,jmax-jmin+1,kmax-kmin+1
        allocate(toto(imax-imin+1,jmax-jmin+1,kmax-kmin+1))
        toto=grid(lmax)%cube(imin:imax,jmin:jmax,kmin:kmax)
        write(*,*),'nx',imax-imin+1,jmax-jmin+1,kmax-kmin+1
        write(20)toto
     else
        if(ny_sample==0)ny_sample=nx_sample
        if(nz_sample==0)nz_sample=ny_sample
        write(20)nx_sample,ny_sample,nz_sample
        write(*,*),'nx',nx_sample,ny_sample,nz_sample
       allocate(toto(1:nx_sample,1:ny_sample,1:nz_sample))
        do i=1,nx_sample
           xx=(dble(i)-0.5)/dble(nx_sample)*dble(imax-imin+1)
           ix=int(xx)
           ddx=xx-ix
           dex=1d0-ddx
           ix=ix+imin
           ix=min(ix,imax)
           ixp1=min(ix+1,imax)
           do j=1,ny_sample
              yy=(dble(j)-0.5)/dble(ny_sample)*dble(jmax-jmin+1)
              iy=int(yy)
              ddy=yy-iy
              dey=1d0-ddy
              iy=iy+jmin
              iy=min(iy,jmax)
              iyp1=min(iy+1,jmax)
              do k=1,nz_sample
                 zz=(dble(k)-0.5)/dble(nz_sample)*dble(kmax-kmin+1)
                 iz=int(zz)
                 ddz=zz-iz
                 dez=1d0-ddz
                 iz=iz+kmin
                 iz=min(iz,kmax)
                 izp1=min(iz+1,kmax)
                 toto(i,j,k)=dex*dey*dez*log10(grid(lmax)%cube(ix  ,iy  ,iz  ))+ &
                      &      ddx*dey*dez*log10(grid(lmax)%cube(ixp1,iy  ,iz  ))+ &
                      &      dex*ddy*dez*log10(grid(lmax)%cube(ix  ,iyp1,iz  ))+ &
                      &      dex*dey*ddz*log10(grid(lmax)%cube(ix  ,iy  ,izp1))+ &
                      &      ddx*ddy*dez*log10(grid(lmax)%cube(ixp1,iyp1,iz  ))+ &
                      &      dex*ddy*ddz*log10(grid(lmax)%cube(ix  ,iyp1,izp1))+ &
                      &      ddx*dey*ddz*log10(grid(lmax)%cube(ixp1,iy  ,izp1))+ &
                      &      ddx*ddy*ddz*log10(grid(lmax)%cube(ixp1,iyp1,izp1))
              end do
           end do
        end do
        toto=10.**(toto)
        write(20)toto
     endif
     close(20)
  endif
  
  if(TRIM(filetype).eq.'grafic')then
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='unformatted')
     dummy=0.0
     write(20)imax-imin+1,jmax-jmin+1,kmax-kmin+1,dummy,dummy,dummy,dummy,dummy,dummy,dummy,dummy
     do iz=kmin,kmax
        write(20)((real(grid(lmax)%cube(ix,iy,iz)),ix=imin,imax),iy=jmin,jmax)
     end do
     close(20)
  endif
  
  if(TRIM(filetype).eq.'vtk')then       ! Bon en fait c'est du gnuplot maintenant mais bon...
!      nomfich=TRIM(outfich)//'_'//TRIM(str(mod(direction+corr,3)))
     nomfich=TRIM(outfich)
     write(*,*)'Writing file '//TRIM(nomfich)
     open(unit=20,file=nomfich,form='formatted')
!      open(unit=30,file='header.gnu',form='formatted')
!      write(30,'("DATASET STRUCTURED_POINTS")')
!     write(30,'("DIMENSIONS ",3(I4,1x))')imax-imin+1,jmax-jmin+1,kmax-kmin+1
!     write(30,'("ORIGIN 0.0 0.0 0.0")')
!     write(30,'("SPACINGS 1.0 1.0 1.0")')
!      write(30,*)'BOXLEN :', boxlen
!      write(30,*)'TIME :', t*scale_t_kyrs, ' Kyears'
!      write(30,'(" ")')
    

norm = 10.*norm/((jmax-jmin)*0.5**lmax*scale_ua*boxlen)
     if  (old == 0) then    ! cube pour gnuplot
     select case(type)
     case(-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19)
!         write(30,'("POINT_DATA ",I10)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
!         write(30,'("SCALARS values float")')
!         write(30,'("LOOKUP_TABLE default")')
        write(20,'(a,es16.6e3)')"t =",t*scale_t_kyrs
              write(20,*)
              write(20,*)
        do k=kmin,kmax
           zz=dble(k-1)/(2**lmax)-0.5          
!           do l=1,2
              do j=jmin,jmax
                 yy=dble(j-1)/(2**lmax)-0.5         
                 do i=imin,imin!imax
                        xx=dble(i-1)/(2**lmax)-0.5           
                        !write(20,'(E9.4)')xx,yy,zz,real((grid(lmax)%cube(i,j,k)))
!                            write(20,*)yy,xx*boxlen*scale_ua,zz*boxlen*scale_ua,&
!                                   &real((grid(lmax)%cube(1,i,k)))
                            write(20,'(10(es16.6e3,1x))')xx,yy*boxlen*scale_ua,zz*boxlen*scale_ua,&
                                   &real((grid(lmax)%cube(1,j,k)))
                 end do
              end do
              write(20,*)
!           end do
        end do
     
     case(17)
!         write(30,'("POINT_DATA ",I10)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
!         write(30,'("SCALARS values float")')
!         write(30,'("LOOKUP_TABLE default")')
              write(20,'(a,es16.6e3)')"t =",t*scale_t_kyrs
              write(20,*)
              write(20,*)
        do k=kmin,kmax
           zz=dble(k-1)/(2**lmax)-0.5          
!           do l=1,2
              do j=jmin,jmax
                 yy=dble(j-1)/(2**lmax)-0.5         
                 do i=imin,imax
                        xx=dble(i-1)/(2**lmax)-0.5           
                        !write(20,'(E9.4)')xx,yy,zz,real((grid(lmax)%cube(i,j,k)))
                        if  (direction==1) then
                            write(20,'(10(es16.6e3,1x))')xx,yy*boxlen*scale_ua,zz*boxlen*scale_ua,&
                                   &real((grid(lmax)%cube(i,j,k)))
                        elseif (direction==2) then
                            write(20,'(10(es16.6e3,1x))')xx,yy*boxlen*scale_ua,zz*boxlen*scale_ua,&
                                   &real((grid(lmax)%cube(k,i,j)))
                        elseif (direction==3) then
                            write(20,'(10(es16.6e3,1x))')xx,yy*boxlen*scale_ua,zz*boxlen*scale_ua,&
                                   &real((grid(lmax)%cube(j,k,i)))
                        endif
                 end do
              end do
              write(20,*)
!           end do
        end do
              write(20,*)
              write(20,*)

     close(20)
        open(unit=21,file='cubev.dat',form='formatted')
!         write(30,*)'VECTORS :'
!         write(30,'("POINT_DATA ",I8)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
!         write(30,'("VECTORS values float")')
        !write(21,*) 0.0, (dble(jmin+10)/(2**lmax)-0.5)*boxlen*scale_ua,(dble(jmin-5)/(2**lmax)-0.5)*boxlen*scale_ua, 0., 20., 0.
        !write(21,*) 
        !write(21,*) 
!        write(20,'("LOOKUP_TABLE default")')
        do k=kmin, kmax,20
           zz=dble(k-1)/(2**lmax)-0.5
           do j=jmin,jmax,20
              yy=dble(j-1)/(2**lmax)-0.5          
              do i=imin,imax,20
                 xx=dble(i-1)/(2**lmax)-0.5          
                 !write(20,'(3F9.3)')xx,yy,zz, xx+real((grid(lmax)%cubevec(i,j,k,1)))/norm& 
                 !write(20,*)xx+1./(2**(lmax+1)),yy+1./(2**(lmax+1)),zz+1./(2**(lmax+1)) &
                if  (direction==1) then
                    !print*,abs(yy+real((grid(lmax)%cubevec(i,j,k,2)))/norm/boxlen/scale_ua)*2**lmax,&
                    !&(grid(ilevel)%jmax-grid(ilevel)%jmin)/2.
                 if  ((abs(yy+0.5-center(2)+real((grid(lmax)%cubevec(i,j,k,2)))/norm/boxlen/scale_ua)*2**lmax &
                      & .le. (grid(ilevel)%jmax-grid(ilevel)%jmin)/2.) .and. &
                      &(abs(zz+0.5-center(3)+real((grid(lmax)%cubevec(i,j,k,3)))/norm/boxlen/scale_ua)*2**lmax &
                      & .le. (grid(ilevel)%jmax-grid(ilevel)%jmin)/2.) ) then
                     write(21,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua &
                      & ,real((grid(lmax)%cubevec(i,j,k,1)))/norm& 
                      & ,real((grid(lmax)%cubevec(i,j,k,2)))/norm& 
                      & ,real((grid(lmax)%cubevec(i,j,k,3)))/norm 
                 endif
                elseif (direction==2) then
!                print*,grid(ilevel)%imax,grid(ilevel)%imin
!                print*,grid(ilevel)%jmax,grid(ilevel)%jmin
!                print*,grid(ilevel)%kmax,grid(ilevel)%kmin
!                stop
                 if  ((abs(yy+0.5-center(3)+real((grid(lmax)%cubevec(k,i,j,3)))&
                      &/norm/boxlen/scale_ua)*2**lmax &
                      & .le. (grid(ilevel)%kmax-grid(ilevel)%kmin)/2.) .and. &
                      &(abs(zz+0.5-center(1)+real((grid(lmax)%cubevec(k,i,j,1)))&
                      &/norm/boxlen/scale_ua)*2**lmax &
                      & .le. (grid(ilevel)%imax-grid(ilevel)%imin)/2.) ) then
                     write(21,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua &
                      & ,real((grid(lmax)%cubevec(k,i,j,2)))/norm& 
                      & ,real((grid(lmax)%cubevec(k,i,j,3)))/norm& 
                      & ,real((grid(lmax)%cubevec(k,i,j,1)))/norm 
                 endif
                elseif (direction==3) then
                 if  ((abs(yy+0.5-center(1)+real((grid(lmax)%cubevec(j,k,i,1)))/norm/boxlen/scale_ua)*2**lmax &
                      & .le. (grid(ilevel)%imax-grid(ilevel)%imin)/2.) .and. &
                      &(abs(zz+0.5-center(2)+real((grid(lmax)%cubevec(j,k,i,2)))/norm/boxlen/scale_ua)*2**lmax &
                      & .le. (grid(ilevel)%jmax-grid(ilevel)%jmin)/2.) ) then
                     write(21,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua &
                      & ,real((grid(lmax)%cubevec(j,k,i,3)))/norm& 
                      & ,real((grid(lmax)%cubevec(j,k,i,1)))/norm& 
                      & ,real((grid(lmax)%cubevec(j,k,i,2)))/norm 
                 endif
                endif
              end do
           end do
        end do
        write(21,*) 
        write(21,*) 
        close(21)        

     case(16)
!         write(30,'("POINT_DATA ",I10)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
!         write(30,'("SCALARS values float")')
!         write(30,'("LOOKUP_TABLE default")')
              write(20,'(a,es16.6e3)')"t =",t*scale_t_kyrs
              write(20,*)
              write(20,*)
        do k=kmin,kmax
           zz=dble(k-1)/(2**lmax)-0.5          
!           do l=1,2
              do j=jmin,jmax
                 yy=dble(j-1)/(2**lmax)-0.5         
                 do i=imin,imax
                        xx=dble(i-1)/(2**lmax)-0.5           
                        !write(20,'(E9.4)')xx,yy,zz,real((grid(lmax)%cube(i,j,k)))
                        if  (direction==1) then
                            write(20,'(10(es16.6e3,1x))')xx,yy*boxlen*scale_ua,zz*boxlen*scale_ua,&
                                   &real((grid(lmax)%cube(i,j,k)))
                        elseif (direction==2) then
                            write(20,'(10(es16.6e3,1x))')xx,yy*boxlen*scale_ua,zz*boxlen*scale_ua,&
                                   &real((grid(lmax)%cube(k,i,j)))
                        elseif (direction==3) then
                            write(20,'(10(es16.6e3,1x))')xx,yy*boxlen*scale_ua,zz*boxlen*scale_ua,&
                                   &real((grid(lmax)%cube(j,k,i)))
                        endif
                 end do
              end do
              write(20,*)
!           end do
        end do
     close(20)
        open(unit=22,file='cubeb.dat',form='formatted')
!         write(30,*)'VECTORS :'
!         write(30,'("POINT_DATA ",I8)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
!         write(30,'("VECTORS values float")')
        !write(22,*) 0.0, (dble(jmin+10)/(2**lmax)-0.5)*boxlen*scale_ua,(dble(jmin-5)/(2**lmax)-0.5)*boxlen*scale_ua, 0., 20., 0.
        !write(22,*) 
        !write(22,*) 
!        write(20,'("LOOKUP_TABLE default")')
        do k=kmin+10,kmax-5,12
           zz=dble(k-1)/(2**lmax)-0.5
           do j=jmin+10,jmax-5,12
              yy=dble(j-1)/(2**lmax)-0.5          
              do i=imin,imax,12
                 xx=dble(i-1)/(2**lmax)-0.5          
                 !write(20,'(3F9.3)')xx,yy,zz, xx+real((grid(lmax)%cubevec(i,j,k,1)))/norm& 
                 !write(20,*)xx+1./(2**(lmax+1)),yy+1./(2**(lmax+1)),zz+1./(2**(lmax+1)) &
                 if  (direction==1) then
                     write(22,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua &
                      & ,real((grid(lmax)%cubevec(i,j,k,1)))/norm& 
                      & ,real((grid(lmax)%cubevec(i,j,k,2)))/norm& 
                      & ,real((grid(lmax)%cubevec(i,j,k,3)))/norm 
                elseif (direction==2) then
                     write(22,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua &
                      & ,real((grid(lmax)%cubevec(k,i,j,2)))/norm& 
                      & ,real((grid(lmax)%cubevec(k,i,j,3)))/norm& 
                      & ,real((grid(lmax)%cubevec(k,i,j,1)))/norm 
                elseif (direction==3) then
                     write(22,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua &
                      & ,real((grid(lmax)%cubevec(j,k,i,3)))/norm& 
                      & ,real((grid(lmax)%cubevec(j,k,i,1)))/norm& 
                      & ,real((grid(lmax)%cubevec(j,k,i,2)))/norm 
                endif
              end do
           end do
        end do
        close(22)        

     end select
     end if
     
     if  (old == 1) then     ! l'original, le cube normal
     select case(type)
     case(-1,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,18,19)
!         write(30,'("POINT_DATA ",I10)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
!         write(30,'("SCALARS values float")')
!         write(30,'("LOOKUP_TABLE default")')
        do k=kmin,kmax
           zz=dble(k-1)/(2**lmax)-0.5           
           do j=jmin,jmax
              yy=dble(j-1)/(2**lmax)-0.5          
              do i=imin,imax
                     xx=dble(i-1)/(2**lmax)-0.5          
                     !write(20,'(E9.4)')xx,yy,zz,real((grid(lmax)%cube(i,j,k)))
                     if  ((direction == 1).or.(direction==0)) write(20,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua, &
                                                     & real((grid(lmax)%cube(i,j,k))),real((grid(lmax)%cubevec(i,j,k,1)))
                     if  (direction == 2) write(20,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua, &
                                                     & real((grid(lmax)%cube(k,i,j)))
                     if  (direction == 3) write(20,'(10(es16.6e3,1x))')xx*boxlen*scale_ua,yy*boxlen*scale_ua,zz*boxlen*scale_ua, & 
                                                     & real((grid(lmax)%cube(j,k,i)))
              end do
           end do
           write(20,*)
        end do
     case(16,17)
!         write(30,'("POINT_DATA ",I8)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)
!         write(30,'("VECTORS values float")')
!        write(20,'("LOOKUP_TABLE default")')
        do k=kmin+1,kmax
           zz=dble(k-1)/(2**lmax)-0.5           
           do j=jmin+1,jmax
              yy=dble(j-1)/(2**lmax)-0.5          
              do i=imin,imax
                 xx=dble(i-1)/(2**lmax)-0.5          
                 !write(20,'(3F9.3)')xx,yy,zz, xx+real((grid(lmax)%cubevec(i,j,k,1)))/norm& 
                 !write(20,*)xx+1./(2**(lmax+1)),yy+1./(2**(lmax+1)),zz+1./(2**(lmax+1)) &
                 write(20,'(10(es16.6e3,1x))')xx,yy,zz &
                      & ,real((grid(lmax)%cubevec(i,j,k,1)))/norm& 
                      & ,real((grid(lmax)%cubevec(i,j,k,2)))/norm& 
                      & ,real((grid(lmax)%cubevec(i,j,k,3)))/norm 
              end do
           end do
        end do
        
     end select
     end if

!      close(30)
     close(20)
  endif
  if(TRIM(filetype).eq.'radmc')then
     select case(type)
     case(1)
        ! Write coordinates in cgs
        nomfich='amr_grid.inp'
        write(*,*)'Writing file '//TRIM(nomfich)
        open(unit=20,file=nomfich,form='formatted')
        write(20,'("1")') ! Default
        write(20,'("0")') ! Amrstyle : uniform grid
        write(20,'(I3)') 10 ! Cartesian coordinates system
        write(20,'("0")') ! Grid info
        write(20,'("1  1  1")') ! 3D coordinates
        write(20,'(3(I3,1x))')imax-imin+1,jmax-jmin+1,kmax-kmin+1 ! Number of grid cells
        
        dx  = 1.0d0/2.0d0**(lmax)!(xxmax-xxmin)*(2**lmax)
        eps = 0.5 - xxmin 
        print*,eps,dx,ul,boxlen,xxmin
        do i=1,imax-imin+2
           pos = (-eps + dx * (i-1)) * ul * boxlen
           write(20,'(E13.7)')pos
        end do
        dy = (yymax-yymin)*(2**lmax)
        eps = 0.5 - yymin 
        do i=1,jmax-jmin+2
           pos = (-eps + dx * (i-1)) * ul * boxlen           
           write(20,'(E13.7)')pos
        end do
        dz = (zzmax-zzmin)*(2**lmax)
        eps = 0.5 - zzmin 
        do i=1,kmax-kmin+2
!           pos = (-eps + eps*dx * (i-1)) * ul * boxlen
           pos = (-eps + dx * (i-1)) * ul * boxlen
           write(20,'(E13.7)')pos
        end do
        close(20)
        
        !Write density
        nomfich='dust_density.inp'
        write(*,*)'Writing file '//TRIM(nomfich)
        open(unit=20,file=nomfich,form='formatted')
        write(20,'("1")') ! Default
        write(20,'(I10)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1) ! Number of cells
        write(20,'("1")') ! Number of dust species
        do k=kmin,kmax
           do j=jmin,jmax
              do i=imin,imax
                 write(20,'(E9.4)')real((grid(lmax)%cube(i,j,k)))/100.0d0
              end do
           end do
        end do
maxtp=0.0
     case(18)        
        nomfich='dust_temperature.inp'
        write(*,*)'Writing file '//TRIM(nomfich)
        open(unit=20,file=nomfich,form='formatted')
        write(20,'("1")') ! Default
        write(20,'(I10)')(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1) ! Number of cells
        write(20,'("1")') ! Number of dust species
        do k=kmin,kmax
           do j=jmin,jmax
              do i=imin,imax
                 write(20,'(E9.4)')real((grid(lmax)%cube(i,j,k)))
                 if(real((grid(lmax)%cube(i,j,k))) .gt. maxtp) maxtp=real((grid(lmax)%cube(i,j,k)))
              end do
           end do
        end do
     
     end select

     close(20)
print*,maxtp
  endif
  
contains
  
  subroutine read_params
    
    implicit none
    
    integer       :: i,n
    integer       :: iargc
    character(len=4)   :: opt
    character(len=128) :: arg
    LOGICAL       :: bad, ok
    real          :: fenetre
    
    n = iargc()
    if (n < 4) then
       print *, 'usage: amr2cube  -inp  input_dir'
       print *, '                 -out  output_file'
       print *, '                 [-xmi xmin] '
       print *, '                 [-xma xmax] '
       print *, '                 [-ymi ymin] '
       print *, '                 [-yma ymax] '
       print *, '                 [-zmi zmin] '
       print *, '                 [-zma zmax] '
       print *, '                 [-lma lmax] '
       print *, '                 [-typ type] '
       print *, '                 [-fil filetype] '
       print *, 'ex: amr2cube -inp output_00001 -out cube.dat'// &
            &   ' -typ 1 -xmi 0.1 -xma 0.7 -lma 12'
       print *, ' '
       print *, ' type :-1 = cpu number'
       print *, ' type : 0 = ref. level (default)'
       print *, ' type : 1-9 = variable number'
       stop
    end if
    
    do i = 1,n,2
       call getarg(i,opt)
       if (i == n) then
          print '("option ",a2," has no argument")', opt
          stop 2
       end if
       call getarg(i+1,arg)
       select case (opt)
       case ('-inp')
          repository = trim(arg)
       case ('-out')
          outfich = trim(arg)
       case ('-fil')
          filetype = trim(arg)
       case ('-xmi')
          read (arg,*) xmin
       case ('-xma')
          read (arg,*) xmax
       case ('-ymi')
          read (arg,*) ymin
       case ('-yma')
          read (arg,*) ymax
       case ('-zmi')
          read (arg,*) zmin
       case ('-zma')
          read (arg,*) zmax
       case ('-lma')
          read (arg,*) lmax
       case ('-nx')
          read (arg,*) nx_sample
       case ('-ny')
          read (arg,*) ny_sample
       case ('-nz')
          read (arg,*) nz_sample
       case ('-typ')
          read (arg,*) type
       case ('-fen')
          read (arg,*) fenetre
       case ('-nor')
          read (arg,*) norm
       case ('-old')
          read (arg,*) old
       case ('-dir')
          read (arg,*) direction
       case default
          print '("unknown option ",a2," ignored")', opt
       end select
    end do
    
    center = (/0.5,0.5,0.5/)
    
    xmin = 0.5 - 0.5**fenetre
    xmax = 0.5 + 0.5**fenetre
    ymin = 0.5 - 0.5**fenetre
    ymax = 0.5 + 0.5**fenetre
    zmin = 0.5 - 0.5**fenetre
    zmax = 0.5 + 0.5**fenetre
    
    return
    
  end subroutine read_params

end program amr2cube

!=======================================================================
subroutine title(n,nchar)
!=======================================================================
  implicit none
  integer::n
  character*5::nchar

  character*1::nchar1
  character*2::nchar2
  character*3::nchar3
  character*4::nchar4
  character*5::nchar5

  if(n.ge.10000)then
     write(nchar5,'(i5)') n
     nchar = nchar5
  elseif(n.ge.1000)then
     write(nchar4,'(i4)') n
     nchar = '0'//nchar4
  elseif(n.ge.100)then
     write(nchar3,'(i3)') n
     nchar = '00'//nchar3
  elseif(n.ge.10)then
     write(nchar2,'(i2)') n
     nchar = '000'//nchar2
  else
     write(nchar1,'(i1)') n
     nchar = '0000'//nchar1
  endif


end subroutine title
!================================================================
!================================================================
!================================================================
!================================================================
subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer     ,INTENT(IN)                     ::bit_length,npoint
  integer     ,INTENT(IN) ,dimension(1:npoint)::x,y,z
  real(kind=8),INTENT(OUT),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask,y_bit_mask,z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i,ip,cstate,nstate,b0,b1,b2,sdigit,hdigit

  if(bit_length>bit_size(bit_length))then
     write(*,*)'Maximum bit length=',bit_size(bit_length)
     write(*,*)'stop in hilbert3d'
     stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint

     ! convert to binary
     do i=0,bit_length-1
        x_bit_mask(i)=btest(x(ip),i)
        y_bit_mask(i)=btest(y(ip),i)
        z_bit_mask(i)=btest(z(ip),i)
     enddo

     ! interleave bits
     do i=0,bit_length-1
        i_bit_mask(3*i+2)=x_bit_mask(i)
        i_bit_mask(3*i+1)=y_bit_mask(i)
        i_bit_mask(3*i  )=z_bit_mask(i)
     end do

     ! build Hilbert ordering using state diagram
     cstate=0
     do i=bit_length-1,0,-1
        b2=0 ; if(i_bit_mask(3*i+2))b2=1
        b1=0 ; if(i_bit_mask(3*i+1))b1=1
        b0=0 ; if(i_bit_mask(3*i  ))b0=1
        sdigit=b2*4+b1*2+b0
        nstate=state_diagram(sdigit,0,cstate)
        hdigit=state_diagram(sdigit,1,cstate)
        i_bit_mask(3*i+2)=btest(hdigit,2)
        i_bit_mask(3*i+1)=btest(hdigit,1)
        i_bit_mask(3*i  )=btest(hdigit,0)
        cstate=nstate
     enddo

     ! save Hilbert key as double precision real
     order(ip)=0.
     do i=0,3*bit_length-1
        b0=0 ; if(i_bit_mask(i))b0=1
        order(ip)=order(ip)+dble(b0)*dble(2)**i
     end do

  end do

end subroutine hilbert3d
!================================================================
!================================================================
!================================================================
!================================================================
function str(k)
!   "Convert an integer to string."
    character(len=20) :: str
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
end function str
