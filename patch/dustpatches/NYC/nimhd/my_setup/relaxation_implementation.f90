subroutine relaxation_implementation(ilevel,nstp)
  use amr_commons
  use hydro_commons
  use units_commons
  use poisson_commons
  implicit none

  integer::ilevel,nstp
  !--------------------------------------------------------------------------
  ! Enforcing relaxation damping to velocity components at early times,     !
  ! of in geometric regions. This now happens after set_uold, so do to uold !
  !--------------------------------------------------------------------------
  integer::i,ivar,idim,ind,ncache,igrid,iskip
  integer::ix,iy,iz,id1,ig1,ih1,id2,ig2,ih2
  integer::ngrid,nx_loc
  integer,dimension(1:nvector),save::ind_grid,ind_cell
  real(dp) ::A,B,C,e_mag,e_kin,enint,d,u,v,w
  real(dp)::dx,scale,v2,rad2,relax_f,x0,xn,yn,zn,r0,cs0,cs2,rout,rho0,H,Bz,Mstar,HoverR,H1
  real(dp):: v_quasikep,relax,rin,rr,radiusin,dx_loc
  real(dp),dimension(ndim)::vel

  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp):: epsilon_0,sum_dust,alpha_disk
  integer::idust
#if NDUST>0
  real(dp),dimension(1:ndust):: dustMRN
#endif
    integer,dimension(1:3,1:2,1:8)::iii,jjj

#if NDUST>0
  epsilon_0 = dust_ratio(1)
#endif
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel
  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
  relax=1000.0/real(nstp,dp)+1.0d0

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  x0 = boxlen/2.0d0
  
  dx=0.5d0**ilevel
  dx_loc=dx*scale
  r0=5.0d0
  rin= 0.2d0
  rout= 5.5d0
  rho0 = 2.3434e-11/scale_d
    HoverR=0.05
  H=HoverR*r0
  Mstar = 1.0d0
 bz=1.d-7
      cs0 =  (H*sqrt(Mstar/r0**3.0))**2.0
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
!!$   ! Gather neighboring grids
!!$     do i=1,ngrid
!!$        igridn(i,0)=ind_grid(i)
!!$     end do
!!$     do idim=1,ndim
!!$        do i=1,ngrid
!!$           ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
!!$           ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
!!$           igridn(i,2*idim-1)=son(ind_left (i,idim))
!!$           igridn(i,2*idim  )=son(ind_right(i,idim))
!!$        end do
!!$     end do
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
      do i=1,ngrid
         xn=(xx(i,1)-x0)
         yn=(xx(i,2)-x0)
         zn=(xx(i,3)-x0)

            RR = sqrt(xn**2.0+yn**2.0)
            radiusin=sqrt(rin**2.0+zn**2.0)
            H1= HoverR*RR
            if (nstp<100) then
            u=0.0d0; v=0.0d0; w=0.0d0
            d=max(uold(ind_cell(i),1),smallr)
            cs2= cs0*(RR/R0)**(-1.0d0)
            sum_dust=0.0d0
            do idust = 1, ndust
               sum_dust=sum_dust+uold(ind_cell(i),firstindex_ndust+idust)/d
            end do
            
            enint= d*cs2*(1.0-sum_dust)/(gamma-1.0d0)
            !uold(ind_cell(i),2)=  uold(ind_cell(i),2)/relax
            !uold(ind_cell(i),3)=  uold(ind_cell(i),3)/relax
            uold(ind_cell(i),4)= uold(ind_cell(i),4)/relax
            if(ndim>0)u=uold(ind_cell(i),2)/d
            if(ndim>1)v=uold(ind_cell(i),3)/d
            if(ndim>2)w=uold(ind_cell(i),4)/d

            e_mag= 0.0_dp
#ifdef SOLVERmhd                           
            A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
            B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
            C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
            e_mag=0.5d0*(A**2+B**2+C**2)
#endif
            e_kin=0.5d0*d*(u**2+v**2+w**2)
#if NENER>0
            do irad=1,nener
               e_mag=e_mag+uold(ind_cell(i),8+irad)
            end do
#endif
            uold(ind_cell(i),5)=e_kin+e_mag+enint
            endif

         
         if (RR<rin) then

            cs2= cs0*(rin/R0)**(-1.0d0)
            alpha_disk= Mstar/cs2

            uold(ind_cell(i),1)= 1.d-20/scale_d

#if NDUST>0
        if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
        do idust =1,ndust


           uold(ind_cell(i), firstindex_ndust+idust)= dust_ratio(idust)/(1.0d0+dust_ratio(idust))*uold(ind_cell(i),1)
           if(uold(ind_cell(i),1)<1e-17/scale_d) uold(ind_cell(i), firstindex_ndust+idust)=0.0d0
           if(mrn)  uold(ind_cell(i), firstindex_ndust+idust)= dustMRN(idust)*uold(ind_cell(i),1)
           sum_dust = sum_dust +uold(ind_cell(i), firstindex_ndust+idust)/uold(ind_cell(i),1)
!!$           v_dust(ind_cell(i),idust,1)=0.0d0
!!$           v_dust(ind_cell(i),idust,2)=0.0d0
!!$           v_dust(ind_cell(i),idust,3)=0.0d0

        end do   
#endif 
            uold(ind_cell(i),2)= 0.0d0
            uold(ind_cell(i),3)= 0.0d0
            uold(ind_cell(i),4)=0.0d0
          e_mag=0.d0
         !if(nstp>100) then
            uold(ind_cell(i),6)     = 0.d0
            uold(ind_cell(i),7)     = 0.d0
            uold(ind_cell(i),8)     = 0.0d0!1./sqrt(0.5*Bz/(rho0*cs2*(1.0-sum_dust)*(rin/r0)**(-1.5d0)))
            uold(ind_cell(i),nvar+1)= 0.d0
            uold(ind_cell(i),nvar+2)= 0.d0
            uold(ind_cell(i),nvar+3)= uold(ind_cell(i),8)
            
         !endif
#ifdef SOLVERmhd                           
            A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
            B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
            C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
            e_mag=0.5d0*(A**2+B**2+C**2)
#endif
#if NENER>0
            do irad=1,nener
               e_mag=e_mag+uold(ind_cell(i),8+irad)
            end do
#endif         
            uold(ind_cell(i),5)= uold(ind_cell(i),1)*cs2*(1.0d0-sum_dust)+e_mag
         endif
 
      enddo

   enddo
    enddo
111 format('   Entering relaxation_implementation for level ',I2)

end subroutine relaxation_implementation




subroutine turb_initialisation(ilevel)
  use amr_commons
  use hydro_commons
  use units_commons
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
  real(dp) ::A,B,C,e_mag,e_kin,enint,d,u,v,w
  real(dp)::dx,scale,v2,rad2,relax_f,x0,xn,yn,r0,cs0,cs,rout,rho0,H,Bz,Mstar
  real(dp)::relax,rin,rr
  real(dp),dimension(ndim)::v_turb

  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp):: epsilon_0,sum_dust,vrms
  integer::idust
#if NDUST>0
  real(dp),dimension(1:ndust):: dustMRN
#endif
#if NDUST>0
  epsilon_0 = dust_ratio(1)
#endif
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  vrms=0.01
  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  x0 = boxlen/2.0d0
  

  r0=5.0d0
  rin= 0.2d0
  rout= 5.5d0
  rho0 = 2.3434e-11/scale_d
  H=0.05*r0
  Bz=1.0d0
  Mstar = 1.0d0

      cs0 =  (H*sqrt(Mstar/r0**3.0))**2.0
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
      do i=1,ngrid
  
         xn=(xx(i,1)-x0)
         yn=(xx(i,2)-x0)
         
         RR = sqrt(xn**2.0+yn**2.0)

         if (RR>rin) then
            d=max(uold(ind_cell(i),1),smallr)
            cs= sqrt(cs0*(rr/R0)**(-1.0d0))
            call get_vturb(vrms,cs,v_turb)
            uold(ind_cell(i),2)=uold(ind_cell(i),2)+uold(ind_cell(i),1)*v_turb(1)
            uold(ind_cell(i),3)=uold(ind_cell(i),3)+uold(ind_cell(i),1)*v_turb(2)
            uold(ind_cell(i),4)=uold(ind_cell(i),4)+uold(ind_cell(i),1)*v_turb(3)
         endif
      enddo

   enddo
    enddo
write (*,*) 'relaxation is on'
111 format('   Entering turbulence init for level ',I2)

end subroutine turb_initialisation
subroutine get_vturb(vrms,cs,v_turb)

  use amr_commons, only:myid,ncpu
  use pm_commons, only:localseed,IRandNumSize,iseed
  use random

  implicit none

  integer :: i
  integer ,dimension(1:ncpu,1:IRandNumSize)    :: allseed
  double precision ::  vrms, cs, vel
  double precision, dimension(3) :: v_turb

  double precision :: u1, v1, u2, v2, theta, phi, x, y, z

#ifdef DEBUGRANDOM
  logical, save :: first=.true.
#endif

  if (localseed(1)==-1) then
     call rans(ncpu,iseed,allseed)
     localseed = allseed(myid,1:IRandNumSize)
  end if

  ! magnitude --> Gressel is v_rms = 0.01 * cs
  if (vrms .lt. 0.d0) then
    call gaussdev(localseed,vel)
    vel = vel * abs(vrms)*cs
  else
    call gaussdev(localseed,vel)
    vel = vel * vrms
  endif

  call ranf(localseed,v_turb(1))
  call ranf(localseed,v_turb(2))
  call ranf(localseed,v_turb(3))

  v_turb = (v_turb - 0.5d0)*2.d0
  v_turb = v_turb/sqrt(sum((v_turb(1:3))**2)) * vel

  ! NEED TO HAVE Cs OR Vrms in code units.

#ifdef DEBUGRANDOM
  if (myid == 1 .and. first) then
    open(42,file="test_gauss.dat",status='unknown')
    do i=1,10000
      call gaussdev(localseed,vel)
      write(42,*) vel
    enddo
    close(42)
    first = .false.
  endif
#endif

end subroutine get_vturb


subroutine mag_initialisation(ilevel)
  use amr_commons
  use hydro_commons
  use units_commons
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
  real(dp) ::A,B,C,e_mag,e_kin,enint,d,u,v,w
  real(dp)::dx,scale,v2,rad2,relax_f,x0,xn,yn,r0,cs0,cs2,rout,rho0,H,Bz,Mstar
  real(dp)::relax,rin,rr
  real(dp),dimension(ndim)::vel

  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  real(dp):: epsilon_0,sum_dust
  integer::idust
#if NDUST>0
  real(dp),dimension(1:ndust):: dustMRN
#endif
#if NDUST>0
  epsilon_0 = dust_ratio(1)
#endif
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  relax=1000.0

  ! Mesh spacing at that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  x0 = boxlen/2.0d0
  

  Bz=1.0d0

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
      do i=1,ngrid
 
         if (RR>rin) then
            uold(ind_cell(i),6)     = 0.d0
            uold(ind_cell(i),7)     = 0.d0
            uold(ind_cell(i),8)     = Bz

      ! Right B fields. Div * B is zero with linear operator (Bxr - Bxl)/dx ...
            uold(ind_cell(i),nvar+1)= 0.d0
            uold(ind_cell(i),nvar+2)= 0.d0
            uold(ind_cell(i),nvar+3)= uold(i,8)
         endif
      enddo

   enddo
    enddo
write (*,*) 'relaxation is on'
111 format('   Entering relaxation_implementation for level ',I2)

end subroutine mag_initialisation
