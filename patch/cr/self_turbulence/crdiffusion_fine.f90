!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crdiff_fine(ilevel,compute)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the thermal conduction scheme.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,compute,igroup
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call crdifffine1(ind_grid,ngrid,ilevel,compute,igroup)
  end do

111 format('   Entering conduction_fine for level ',i2)

end subroutine crdiff_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine crdifffine1(ind_grid,ncache,ilevel,compute,igroup)
  use amr_commons
  use hydro_commons
  use poisson_commons
!  use radiation_parameters,ONLY:dt_imp
  use radiation_parameters,ONLY:mu_gas
  use cooling_module
  implicit none
  integer::ilevel,ncache,compute,igroup
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolate from
  ! coarser level missing grid variables. It then calls the
  ! thermal conduction solver that compute energy flux. This flux is zeroed at 
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated 
  ! and stored in array unew(:), both at the current level and at the 
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uold_loc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+7),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::facdx
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim),save::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim),save::tmp

  !real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::xloc=0.0d0
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2),save::residu_loc=0.0d0

  integer,dimension(1:nvector),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist
  real(dp),dimension(1:nvector),save:: residu

  integer::neul=5
  integer::ind_buffer1,ind_buffer2,ind_buffer3
  integer::ind_father1,ind_father2,ind_father3
  integer::i,j,ivar,idim,irad,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dflux_x,dflux_y,dflux_z
  real(dp)::dx,scale,oneontwotondim,length,twodx,twodx2,bnorm2
  real(dp)::dflux,weight
  real(dp)::eth,ekin,emag,erad,bb

  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_omega,scale_tau,scale_kappa
  real(dp)::vx,vy,vz,dens,bnorm,bx,by,bz,vnorm,va,Ma,Dpara,kperp,kpar,deCRL,deCRR,deCRC,deCR
  integer::ind_cr1,l
  real(dp)::kparasub,kperpsub,Temp,gradPcr, gradPcr_x, gradPcr_y, gradPcr_z

!!$  scale_tau=1d0/(scale_t/scale_l**3) ! Time in code units
  ind_cr1=9
  if(twotemp)ind_cr1=10
  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_kappa = scale_l**2/scale_t
  kpar=Dcr/scale_kappa
  kpar=1.d10*flinj*boxlen*scale_l/scale_kappa ! Linj*c/3

  oneontwotondim = 1.d0/dble(twotondim)

  residu     = 0.0d0
  residu_loc = 0.0d0

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale
  twodx =dx*dx
  twodx2=twodx**2
 
 ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max
     
     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
           nbuffer=nbuffer+1
           ind_nexist(nbuffer)=i
           ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do

     ! If not, interpolate variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar+3
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do
           do i=1,nbuffer
              ind1(i,j)=son(ibuffer_father(i,j))
           end do
        end do
        call interpol_hydro_cond(u1,ind1,u2,nbuffer)
     end if

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do
        
        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2
        
        ! Gather hydro variables
        do ivar=1,nvar+3
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
              uold_loc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
              facdx(ind_exist(i),i3,j3,k3)=1.0d0

              if(compute==2 .and. ivar==2)then 
                 uloc(ind_exist(i),i3,j3,k3,ivar)=unew(ind_cell(i),ivar)
                 if(son(ind_cell(i))>0) uloc(ind_exist(i),i3,j3,k3,ivar)=0.0d0! neighbor cell at ilevel+1, put value to zero in vector because it is considered as a Dirichlet BC, i.e. in the RHS
              endif

           end do
           do i=1,nbuffer
!!$              if(interpol_type_cond.gt.0)then
!!$                 facdx(ind_nexist(i),i3,j3,k3)=1.0d0
!!$              else
!!$                 facdx(ind_nexist(i),i3,j3,k3)=1.5d0
!!$              endif
              uloc (ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
              uold_loc (ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
              if(compute==2 .and. ivar==2)uloc(ind_nexist(i),i3,j3,k3,ivar)=0.0! neighbor cell at ilevel-1, put value to zero in vector because it is considered as a Dirichlet BC, i.e. in the RHS
           end do
        end do
        
        do i=1,nbuffer
           if(interpol_type_cond.gt.0)then
              facdx(ind_nexist(i),i3,j3,k3)=1.0d0
           else
              facdx(ind_nexist(i),i3,j3,k3)=1.5d0
           endif
        end do

      if(alfven_diff_coeff)then
        ! Compute Alfvenic Mach number = V/V_A and store it into uloc(:,3)
        do i=1,nexist
              dens   = uold(ind_cell(i),1)
              vx     = uold(ind_cell(i),2)/dens
              vy     = uold(ind_cell(i),3)/dens
              vz     = uold(ind_cell(i),4)/dens
              vnorm  = (vx**2+vy**2+vz**2)**0.5
              bx     = 0.5*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
              by     = 0.5*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
              bz     = 0.5*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
              bnorm  = (bx**2+by**2+bz**2)**0.5
              
              va = bnorm/(dens**0.5)
              Ma = vnorm/va
              !==========================================================
              ! Ma>1 : lambda = Linj/Ma^3, isotrop, D=c*lambda/3
              ! Ma<1 : lambda_perp = lambda_para*Ma^4 and lambda_para=Linj/Ma^2
              !==========================================================
              if(Ma > 1)then
                 Dpara = kpar/(Ma**3)
                 kperp = 1.0d0
              else
                 Dpara = kpar*(Ma**2)
                 kperp=Ma**4
              end if
           !else
              Dpara=Dcr/scale_kappa
              kperp=k_perp
           !endif
           uloc(ind_exist(i),i3,j3,k3,3)=Dpara
           uloc(ind_exist(i),i3,j3,k3,4)=kperp ! Dperp=kperp*Dpara
        end do
        do i=1,nbuffer
           if(alfven_diff_coeff)then
              dens= u2(i,ind_son,1)
              vx = u2(i,ind_son,2)/dens
              vy = u2(i,ind_son,3)/dens
              vz = u2(i,ind_son,4)/dens
              vnorm  = (vx**2+vy**2+vz**2)**0.5
              bx = 0.5*(u2(i,ind_son,6)+u2(i,ind_son,nvar+1))
              by = 0.5*(u2(i,ind_son,7)+u2(i,ind_son,nvar+2))
              bz = 0.5*(u2(i,ind_son,8)+u2(i,ind_son,nvar+3))
              bnorm = (bx**2+by**2+bz**2)**0.5
              va = bnorm/dens**0.5
              Ma = vnorm/va
              !==========================================================
              ! Ma>1 : lambda = Linj/Ma^3, isotrop, D=c*lambda/3
              ! Ma<1 : lambda_perp = lambda_para*Ma^4 and lambda_para=Linj/Ma^2
              !==========================================================
              if(Ma > 1)then
                 Dpara = kpar/Ma**3
                 kperp = 1.0d0
              else
                 Dpara = kpar*Ma**2
                 kperp=Ma**4
              end if
              
           else
              Dpara=Dcr/scale_kappa
              kperp=k_perp
           end if

           uloc(ind_nexist(i),i3,j3,k3,3)=Dpara
           uloc(ind_nexist(i),i3,j3,k3,4)=kperp ! Dperp=kperp*Dpara
        end do
      endif
     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids

  if(variable_diff_coeff)then
   ! Compute diffusion coefficient due to streaming instability and store it into uloc(:,5,6)

#if NDIM==1
   ! Yo: It does not work with NCR>1
   j1=1;k1=1
   do l=1,ncache
   do i1=0,3
       dens   = uold_loc(l,i1,j1,k1,1)
       vx     = uold_loc(l,i1,j1,k1,2)/dens
       vy     = uold_loc(l,i1,j1,k1,3)/dens
       vz     = uold_loc(l,i1,j1,k1,4)/dens
       bx     = 0.5*(uold_loc(l,i1,j1,k1,6) + uold_loc(l,i1,j1,k1,nvar+1))
       by     = 0.5*(uold_loc(l,i1,j1,k1,7) + uold_loc(l,i1,j1,k1,nvar+2))
       bz     = 0.5*(uold_loc(l,i1,j1,k1,8) + uold_loc(l,i1,j1,k1,nvar+3))

       bnorm  = sqrt(bx**2+by**2+bz**2)

       va = bnorm/(dens**0.5)
       
       ! compute temperature
       ekin=0.5*dens*(vx**2+vy**2+vz**2)
       emag=0.0d0
       do idim=1,3
           bb=uold_loc(l,i1,j1,k1,5+idim)+uold_loc(l,i1,j1,k1,nvar+idim)
           emag=emag+bb*bb
       end do
       emag=emag*0.125d0
       erad=0.0d0
#if NENER>0
       do irad=1,nener
           erad=erad+uold_loc(l,i1,j1,k1,8+irad)
       end do
#endif
       
       eth=uold_loc(l,i1,j1,k1,5)-ekin-emag-erad
       
       temp=eth*(gamma-1.0d0)/dens*scale_t2*mu_gas    ! Warning T!!!
       
       gradPcr_x=(uold_loc(l,i1+1,j1,k1,igroup)-uold_loc(l,i1-1,j1,k1,igroup))/twodx
       gradPcr = (gradPcr_x*bx)*(gamma_rad(irad)-1.0d0)/bnorm

       call subgridcr_diffusion(dens, Temp, bnorm, gradPcr, kparasub, kperpsub)

       uloc(l,i1,j1,k1,nvar+5)=kparasub
       uloc(l,i1,j1,k1,nvar+6)=kperpsub

       !print *,"kparasub = ",kparasub," kperpsub = ",kperpsub
   enddo
   end do
#endif  
#if NDIM==2
   ! Yo: It does not work with NCR>1
   k1=1
   do l=1,ncache
   do j1=0,3
       do i1=0,3
           dens   = uold_loc(l,i1,j1,k1,1)
           vx     = uold_loc(l,i1,j1,k1,2)/dens
           vy     = uold_loc(l,i1,j1,k1,3)/dens
           vz     = uold_loc(l,i1,j1,k1,4)/dens
           bx     = 0.5*(uold_loc(l,i1,j1,k1,6) + uold_loc(l,i1,j1,k1,nvar+1))
           by     = 0.5*(uold_loc(l,i1,j1,k1,7) + uold_loc(l,i1,j1,k1,nvar+2))
           bz     = 0.5*(uold_loc(l,i1,j1,k1,8) + uold_loc(l,i1,j1,k1,nvar+3))
           
           bnorm  = sqrt(bx**2+by**2+bz**2)
           
           va = bnorm/(dens**0.5)
           
           ! compute temperature
           ekin=0.5*dens*(vx**2+vy**2+vz**2)
           emag=0.0d0
           do idim=1,3
               bb=uold_loc(l,i1,j1,k1,5+idim)+uold_loc(l,i1,j1,k1,nvar+idim)
               emag=emag+bb*bb
           end do
           emag=emag*0.125d0
           erad=0.0d0
#if NENER>0
           do irad=1,nener
               erad=erad+uold_loc(l,i1,j1,k1,8+irad)
           end do
#endif
           
           eth=uold_loc(l,i1,j1,k1,5)-ekin-emag-erad
           
           temp=eth*(gamma-1.0d0)/dens*scale_t2*mu_gas    ! Warning T!!!


           
           gradPcr_x=(uold_loc(l,i1+1,j1  ,k1  ,igroup)&
               &    -uold_loc(l,i1-1,j1  ,k1  ,igroup))/twodx
           gradPcr_y=(uold_loc(l,i1  ,j1+1,k1  ,igroup)&
               &    -uold_loc(l,i1  ,j1-1,k1  ,igroup))/twodx
           gradPcr = (gradPcr_x*bx+gradPcr_y*by)*(gamma_rad(irad)-1.0d0)/bnorm
           
           !write(*,*) "Bnorm_outside_before_2 = ",bnorm
           call subgridcr_diffusion(dens, Temp, bnorm, gradPcr, kparasub, kperpsub)
           !write(*,*) "Bnorm_outside_after_2 = ",bnorm

           uloc(l,i1,j1,k1,nvar+5)=kparasub
           uloc(l,i1,j1,k1,nvar+6)=kperpsub

           !print *,"kparasub = ",kparasub," kperpsub = ",kperpsub
           
       end do
   end do
   end do
#endif  
#if NDIM==3
   ! Yo: It does not work with NCR>1
   do l=1,ncache
   do k1=0,3
   do j1=0,3
       do i1=0,3
           dens   = uold_loc(l,i1,j1,k1,1)
           vx     = uold_loc(l,i1,j1,k1,2)/dens
           vy     = uold_loc(l,i1,j1,k1,3)/dens
           vz     = uold_loc(l,i1,j1,k1,4)/dens
           bx     = 0.5*(uold_loc(l,i1,j1,k1,6) + uold_loc(l,i1,j1,k1,nvar+1))
           by     = 0.5*(uold_loc(l,i1,j1,k1,7) + uold_loc(l,i1,j1,k1,nvar+2))
           bz     = 0.5*(uold_loc(l,i1,j1,k1,8) + uold_loc(l,i1,j1,k1,nvar+3))
           
           bnorm  = sqrt(bx**2+by**2+bz**2)


           !write(*,*) "bnorm = ",bnorm," Gauss"
           !write(*,*) "density = ",dens," g/cm^3"
           !write(*,*) "vx = ",vx," cm/s"
           !write(*,*) "vy = ",vy," cm/s"
           !write(*,*) "vz = ",vz," cm/s"
           
           va = bnorm/(dens**0.5)
           
           ! compute temperature
           ekin=0.5*dens*(vx**2+vy**2+vz**2)
           emag=0.0d0
           do idim=1,3
               bb=uold_loc(l,i1,j1,k1,5+idim)+uold_loc(l,i1,j1,k1,nvar+idim)
               emag=emag+bb*bb
           end do
           emag=emag*0.125d0
           erad=0.0d0
#if NENER>0
           do irad=1,nener
               erad=erad+uold_loc(l,i1,j1,k1,8+irad)
           end do
#endif
           
           eth=uold_loc(l,i1,j1,k1,5)-ekin-emag-erad

           !write(*,*) "eth = ",eth,", uold_loc(l,i1,j1,k1,5) = ",uold_loc(l,i1,j1,k1,5),", ekin = ",ekin,", emag = ",emag,", erad = ",erad
           
           temp=eth*(gamma-1.0d0)/dens*scale_t2*mu_gas    ! Warning T!!!

           !write(*,*) "temp : ",temp, eth, scale_t2, mu_gas
           !write(*,*) "Pcr = ",uold_loc(l,i1,j1  ,k1  ,igroup)

           !write (*,*) "igroup = ",igroup," nvar - igroup = ",nvar-igroup
           gradPcr_x=(uold_loc(l,i1+1,j1  ,k1  ,igroup)&
               &    -uold_loc(l,i1-1,j1  ,k1  ,igroup))/twodx
           gradPcr_y=(uold_loc(l,i1  ,j1+1,k1  ,igroup)&
               &    -uold_loc(l,i1  ,j1-1,k1  ,igroup))/twodx
           gradPcr_z=(uold_loc(l,i1  ,j1  ,k1+1,igroup)&
               &    -uold_loc(l,i1  ,j1  ,k1-1,igroup))/twodx
           gradPcr = (gradPcr_x*bx+gradPcr_y*by+gradPcr_z*bz)*(gamma_rad(irad)-1.0d0)/bnorm

           !if (gradPcr .ne. 0) then
           !write(*,*) "gradPcr_x = ",gradPcr_x," gradPcr_y = ",gradPcr_y," gradPcr_z = ",gradPcr_z," gradPcr =",gradPcr
           !end if 

           !write(*,*) "out : ",temp
           !write(*,*) gradPcr
           !write(*,*) "out : ",dens, bnorm
           !write(*,*) "Bnorm_outside_before_3 = ",bnorm
           !write(*,*) gradPcr_x, gradPcr_y, gradPcr_z, uold_loc(l,i1+1,j1  ,k1  ,igroup), uold_loc(l,i1-1,j1  ,k1  ,igroup)
           
           if (gradPcr .ne. 0.) then

           call subgridcr_diffusion(dens, temp, bnorm, gradPcr, kparasub, kperpsub)
           !write(*,*) "Bnorm_outside_after_3 = ",bnorm
           !write(*,*) kparasub
           !write(*,*) "T = ",temp," Kelvins"
           uloc(l,i1,j1,k1,nvar+5)=kparasub
           uloc(l,i1,j1,k1,nvar+6)=kperpsub

           !----------------------------------------------------------------------------------------------------------------------------------
           !write(*,*) "gamma = ",gamma
           !write(*,*) "eth = ",eth
           !write(*,*) "dens = ",dens
           !write(*,*) "scale_t2 = ",scale_t2
           !write(*,*) "mu_gas = ",mu_gas

           !write(*,*) "kparasub = ",kparasub,", scale_kappa = ",scale_kappa
           !write(*,*) "kparasub = ",kparasub*scale_kappa," DCR = ",DCR
           !write(*,*) "gradPcr = ", gradPcr !* scale_d*scale_v**2/scale_l
           
           else 
           ! In the case where the CRs gradient is null, the self-generated diffusion coefficient is 
           ! infinite (high) in order to make the background diffusion coefficient dominate 
               uloc(l,i1,j1,k1,nvar+5)=Dcr*1d4/scale_kappa
               uloc(l,i1,j1,k1,nvar+6)=Dcr*1d4/scale_kappa
           endif 

           !if (kparasub < 0. .or. kperpsub < 0.) then 
           !print *,"kparasub = ",kparasub," kperpsub = ",kperpsub
           !endif
           
           
           
       end do
   end do
   end do
end do
   
#endif  
end if


if(streaming_diffusion)then

#if NDIM==1
! Yo: It does not work with NCR>1
j1=1;k1=1
do l=1,ncache
do i1=0,3
   deCRL=abs(uold_loc(l,i1  ,j1,k1,indcr1)-uold_loc(l,i1-1,j1,k1,indcr1))/dx
   deCRR=abs(uold_loc(l,i1+1,j1,k1,indcr1)-uold_loc(l,i1  ,j1,k1,indcr1))/dx
   deCRC=abs(uold_loc(l,i1+1,j1,k1,indcr1)-uold_loc(l,i1-1,j1,k1,indcr1))/twodx
   deCR=MAX(deCRC,deCRL,deCRR)
   dens   = uold_loc(l,i1,j1,k1,1)
   bx     = 0.5d0*(uold_loc(l,i1,j1,k1,6)+uold_loc(l,i1,j1,k1,nvar+1))
   by     = 0.5d0*(uold_loc(l,i1,j1,k1,7)+uold_loc(l,i1,j1,k1,nvar+2))
   bz     = 0.5d0*(uold_loc(l,i1,j1,k1,8)+uold_loc(l,i1,j1,k1,nvar+3))
   bnorm2  = bx**2+by**2+bz**2
   va = sqrt(bnorm2/dens)

   deCR  =MAX(deCR,va*uold_loc(l,i1,j1,k1,indcr1)/Dmax)
   length=MIN(uold_loc(l,i1,j1,k1,indcr1)/deCR,nlength_str*dx)
   uloc(l,i1,j1,k1,7)=va*length/scale_kappa*gamma_rad(1)
enddo
enddo
#endif     
#if NDIM==2
! Yo: It does not work with NCR>1
k1=1
do l=1,ncache
do j1=0,3
do i1=0,3
   deCR= ( (uold_loc(l,i1+1,j1,k1,ind_cr1)-uold_loc(l,i1-1,j1,k1,ind_cr1))**2 &
       & +(uold_loc(l,i1,j1+1,k1,ind_cr1)-uold_loc(l,i1,j1-1,k1,ind_cr1))**2 &
       & ) / twodx2
   deCR=sqrt(deCR)
   dens   = uold_loc(l,i1,j1,k1,1)
   bx     = 0.5d0*(uold_loc(l,i1,j1,k1,6)+uold_loc(l,i1,j1,k1,nvar+1))
   by     = 0.5d0*(uold_loc(l,i1,j1,k1,7)+uold_loc(l,i1,j1,k1,nvar+2))
   bz     = 0.5d0*(uold_loc(l,i1,j1,k1,8)+uold_loc(l,i1,j1,k1,nvar+3))
   bnorm2  = bx**2+by**2+bz**2
   va = sqrt(bnorm2/dens)

   deCR  =MAX(deCR,va*uold_loc(l,i1,j1,k1,ind_cr1)/Dmax)
   length=MIN(uold_loc(l,i1,j1,k1,ind_cr1)/deCR,nlength_str*dx)
   uloc(l,i1,j1,k1,7)=va*length/scale_kappa*gamma_rad(1)
enddo
enddo
enddo
#endif     
#if NDIM==3
! Yo: It does not work with NCR>1
do l=1,ncache
do k1=0,3
do j1=0,3
do i1=0,3
   deCR= ( (uold_loc(l,i1+1,j1  ,k1  ,ind_cr1)-uold_loc(l,i1-1,j1  ,k1  ,ind_cr1))**2 &
       & +(uold_loc(l,i1  ,j1+1,k1  ,ind_cr1)-uold_loc(l,i1  ,j1-1,k1  ,ind_cr1))**2 &
       & +(uold_loc(l,i1  ,j1  ,k1+1,ind_cr1)-uold_loc(l,i1  ,j1  ,k1-1,ind_cr1))**2 &
       & ) / twodx2
   deCR=sqrt(deCR)
   dens   = uold_loc(l,i1,j1,k1,1)
   bx     = 0.5d0*(uold_loc(l,i1,j1,k1,6)+uold_loc(l,i1,j1,k1,nvar+1))
   by     = 0.5d0*(uold_loc(l,i1,j1,k1,7)+uold_loc(l,i1,j1,k1,nvar+2))
   bz     = 0.5d0*(uold_loc(l,i1,j1,k1,8)+uold_loc(l,i1,j1,k1,nvar+3))
   bnorm2  = bx**2+by**2+bz**2
   va = sqrt(bnorm2/dens)

   deCR  =MAX(deCR,va*uold_loc(l,i1,j1,k1,ind_cr1)/Dmax)
   length=MIN(uold_loc(l,i1,j1,k1,ind_cr1)/deCR,nlength_str*dx)
   uloc(l,i1,j1,k1,nvar+7)=va*length/scale_kappa*gamma_rad(1)
   !write(*,*) uloc(l,i1,j1,k1,7)
   
enddo
enddo
enddo
enddo
#endif     

endif


  !-----------------------------------------------
  ! Compute energy flux due to thermal conduction
  !-----------------------------------------------
  call crdiff_split(uloc,flux,dx,dx,dx,dt_imp,ncache,compute,facdx,igroup)

  !-----------------------------------------------------
  ! update at level ilevel
  !-----------------------------------------------------
  i0=0; j0=0; k0=0
  if(idim==1)i0=1
  if(idim==2)j0=1
  if(idim==3)k0=1
  do k2=k2min,k2max
     do j2=j2min,j2max
        do i2=i2min,i2max
           ind_son=1+i2+2*j2+4*k2
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
           end do
           i3=1+i2
           j3=1+j2
           k3=1+k2
           ! Update conservative variables new state vector
           do i=1,ncache
              if(son(ind_cell(i))==0)then        
                 if(compute==1)then
                    residu_loc(i,i3   ,j3   ,k3   )= 0.0d0!uold(ind_cell(i),igroup)
                 else if(compute==2)then ! compute Ap 
                    residu_loc(i,i3   ,j3   ,k3   )=unew(ind_cell(i),2)
                 else if(compute==3)then
                    residu_loc(i,i3   ,j3   ,k3   )=1.0d0
                 end if
              
              endif
           end do
        end do
     end do
  end do

  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do i=1,ncache
           if(son(ind_cell(i))==0)then
              residu_loc(i,i3   ,j3   ,k3   )=residu_loc(i,i3   ,j3   ,k3   )+ &
                   & (flux(i,i3   ,j3   ,k3   ,5,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,5,idim))
           endif
        end do
     end do
     end do
     end do
  end do

  i0=0; j0=0; k0=0
  if(idim==1)i0=1
  if(idim==2)j0=1
  if(idim==3)k0=1
  do k2=k2min,k2max
     do j2=j2min,j2max
        do i2=i2min,i2max
           ind_son=1+i2+2*j2+4*k2
           iskip=ncoarse+(ind_son-1)*ngridmax
           do i=1,ncache
              ind_cell(i)=iskip+ind_grid(i)
           end do
           i3=1+i2
           j3=1+j2
           k3=1+k2
           ! Update conservative variables new state vector
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 
                 if(compute==1)then
                    unew(ind_cell(i),1) = -residu_loc(i,i3   ,j3   ,k3  ) ! r0
                    unew(ind_cell(i),2) = -residu_loc(i,i3   ,j3   ,k3  ) ! p0
                    if(variable_diff_coeff)then
                     if(alfven_diff_coeff)then
                         Dcrdiff(ind_cell(i),1)=0.
                         Dcrdiff(ind_cell(i),2)=0. !Do somethin !!!!!
                     else
                         Dcrdiff(ind_cell(i),1)=1d0/(1d0/(DCR/scale_kappa)+1d0/uloc(i,i3,j3,k3,nvar+5))*scale_kappa
                         Dcrdiff(ind_cell(i),2)=1d0/(1d0/(k_perp*DCR/scale_kappa)+1d0/uloc(i,i3,j3,k3,nvar+6))*scale_kappa
                         
                         !if (DCR > uloc(i,i3,j3,k3,nvar+5)*scale_kappa) then
                         !write(*,*) "DCR = ",DCR," Dself,par = ",uloc(i,i3,j3,k3,nvar+5)*scale_kappa !---------------------------------------------------------------
                         !endif

                         !Dcrdiff(ind_cell(i),1)=DCR
                         !Dcrdiff(ind_cell(i),2)=k_perp*DCR
                     endif

                     !write(*,*) "Dcrdiff(ind_cell(i),1) = ",Dcrdiff(ind_cell(i),1)!uloc(i,i3,j3,k3,nvar+5)*scale_kappa !---------------------------------------------------------------
                     !write(*,*) "Dcrdiff(ind_cell(i),2) = ",Dcrdiff(ind_cell(i),2)!uloc(i,i3,j3,k3,nvar+6)*scale_kappa !---------------------------------------------------------------
                 end if
                 else if(compute==2)then ! compute Ap 
                    unew(ind_cell(i),3) = residu_loc(i,i3   ,j3   ,k3   )
                 else if(compute==3)then ! Compute 1/A
                    unew(ind_cell(i),4) = 1.0d0!/residu(i)
                 end if
              endif
           end do
        end do
     end do
  end do


end subroutine crdifffine1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_hydro_cond(u1,ind1,u2,nn)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3)::u1
  integer ,dimension(1:nvector,0:twondim)           ::ind1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3)::u2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are:
  ! interpol_var=0: rho, rho u and E
  ! interpol_var=1: rho, rho u and rho epsilon
  ! The interpolation method is:
  ! interpol_type=0 straight injection
  ! interpol_type=1 linear interpolation with MinMod slope
  ! interpol_type=2 linear interpolation with Monotonized Central slope
  ! interpol_type=3 linear interpolation without limiters
  !----------------------------------------------------------
  integer::i,j,ivar,idim,ind,ix,iy,iz,neul=5,irad

  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,0:twondim),save::a
  real(dp),dimension(1:nvector,1:ndim),save::w
  real(dp),dimension(1:nvector),save::ekin,emag,erad
  real(dp),dimension(1:nvector,0:twondim  ,1:6),save::B1
  real(dp),dimension(1:nvector,1:twotondim,1:6),save::B2

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)
  end do

  ! If necessary, convert father total energy into internal energy
  if(interpol_var_cond==1)then
     do j=0,twondim
        ekin(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u1(i,j,idim+1)**2/max(u1(i,j,1),smallr)
           end do
        end do
        emag(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              emag(i)=emag(i)+0.125d0*(u1(i,j,idim+neul)+u1(i,j,idim+nvar))**2
           end do
        end do
        erad(1:nn)=0.0d0
#if NENER>0
        do irad=1,nener
           do i=1,nn
              erad(i)=erad(i)+u1(i,j,8+irad)
           end do
        end do
#endif
        do i=1,nn
           u1(i,j,neul)=u1(i,j,neul)-ekin(i)-emag(i)-erad(i)
        end do
     end do
  end if


  !------------------------------------------------
  ! Loop over cell-centered interpolation variables
  !------------------------------------------------
  do ivar=1,nvar
  if(ivar<=neul.or.ivar>neul+ndim)then

     ! Load father variable
     do j=0,twondim
        do i=1,nn 
           a(i,j)=u1(i,j,ivar)
        end do
     end do

     ! Reset gradient
     w(1:nn,1:ndim)=0.0D0

     ! Compute gradient with chosen limiter
     if(interpol_type_cond==1)call compute_limiter_minmod(a,w,nn)
     if(interpol_type_cond==2)call compute_limiter_central(a,w,nn)
     if(interpol_type_cond==3)call compute_central(a,w,nn)

     ! Interpolate over children cells
     do ind=1,twotondim
        u2(1:nn,ind,ivar)=a(1:nn,0)
        do idim=1,ndim
           do i=1,nn
              u2(i,ind,ivar)=u2(i,ind,ivar)+w(i,idim)*xc(ind,idim)
           end do
        end do
     end do

  end if
  end do
  ! End loop over cell-centered variables

  ! Update cell centered magnetic field also in redundant array
#if NDIM<2
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,2+nvar)=u2(i,ind,2+neul)
     end do
  end do
#endif
#if NDIM<3
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,3+nvar)=u2(i,ind,3+neul)
     end do
  end do
#endif

  !------------------------------------------------
  ! Loop over face-centered interpolation variables
  !------------------------------------------------
  do j=0,twondim
     do i=1,nn
        B1(i,j,1)=u1(i,j,neul+1)
        B1(i,j,2)=u1(i,j,neul+2)
        B1(i,j,3)=u1(i,j,neul+3)
        B1(i,j,4)=u1(i,j,nvar+1)
        B1(i,j,5)=u1(i,j,nvar+2)
        B1(i,j,6)=u1(i,j,nvar+3)
     end do
  end do
  call interpol_mag(B1,ind1,B2,nn)
  do ind=1,twotondim
     do i=1,nn
        u2(i,ind,neul+1)=B2(i,ind,1)
        u2(i,ind,nvar+1)=B2(i,ind,4)
#if NDIM>1        
        u2(i,ind,neul+2)=B2(i,ind,2)
        u2(i,ind,nvar+2)=B2(i,ind,5)
#endif
#if NDIM>2
        u2(i,ind,neul+3)=B2(i,ind,3)
        u2(i,ind,nvar+3)=B2(i,ind,6)
#endif
     end do
  end do

  ! If necessary, convert children internal energy into total energy
  if(interpol_var_cond==1)then
     do ind=1,twotondim
        ekin(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              ekin(i)=ekin(i)+0.5d0*u2(i,ind,idim+1)**2/max(u2(i,ind,1),smallr)
           end do
        end do
        emag(1:nn)=0.0d0
        do idim=1,3
           do i=1,nn
              emag(i)=emag(i)+0.125d0*(u2(i,ind,idim+neul)+u2(i,ind,idim+nvar))**2
           end do
        end do
        erad(1:nn)=0.0d0
#if NENER>0
        do irad=1,nener
           do i=1,nn
              erad(i)=erad(i)+u2(i,ind,8+irad)
           end do
        end do
#endif
       do i=1,nn
           u2(i,ind,neul)=u2(i,ind,neul)+ekin(i)+emag(i)+erad(i)
        end do
     end do
  end if

end subroutine interpol_hydro_cond
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_mag_cond(B1,ind1,B2,nn)
  use amr_commons
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim  ,1:6)::B1
  integer ,dimension(1:nvector,0:twondim)      ::ind1
  real(dp),dimension(1:nvector,1:twotondim,1:6)::B2
  !----------------------------------------------------------
  ! This routine performs a prolongation (interpolation)
  ! operation for newly refined cells or buffer cells.
  ! The interpolated variables are Bx, By and Bz.
  ! Divergence free is garanteed.
  ! The scheme is the one invented by Toth and Balsara.
  ! interpol_mag_type=0: straight injection
  ! interpol_mag_type=1: linear interpolation with MinMod slope
  ! interpol_mag_type=2: linear interpolation with Monotonized Central slope
  ! interpol_mag_type=3: linear interpolation without limiters
  !----------------------------------------------------------
  integer::i,j,k,ind,l,idim,imax,jmax,kmax
  real(dp),dimension(1:nvector,-1:1,0:1,0:1),save::u
  real(dp),dimension(1:nvector,0:1,-1:1,0:1),save::v
  real(dp),dimension(1:nvector,0:1,0:1,-1:1),save::w

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Compute interpolated fine B over coarse side faces
  call interpol_faces_cond(B1,u,v,w,nn)
  
  ! Get fine B from refined faces, if any
  call copy_from_refined_faces(B1,ind1,u,v,w,nn)
 
  ! Compute interpolated fine B inside coarse cell.
  call cmp_central_faces(u,v,w,nn)

  ! Scatter results
  do i=0,imax
  do j=0,jmax
  do k=0,kmax
     ind=1+i+2*j+4*k
     do l=1,nn
        B2(l,ind,1)=u(l,i-1,j,k)
        B2(l,ind,2)=v(l,i,j-1,k)
        B2(l,ind,3)=w(l,i,j,k-1)
        B2(l,ind,4)=u(l,i,j,k)
        B2(l,ind,5)=v(l,i,j,k)
        B2(l,ind,6)=w(l,i,j,k)
     end do
  end do
  end do
  end do

end subroutine interpol_mag_cond
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine interpol_faces_cond(b1,u,v,w,nn)
  use amr_commons
  use hydro_commons, ONLY: interpol_mag_type_cond
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:twondim,1:6)::b1
  real(dp),dimension(1:nvector,-1:1,0:1,0:1)::u
  real(dp),dimension(1:nvector,0:1,-1:1,0:1)::v
  real(dp),dimension(1:nvector,0:1,0:1,-1:1)::w

  ! TVD interpolation from coarse faces
  integer::i,j,k,l,imax,jmax,kmax
  real(dp),dimension(1:nvector,0:4),save::b
  real(dp),dimension(1:nvector,1:2),save::s

  imax=1; jmax=0; kmax=0
#if NDIM>1
  jmax=1
#endif
#if NDIM>2
  kmax=1
#endif

  ! Left face along direction x (interpolate Bx)
  do l=1,nn
     b(l,0)=b1(l,0,1)
  end do
#if NDIM>1     
  do l=1,nn
     b(l,1)=b1(l,3,1)
     b(l,2)=b1(l,4,1)
  end do
#endif
#if NDIM>2
  do l=1,nn
     b(l,3)=b1(l,5,1)
     b(l,4)=b1(l,6,1)
  end do
#endif

  s(1:nn,1:2)=0.0
#if NDIM==2
  if(interpol_mag_type_cond>0)call compute_1d_tvd_cond(b,s,nn)
#endif
#if NDIM==3
  if(interpol_mag_type_cond>0)call compute_2d_tvd_cond(b,s,nn)
#endif
  do j=0,jmax
  do k=0,kmax
     do l=1,nn
        u(l,-1,j,k)=b(l,0)+0.5*s(l,1)*(dble(j)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

  ! Right face along direction x (interpolate Bx)
  do l=1,nn
     b(l,0)=b1(l,0,4)
  end do
#if NDIM>1     
  do l=1,nn
     b(l,1)=b1(l,3,4)
     b(l,2)=b1(l,4,4)
  end do
#endif
#if NDIM>2
  do l=1,nn
     b(l,3)=b1(l,5,4)
     b(l,4)=b1(l,6,4)
  end do
#endif

  s(1:nn,1:2)=0.0
#if NDIM==2
  if(interpol_mag_type_cond>0)call compute_1d_tvd_cond(b,s,nn)
#endif
#if NDIM==3
  if(interpol_mag_type_cond>0)call compute_2d_tvd_cond(b,s,nn)
#endif
  do j=0,jmax
  do k=0,kmax
     do l=1,nn
        u(l,+1,j,k)=b(l,0)+0.5*s(l,1)*(dble(j)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

#if NDIM>1
  ! Left face along direction y (interpolate By)
  do l=1,nn
     b(l,0)=b1(l,0,2)
  end do
  do l=1,nn
     b(l,1)=b1(l,1,2)
     b(l,2)=b1(l,2,2)
  end do
#if NDIM>2     
  do l=1,nn
     b(l,3)=b1(l,5,2)
     b(l,4)=b1(l,6,2)
  end do
#endif

  s(1:nn,1:2)=0.0
#if NDIM==2
  if(interpol_mag_type_cond>0)call compute_1d_tvd_cond(b,s,nn)
#endif
#if NDIM==3
  if(interpol_mag_type_cond>0)call compute_2d_tvd_cond(b,s,nn)
#endif
  do i=0,imax
  do k=0,kmax
     do l=1,nn
        v(l,i,-1,k)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do

  ! Right face along direction y (interpolate By)
  do l=1,nn
     b(l,0)=b1(l,0,5)
  end do
  do l=1,nn
     b(l,1)=b1(l,1,5)
     b(l,2)=b1(l,2,5)
  end do
#if NDIM>2     
  do l=1,nn
     b(l,3)=b1(l,5,5)
     b(l,4)=b1(l,6,5)
  end do
#endif

  s(1:nn,1:2)=0.0
#if NDIM==2
  if(interpol_mag_type_cond>0)call compute_1d_tvd_cond(b,s,nn)
#endif
#if NDIM==3
  if(interpol_mag_type_cond>0)call compute_2d_tvd_cond(b,s,nn)
#endif
  do i=0,imax
  do k=0,kmax
     do l=1,nn
        v(l,i,+1,k)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(k)-0.5)
     end do
  end do
  end do
#endif

#if NDIM>2
  ! Left face along direction z (interpolate Bz)
  do l=1,nn
     b(l,0)=b1(l,0,3)
     b(l,1)=b1(l,1,3)
     b(l,2)=b1(l,2,3)
     b(l,3)=b1(l,3,3)
     b(l,4)=b1(l,4,3)
  end do

  s(1:nn,1:2)=0.0
  if(interpol_mag_type_cond>0)call compute_2d_tvd_cond(b,s,nn)
  do i=0,1
     do j=0,1
        do l=1,nn
           w(l,i,j,-1)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(j)-0.5)
        end do
     end do
  end do

  ! Right face along direction z (interpolate Bz)
  do l=1,nn
     b(l,0)=b1(l,0,6)
     b(l,1)=b1(l,1,6)
     b(l,2)=b1(l,2,6)
     b(l,3)=b1(l,3,6)
     b(l,4)=b1(l,4,6)
  end do

  s(1:nn,1:2)=0.0
  if(interpol_mag_type_cond>0)call compute_2d_tvd_cond(b,s,nn)
  do i=0,1
     do j=0,1
        do l=1,nn
           w(l,i,j,+1)=b(l,0)+0.5*s(l,1)*(dble(i)-0.5)+0.5*s(l,2)*(dble(j)-0.5)
        end do
     end do
  end do
#endif

end subroutine interpol_faces_cond
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_2d_tvd_cond(b,s,nn)
  use amr_commons, ONLY: nvector
  use hydro_commons, ONLY: interpol_mag_type_cond
  use const
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:4)::b
  real(dp),dimension(1:nvector,1:2)::s
  
  integer::i
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop

  if(interpol_mag_type_cond==3)then
     do i=1,nn
        dlft = half*(b(i,0) - b(i,1))
        drgt = half*(b(i,2) - b(i,0))
        s(i,1) = dlft+drgt
     end do
     do i=1,nn
        dlft = half*(b(i,0) - b(i,3))
        drgt = half*(b(i,4) - b(i,0))
        s(i,2) = dlft+drgt
     end do
     return
  endif

  do i=1,nn
     dlft = interpol_mag_type_cond*(b(i,0) - b(i,1))
     drgt = interpol_mag_type_cond*(b(i,2) - b(i,0))
     dcen = half*(dlft+drgt)/interpol_mag_type_cond
     dsgn = sign(one, dcen)
     slop = min(abs(dlft),abs(drgt))
     dlim = slop
     if((dlft*drgt)<=zero)dlim=zero
     s(i,1) = dsgn*min(dlim,abs(dcen))
  end do

  do i=1,nn
     dlft = interpol_mag_type_cond*(b(i,0) - b(i,3))
     drgt = interpol_mag_type_cond*(b(i,4) - b(i,0))
     dcen = half*(dlft+drgt)/interpol_mag_type_cond
     dsgn = sign(one, dcen)
     slop = min(abs(dlft),abs(drgt))
     dlim = slop
     if((dlft*drgt)<=zero)dlim=zero
     s(i,2) = dsgn*min(dlim,abs(dcen))
  end do


end subroutine compute_2d_tvd_cond
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine compute_1d_tvd_cond(b,s,nn)
  use amr_commons, ONLY: nvector
  use hydro_commons, ONLY: interpol_mag_type_cond
  use const
  implicit none
  integer::nn
  real(dp),dimension(1:nvector,0:4)::b
  real(dp),dimension(1:nvector,1:2)::s
  
  integer::i
  real(dp)::dsgn, dlim, dcen, dlft, drgt, slop

  if(interpol_mag_type_cond==3)then
     do i=1,nn
        dlft = half*(b(i,0) - b(i,1))
        drgt = half*(b(i,2) - b(i,0))
        s(i,1) = dlft+drgt
     end do
     return
  endif
  do i=1,nn
     dlft = interpol_mag_type_cond*(b(i,0) - b(i,1))
     drgt = interpol_mag_type_cond*(b(i,2) - b(i,0))
     dcen = half*(dlft+drgt)/interpol_mag_type_cond
     dsgn = sign(one, dcen)
     slop = min(abs(dlft),abs(drgt))
     dlim = slop
     if((dlft*drgt)<=zero)dlim=zero
     s(i,1) = dsgn*min(dlim,abs(dcen))
  end do

end subroutine compute_1d_tvd_cond
!###########################################################
!###########################################################
!###########################################################
!###########################################################

! -------------- Subroutines for the calculation of the sub-resolution diffusion coefficient ----------------
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine subgridcr_diffusion(rho_sim, T_sim, B0_sim, gradPcr, kpara, kperp)
   use cooling_module
   
   ! - Some constants ---------------------------------------
   real::mp, e, GeV, pi 
   
   ! - Variables for NN interpolation -----------------------
   integer::k, n, kvalues  
   real(dp)::X_test, Y_test, rho_test, T_test
   real(dp)::calc_value1, calc_value2, calc_value3, betacr 
   real(dp),dimension(1:5)::mi,ni,mn,nn,Xion,T,rho
   real(dp),dimension(1:5)::X_train, Y_train, distance  
   real(dp),dimension(1:5):: f_value1, f_value2, f_value3
   real(dp),dimension(1:3)::dist_k, x_k, y_k, val1, val2, val3
   
   ! - Variables from RAMSES-ISM ----------------------------
   real(dp)::rho_sim, T_sim, B0_sim, gradPcr
   real(dp)::ECR_sim
   ! - Variables for RAMSES-ISM -----------------------------
   real(dp)::Gamma_in,rho_i,rho_n,rmi,rmn,chi,rXion,xi_n,rnn,rni,rg
   real(dp)::nuin,Va,Vai,Va_temp, Va_final,Turb
   parameter(k_reel = selected_real_kind(6,90))
   real(dp),intent(out)::kpara,kperp
   real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,scale_omega,scale_tau,scale_kappa
   
   call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
   scale_kappa = scale_l**2/scale_t
   scale_m = scale_d*scale_l**3
   scale_n = scale_l**3

   !write(*,*) "scale_l = ",scale_l
   !write(*,*) "scale_t = ",scale_t
   !write(*,*) "scale_d = ",scale_v

   
   ! - Some constants 
   mp = 1.6726d-24 ! Proton mass [g]
   e = 4.80326d-10 ! Electron charge [statcoulomb]
   GeV = 0.00160218 ! GeV in Erg
   pi = acos(-1.0d0)
   
   rho_sim = rho_sim * scale_d
   gradPcr = gradPcr * scale_d*scale_v**2/scale_l
   
   !if (gradPcr .ne. 0.) then 
   !    write(*,*) "gradPcr = ",gradPcr," erg/cm^3"
   !endif
   
   
   !write(*,*) "Bnorm_inside_before = ",B0_sim
   B0_sim  = B0_sim  * sqrt(4.*pi*scale_d*scale_v**2 )
   !write(*,*) "pi = ",pi," scale_d = ",scale_d," scale_v = ",scale_v
   !write(*,*) "Bnorm_inside_after = ",B0_sim
   
   !write(*,*) "In : ",T_sim
   
   
   
   ! - kNN interpolation algorithm parameters ----------------
   k = 5  ! Number of nearest neighbors 
   n = 20 ! Distance ponderation parameter 
   
   ! - Variables from RAMSES-ISM ----------------------------
   !!$rho_sim = 1.d-21 ! g.cm^-3 
   !!$T_sim   = 50.    ! Kelvin
   !!$B0_sim = 10.d-6   ! Gauss 
   ECR_sim = 1.*GeV ! Erg     
   !!$gradPcr = 1d-29  ! Erg/cm^4
   
   ! - Observationnal and simulation values ------------------
   T    = [  8000.,     50.,     50.,     30.,     20.]
   mi   = [ 1.0*mp, 12.0*mp, 12.0*mp, 29.0*mp, 29.0*mp]
   mn   = [1.21*mp, 1.21*mp, 1.67*mp, 2.12*mp, 2.12*mp]
   ni   = [  0.007,  2.4e-2,    0.03,    3e-2,    0.01]
   Xion = [   0.02,    8e-4,    1e-4,    1e-5,    1e-6]
   
   do kvalues=1,5
       rho(kvalues) = (mi(kvalues)*ni(kvalues) + mn(kvalues)*ni(kvalues)*(Xion(kvalues)**(-1) - 1.))
       X_train(kvalues) = log10(T(kvalues))
       Y_train(kvalues) = log10(rho(kvalues))
   end do 
   
   X_test = log10(T_sim)
   Y_test = log10(rho_sim)
   
   ! - Algorithm ---------------------------------------------
   do kvalues=1,5
       f_value1(kvalues) = log10(Xion(kvalues))
       f_value2(kvalues) = mn(kvalues)
       f_value3(kvalues) = mi(kvalues)
   end do 
   
   !write(*,*) X_test
   call getValue_linear(X_train, Y_train, distance, f_value1, f_value2, f_value3, X_test, Y_test, k , dist_k, x_k, y_k, val1, val2, val3, n, calc_value1, calc_value2, calc_value3)
   !write(*,*) X_test
   rXion = 10.**(calc_value1) ! From log axis to lin axis 
   rmn = calc_value2
   rmi = calc_value3
   
   !write(*,*) rXion, rmn/mp, rmi/mp
   
   ! - Basic calculations 
   chi   = (rmn/rmi)*(rXion**(-1) - 1.) 
   rho_n = rho_sim/(1. + chi**(-1))
   rho_i = chi**(-1)/(1. + chi**(-1))*rho_sim
   xi_n  = 1./(1 + chi**(-1))
   rnn   = rho_n/rmn 
   rni   = rho_i/rmi 
   rg    = ECR_sim/(e*B0_sim)
   
   !write (*,*) "ECR_sim = ",ECR_sim," e = ",e," B0_sim = ",B0_sim 
   
   ! - nu_in and Gamma_in calculation 
   Gamma_in = 0.
   nuin  = 2*rnn*8.4d-9*(T_sim/1.d4)**(0.4) 
   Va  = B0_sim/sqrt(4.*pi*rho_sim)
   Vai = B0_sim/sqrt(4.*pi*rho_i)
   Va_temp = B0_sim/sqrt(4.*pi*rho_i)
   !if (ECR_sim < e*B0_sim*Vai/nuin) then 
   !    Gamma_in = - nuin/2. 
   !    Va_final = Vai 
   !end if 
   Va_temp = B0_sim/(4.*pi*rho_sim)
   !if (ECR_sim > e*B0*Va*chi/nuin) then 
   !    Gamma_in = -xi_n*Va**2*e**2*B0_sim**2*ECR_sim**(-2)/(2.*chi*nuin)
   !    Va_final = Va
   !end if 
   
   Gamma_in = -nuin/2.
   Va_final = Vai
   
   ! - Energy density of waves (Turbulence level)
   Turb = abs(Va_final*abs(gradPcr)/(2.*Gamma_in)/(B0_sim**2/(8.*pi)))
   !write(*,*) "Turb = ",Turb," Va_final = ",Va_final," gradPcr = ",gradPcr," Gamma_in = ",Gamma_in," B0_sim = ",B0_sim," pi = ",pi
   !write(*,*) Turb
   
   !write(*,*) rnn, rho_n, rmn/mp
   !write(*,*) e*B0*Va*chi/nuin/GeV, ECR_sim/GeV, e*B0_sim*Vai/nuin/GeV 
   !write(*,*) "Va =",Va, "Vai = ",Vai, "Va_final = ",Va_final
   
   ! - Diffusion coefficients 
   !write(*,*)  4.*pi*rg*betacr*clight/(3.*Turb)
   betacr = sqrt(1 - (1/(1 + ECR_sim/(mp*clight)**2))**2)
   kpara = 4.*pi*rg*betacr*clight/(3.*Turb) / scale_kappa
   
   !write(*,*) "betacr = ",betacr
   !kperp = kpara*Turb**2
   kperp = kpara*1.d-2 ! In order to avoid scheme divergences 
   
   !write(*,*) "-----------------------------------------------------------------"
   !write(*,*) "Tested values : rho = ",rho_sim," g.cm^-3, T = ",T_sim," Kelvin" 
   !write(*,*) "Result : Xion = ",rXion," mi = ",rmi/mp," mp, mn = ",rmn/mp," mp" 
   !write(*,*) "chi = ",chi,", rho_n = ",rho_n," g.cm^-3, rho_i = ",rho_i," g.cm^-3" 
   !write(*,*) "Xi_n = ",xi_n,", nn = ",rnn," cm^-3, ni = ",rni," cm^-3" 
   !write(*,*) "rg = ",rg," cm, nu_in = ",nuin," s^-1, Va_final",Va_final," cm.s^-1" 
   !write(*,*) "Gamma_in = ",Gamma_in," s^-1, I(k) = ",Turb
   !write(*,*) "kpara = ",kpara*scale_kappa," cm^2.s^-1, kperp = ",kperp*scale_kappa," cm^2.s^-1" 
   !write(*,*) "-----------------------------------------------------------------"
   
   end subroutine subgridcr_diffusion
   !###########################################################
   !###########################################################
   !###########################################################
   !###########################################################
   subroutine getValue_linear(X_train, Y_train, distance, value1, value2, value3, X_test, Y_test, k , dist_k, x_k, y_k, val1, val2, val3, n, calc_value1, calc_value2, calc_value3)
       use cooling_module
   implicit none
   ! - Input values -----------------------------------------
   integer::k, n 
   real(dp),intent(in)::X_test, Y_test
   real(dp),dimension(1:5)::X_train, Y_train, value1, value2, value3, distance 
   real(dp),dimension(1:k)::x_k, y_k, val1, val2, val3, dist_k
   
   ! - Output values ----------------------------------------
   real(dp)::calc_value1, calc_value2, calc_value3 
   
   ! - InCode values ----------------------------------------
   real(dp)::distinv
   real(dp)::temp_value1
   real(dp)::temp_value2
   real(dp)::temp_value3
   integer::ivalues
   
   distinv = 0.0D0 
   temp_value1 = 0.0D0
   temp_value2 = 0.0D0
   temp_value3 = 0.0D0 
   
   
   !write(*,*) "Valeur !! : ",X_test
   !write(*,*) "Before : ",calc_value1, calc_value2, calc_value3
   
   call selectKNN(X_train, Y_train, X_test, Y_test, k, value1, value2, value3, dist_k, x_k, y_k, val1, val2, val3)
   
   
   
   if (k.EQ.1) then 
       calc_value1 = val1(1)
       calc_value2 = val2(1)
       calc_value3 = val3(1)
   end if 
   if (k > 1) then 
       do ivalues=1,k 
           distinv = distinv + 1./(dist_k(ivalues)**n) 
           temp_value1 = temp_value1 + val1(ivalues)/(dist_k(ivalues)**n)
           temp_value2 = temp_value2 + val2(ivalues)/(dist_k(ivalues)**n)
           temp_value3 = temp_value3 + val3(ivalues)/(dist_k(ivalues)**n)
       end do 
   end if 
   
   
   calc_value1 = temp_value1/distinv 
   calc_value2 = temp_value2/distinv 
   calc_value3 = temp_value3/distinv 
   
   !write(*,*) "After : ",calc_value1, calc_value2, calc_value3
   
   
   
   end subroutine getValue_linear 
   !###########################################################
   !###########################################################
   !###########################################################
   !###########################################################
   subroutine selectKNN(X_train, Y_train, X_test, Y_test, k, value1, value2, value3, dist_k, x_k, y_k, val1, val2, val3) 
       use cooling_module
   implicit none
   ! - Input values ------------------------------------------
   integer::k, ivalues 
   real(dp),dimension(1:5)::X_train, Y_train, value1, value2, value3, distance  
   real(dp)::X_test, Y_test
   
   ! - Output values ----------------------------------------- 
   real(dp),dimension(1:k)::dist_k, x_k, y_k, val1, val2, val3 
   
   
   call getDistance(X_train, Y_train, X_test, Y_test, value1, value2, value3, distance)
   
   do ivalues=1,k 
       dist_k(ivalues) = distance(ivalues)
       x_k(ivalues) = X_train(ivalues)
       y_k(ivalues) = Y_train(ivalues)
       val1(ivalues) = value1(ivalues)
       val2(ivalues) = value2(ivalues)
       val3(ivalues) = value3(ivalues)
   
       !write(*,*) dist_k(ivalues)
       !write(*,*) val1(ivalues), val2(ivalues), val3(ivalues)
   end do 
   
   end subroutine selectKNN
   !###########################################################
   !###########################################################
   !###########################################################
   !###########################################################
   subroutine getDistance(X_train, Y_train, X_test, Y_test, value1, value2, value3, distance) 
       use cooling_module
   implicit none 
   ! - Input values ------------------------------------------
   integer::ivalues,jvalues
   real(dp),dimension(1:5)::X_train, Y_train, value1, value2, value3  
   real(dp)::X_test,Y_test
   
   ! - Output values -----------------------------------------
   real(dp),dimension(1:5)::distance
   integer,dimension(1:5)::order 
   real(dp)::buffer_value
   
   do ivalues = 1,5 
       distance(ivalues) = sqrt((X_train(ivalues)-X_test)**2. + (Y_train(ivalues)-Y_test)**2.)
       !write(*,*) ivalues," : ",distance(ivalues)
       !write(*,*) "Tsim = ",X_test,"  Ttrain = ",X_train(ivalues)
       !write(*,*) "rhosim = ",Y_test," rhotrain = ",Y_train(ivalues)
       !write(*,*) distance(ivalues)
       order(ivalues) = ivalues
   end do
   
   do ivalues = 1,5
       do jvalues = 1,5
           if (ivalues.NE.jvalues) then 
               if ((distance(ivalues) <= distance(jvalues)).AND.(ivalues > jvalues)) then
   
                   buffer_value = distance(ivalues)
                   distance(ivalues) = distance(jvalues)
                   distance(jvalues) = buffer_value
   
                   buffer_value = X_train(ivalues)
                   X_train(ivalues) = X_train(jvalues)
                   X_train(jvalues) = buffer_value
   
                   buffer_value = Y_train(ivalues)
                   Y_train(ivalues) = Y_train(jvalues)
                   Y_train(jvalues) = buffer_value
   
                   buffer_value = order(ivalues)
                   order(ivalues) = order(jvalues)
                   order(jvalues) = buffer_value
               end if 
           endif
       end do
   end do
   end subroutine getDistance