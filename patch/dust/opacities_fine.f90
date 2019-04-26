subroutine set_opacities(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the dust diffusion scheme.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated. 
  !--------------------------------------------------------------------------
  integer::i,ivar,igrid,ncache,ngrid,ind,iskip,icpu,idust,idim
  integer,dimension(1:nvector),save::ind_grid
  logical:: d_cycle_ok
  integer :: icycle, ncycle

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call opafine1(ind_grid,ngrid,ilevel)
  end do
  
   if(simple_boundary)call make_boundary_hydro(ilevel)


111 format('   Entering dust_diffusion_fine for level ',i2)


end subroutine set_opacites
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine opafine1(ind_grid,ncache,ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  use cooling_module,ONLY:kB,mH
  use cloud_module
  use radiation_parameters
  use units_commons, only : scale_m


  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !---------------------------------------------------------------------!
  ! This routine gathers first hydro variables from neighboring grids  -!
  ! to set initial conditions in a 6x6x6 grid. It interpolate from     -!
  ! coarser level missing grid variables. It then calls the            -!
  ! dust diffusion solver that compute the flux. This flux is zeroed at-!
  ! coarse-fine boundaries, since contribution from finer levels has   -!
  ! already been taken into account. Conservative variables are updated-!
  ! and stored in array unew(:), both at the current level and at the  -!
  ! coarser level if necessary.                                        -!
  !---------------------------------------------------------------------!

  real(dp),dimension(1:nvector,0:twondim  ,1:nvar+3),save::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar+3),save::u2
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3),save::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:2),save::opaloc=0.0d0

  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2),save::ok
  integer ,dimension(1:nvector,1:threetondim     ),save::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       ),save::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         ),save::ibuffer_father
  integer ,dimension(1:nvector,0:twondim         ),save::ind1
  integer ,dimension(1:nvector                   ),save::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::idust,ht
  integer::i,j,ivar,idim,irad,ind_son,ind_father,iskip,nbuffer,ibuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  real(dp)::sum_dust,sum_dust_new,sum_dust_old
  real(dp)::d,u,v,w,A,B,C,enint,e_kin,e_mag,pressure,cs, temp
  real(dp)::rho_gas, pi, t_stop,t_stop_floor,dens_floor,d0,r0

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

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
        call interpol_hydro(u1,ind1,u2,nbuffer)
     endif

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
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

  

        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do
  ! End loop over neighboring grids


  call cmpopa(uloc,opaloc,dx,dx,dx,dtnew(ilevel),ncache)

   !--------------------------------------------------------
   !Udate at level ilevel for the dust velocity
   !--------------------------------------------------------
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
        do i=1,ncache
           if(son(ind_cell(i))==0)then
              do iopa=1,2
                 opacities_PR(ind_cell(i),idust,idim)=opaloc(i,i3,j3,k3,iopa)
              enddo
           end if
     end do
  end do
  end do
  end do
end do


 
end subroutine opafine1


!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine cmpopa(uin,opout,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use hydro_commons
  use units_commons
  use const
  use cloud_module
  use cooling_module,ONLY:kB,mH,clight
  use radiation_parameters  
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:2)::opout
  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::qin 
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
    real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 

  ! declare local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::idust,idim, ifreq
  real(dp)  ::pi, sum_dust,kappagas
  real(dp), dimension(100) :: nufreq,ones
  real(dp), dimension (100):: kappa_nu
  real(dp), dimension(1:ndust) ::d_grain,l_grain,isnot_charged
  epsilon_0 = dust_ratio(1)
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  pi =3.14159265358979323846_dp
  if(mrn.eqv..true.) then
     isnot_charged=0.0d0    
     call size_dust(l_grain)
     do idust=1,ndust
       l_grain(idust) = l_grain(idust)
       d_grain(idust)=grain_dens(idust)
    end do
  else
     do idust=1,ndust
       d_grain(idust)=grain_dens(idust)
       l_grain(idust)=grain_size(idust)

    end do
 endif
 kappagas=0.01d0
 do ifreq = 1, 100
    ones(ifreq)= 1.0d0
    nufreq(ifreq) = 10.0d0**(log10(numax/numin)*real(ifreq,dp)/real(100,dp)+log10(numin))
 enddo
 call ctoprimopa(uin,qin,bf,gravin,dt,ngrid)

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
              kappa_nu=0.0d0
              sum_dust=0.0d0
#if NDUST>0              
              do idust =1,ndust
                 sum_dust= sum_dust + qin(l,i,j,k,firstindex_ndust+idust)
                 kappa_nu=kappa_nu+qin(l,i,j,k,firstindex_ndust+idust)/((4.0d0/3.0d0)*pi*d_grain(idust)*l_grain(idust)**3.0)*merge(2.0*pi*l_grain(idust)*nufreq/clight,ones,2.0*pi*l_grain(idust)*nufreq/clight<ones)
              enddo
#endif
              kappa_nu = kappa_nu + (1.0-sum_dust)*kappagas !note that this is cgs units
        end do
     end do
  end do
end do
end do


end subroutine cmpopa

subroutine ctoprimopa(uin,q,bf,gravin,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use const
  use radiation_parameters,only:small_er
  implicit none

  integer ::ngrid
  real(dp)::dt
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::q  
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf  

  integer ::i, j, k, l, idim
  real(dp)::eint, smalle, smallp, etot
  real(dp),dimension(1:nvector),save::eken,emag,erad

  ! EOS
  real(dp)  :: pp_eos

#if NENER>0
  integer::irad
#endif
#if NVAR>8+NENER
  integer::n
#endif
  real(dp):: sum_dust
#if NDUST>0
  integer:: idust
#endif  
  
  smalle = smallc**2/gamma/(gamma-one)
  smallp = smallr*smallc**2/gamma

  ! Store face centered magnetic field
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2+1
           DO l = 1, ngrid
              if(i<=iu2)then
                 bf(l,i,j,k,1) = uin(l,i,j,k,6)
              else
                 bf(l,i,j,k,1) = uin(l,i-1,j,k,nvar+1)
              endif
           END DO
        end do
     end do
  end do
  do k = ku1, ku2
     do j = ju1, ju2+1
        do i = iu1, iu2
           DO l = 1, ngrid
              if(j<=ju2)then
                 bf(l,i,j,k,2) = uin(l,i,j,k,7)
              else
                 bf(l,i,j,k,2) = uin(l,i,j-1,k,nvar+2)
              endif
           END DO
        end do
     end do
  end do
  do k = ku1, ku2+1
     do j = ju1, ju2
        do i = iu1, iu2
           DO l = 1, ngrid
              if(k<=ku2)then
                 bf(l,i,j,k,3) = uin(l,i,j,k,8)
              else
                 bf(l,i,j,k,3) = uin(l,i,j,k-1,nvar+3)
              endif
           END DO
        end do
     end do
  end do

  ! Convert to primitive variable
  do k = ku1, ku2
     do j = ju1, ju2
        do i = iu1, iu2

           ! Compute density
           do l = 1, ngrid
              q(l,i,j,k,1) = max(uin(l,i,j,k,1),smallr)
           end do
           ! Debug
           if(debug)then
              do l = 1, ngrid
                 if(uin(l,i,j,k,1).le.smallr)then
                    write(*,*)'negative density'
                    write(*,*)uin(l,i,j,k,1)
                    stop
                 end if
              end do
           end if

           ! Compute velocities
           do l = 1, ngrid
              q(l,i,j,k,2) = uin(l,i,j,k,2)/q(l,i,j,k,1)
              q(l,i,j,k,3) = uin(l,i,j,k,3)/q(l,i,j,k,1)
              q(l,i,j,k,4) = uin(l,i,j,k,4)/q(l,i,j,k,1)
           end do

           ! Compute cell centered magnetic field
           DO l = 1, ngrid
              q(l,i,j,k,6) = (uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))*half
              q(l,i,j,k,7) = (uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))*half
              q(l,i,j,k,8) = (uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))*half
           END DO

           ! Compute specific kinetic energy and magnetic energy
           do l = 1, ngrid
              eken(l) = half*(q(l,i,j,k,2)**2+q(l,i,j,k,3)**2+q(l,i,j,k,4)**2)
              emag(l) = half*(q(l,i,j,k,6)**2+q(l,i,j,k,7)**2+q(l,i,j,k,8)**2)
           end do

           ! Compute non-thermal pressure
           erad = zero
#if NENER>0
           do irad = 1,nent
              do l = 1, ngrid
                 q(l,i,j,k,8+irad) = (gamma_rad(irad)-one)*uin(l,i,j,k,8+irad)
                 erad(l) = erad(l)+uin(l,i,j,k,8+irad)
              end do
           enddo
           do irad = 1,ngrp
              do l = 1, ngrid
                 q(l,i,j,k,firstindex_er+irad) = uin(l,i,j,k,firstindex_er+irad)
                 erad(l) = erad(l)+uin(l,i,j,k,firstindex_er+irad)
              end do
           enddo
#endif

  ! Passive scalar (and extinction and internal energy and rad fluxes in M1) !!!!!!
           
           ! Compute thermal pressure through EOS
           do l = 1, ngrid
              sum_dust=0.0d0
#if NDUST>0              
              do idust = 1, ndust
                 sum_dust=sum_dust+q(l,i,j,k,firstindex_ndust+idust)
              end do
#endif  
              etot = uin(l,i,j,k,5) - emag(l) -erad(l)
              eint = etot-eken(l)*q(l,i,j,k,1)
              if(energy_fix)eint=uin(l,i,j,k,nvar)
            
              call pressure_eos((1.0d0-sum_dust)*q(l,i,j,k,1),eint,pp_eos)
              q(l,i,j,k,5)=MAX(pp_eos,smallp)
           end do

           ! Gravity predictor step
           do idim = 1, ndim
              do l = 1, ngrid
                 q(l,i,j,k,idim+1) = q(l,i,j,k,idim+1) !+ gravin(l,i,j,k,idim)*dt*half
              end do
           end do

        end do
     end do
  end do
#if NVAR>8+NENER
  do n = firstindex_pscal+1, firstindex_pscal+npscal
     do k = ku1, ku2
        do j = ju1, ju2
           do i = iu1, iu2
              do l = 1, ngrid
                 q(l,i,j,k,n) = uin(l,i,j,k,n)/max(uin(l,i,j,k,1),smallr)
              end do
           end do
        end do
     end do
  end do
#endif
 
end subroutine ctoprimopa
