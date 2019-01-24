

subroutine set_vdust(ilevel)
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
     call vdustfine1(ind_grid,ngrid,ilevel)
  end do
  

  


111 format('   Entering dust_diffusion_fine for level ',i2)


end subroutine set_vdust
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine vdustfine1(ind_grid,ncache,ilevel)
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
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim),save::vloc=0.0d0
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim),save::gloc=0.0d0

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
  integer::  ncycle,icycle
  real(dp):: dt_dustcycle
  logical :: d_cycle_ok
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

        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              do i=1,nbuffer
                 gloc(ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if
  

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


  !call dustdiff_split(uloc,flux,dx,dx,dx,dtnew(ilevel),ncache)
  call cmpvdust(uloc,vloc,dx,dx,dx,dtnew(ilevel),ncache)

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
              do idust=1,ndust
                 v_dust(ind_cell(i),idust,idim)=vloc(i,i3,j3,k3,idust,idim)
              enddo
           end if
     end do
  end do
  end do
  end do
end do


 
end subroutine vdustfine1


!###########################################################
!###########################################################
!###########################################################
!###########################################################

subroutine cmpvdust(uin,vout,dx,dy,dz,dt,ngrid)
  use amr_parameters
  use hydro_parameters
  use hydro_commons
  use units_commons
  use const
  implicit none

  integer ::ngrid
  real(dp)::dx, dy, dz, dt

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar+3)::uin
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndust,1:ndim)::vout
  ! Primitive variables
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::qin 
  real(dp),dimension(1:nvector,iu1:iu2+1,ju1:ju2+1,ku1:ku2+1,1:3)::bf
    real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gravin 

  ! declare local variables
  integer ::i, j, k, l
  integer ::ilo,ihi,jlo,jhi,klo,khi
  integer ::idust,idim
  real(dp) :: d,u,v,w,e_mag,e_kin, sum_dust, ening, pressure, cs,A,B,C,wnorm,vmax
  real(dp)::  dAy, dAz,dBx,dBz,dCx,dCy
  real(dp):: fx,fy,fz
  real(dp),dimension(1:ndim) :: fpress
  real(dp),dimension(1:ndust)  :: t_stop
  real(dp)  ::pi,tstop_tot,t_stop_floor,dens_floor
  real(dp), dimension(1:ndust) ::d_grain,l_grain,isnot_charged
  real(dp),dimension(1:ndim):: ew
  ilo=MIN(1,iu1+1); ihi=MAX(1,iu2-1)
  jlo=MIN(1,ju1+1); jhi=MAX(1,ju2-1)
  klo=MIN(1,ku1+1); khi=MAX(1,ku2-1)
  !vmax=vdust_max/scale_v
  pi =3.14159265358979323846_dp
  if(mrn.eqv..true.) then
 isnot_charged=0.0d0    
     call size_dust(l_grain)
     do idust=1,ndust
       l_grain(idust) = l_grain(idust)/scale_l
       d_grain(idust)=grain_dens(idust)/scale_d
       isnot_charged(idust)=0.0d0
    end do
  else
     do idust=1,ndust
       d_grain(idust)=grain_dens(idust)/scale_d
       l_grain(idust)=grain_size(idust)/scale_l
       isnot_charged(idust)=0.0d0

    end do
 endif
 call ctoprimdust(uin,qin,bf,gravin,dt,ngrid)

  do k = klo, khi
     do j = jlo, jhi
        do i = ilo, ihi
           do l = 1, ngrid
                 d= max(qin(l,i,j,k,1),smallr)
                 tstop_tot=0.0d0
                 t_stop=0.0d0
                 cs =sqrt(gamma*qin(l,i,j,k,5)/d)
                 do idust = 1,ndust
                    t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/(d- uin(l,i,j,k,firstindex_ndust+idust))
                    if(K_drag)  t_stop(idust) = uin(l,i,j,k,firstindex_ndust+idust)/K_dust(idust)
                    if(dust_barr) t_stop (idust)= 0.1_dp
                    tstop_tot= tstop_tot-t_stop(idust)*qin(l,i,j,k,firstindex_ndust+idust)
                 end do
                 !magnetic field and required derivatives to get rotB
                 dAy=(0.5d0*(uin(l,i,j+1,k,6)+uin(l,i,j+1,k,nvar+1))-0.5d0*(uin(l,i,j-1,k,6)+uin(l,i,j-1,k,nvar+1)))*0.5d0/dy
                 dAz=(0.5d0*(uin(l,i,j,k+1,6)+uin(l,i,j,k+1,nvar+1))-0.5d0*(uin(l,i,j,k-1,6)+uin(l,i,j,k-1,nvar+1)))*0.5d0/dz
                 A=0.5d0*(uin(l,i,j,k,6)+uin(l,i,j,k,nvar+1))
                 
                 dBx=(0.5d0*(uin(l,i+1,j,k,7)+uin(l,i+1,j,k,nvar+2))-0.5d0*(uin(l,i-1,j,k,7)+uin(l,i-1,j,k,nvar+2)))*0.5d0/dx
                 dBz=(0.5d0*(uin(l,i,j,k+1,7)+uin(l,i,j,k+1,nvar+2))-0.5d0*(uin(l,i,j,k-1,7)+uin(l,i,j,k-1,nvar+2)))*0.5d0/dz
                 B=0.5d0*(uin(l,i,j,k,7)+uin(l,i,j,k,nvar+2))
                 
                 dCy=(0.5d0*(uin(l,i,j+1,k,8)+uin(l,i,j+1,k,nvar+3))-0.5d0*(uin(l,i,j-1,k,8)+uin(l,i,j-1,k,nvar+3)))*0.5d0/dy
                 dCx=(0.5d0*(uin(l,i+1,j,k,8)+uin(l,i+1,j,k,nvar+3))-0.5d0*(uin(l,i-1,j,k,8)+uin(l,i-1,j,k,nvar+3)))*0.5d0/dx
                 C=0.5d0*(uin(l,i,j,k,8)+uin(l,i,j,k,nvar+3))
                 
                 !pressure force
                 fpress(1)=-(qin(l,i+1,j,k,5)-qin(l,i-1,j,k,5))*0.5d0/dx/(d*(1.0d0-sum_dust))
                 fpress(2)=-(qin(l,i,j+1,k,5)-qin(l,i,j-1,k,5))*0.5d0/dy/(d*(1.0d0-sum_dust))
                 fpress(3)=-(qin(l,i,j,k+1,5)-qin(l,i,j,k-1,5))*0.5d0/dz/(d*(1.0d0-sum_dust))
                 !magnetic force
                 fx=((dAz-dCx)*C-(dBx-dAy)*A)/(d*(1.0d0-sum_dust))
                 fy=((dBx-dAy)*A-(dCy-dBz)*C)/(d*(1.0d0-sum_dust))
                 fz=((dCy-dBz)*B-(dAz-dCx)*A)/(d*(1.0d0-sum_dust))
                 
                 do idust = 1,ndust
                    t_stop(idust) = t_stop(idust)+tstop_tot
                    vout(l,i,j,k,idust,1)= t_stop(idust)*(1.0d0-sum_dust)*(-fpress(1)-fx*isnot_charged(idust))
                    vout(l,i,j,k,idust,2)= t_stop(idust)*(1.0d0-sum_dust)*(-fpress(2)-fy*isnot_charged(idust))
                    vout(l,i,j,k,idust,3)= t_stop(idust)*(1.0d0-sum_dust)*(-fpress(3)-fz*isnot_charged(idust))
   
#if NDIM==1
                 wnorm= sqrt(vout(l,i,j,k,idust,1)**2.0)
#endif
#if NDIM==2
                 wnorm= sqrt(vout(l,i,j,k,idust,1)**2.0+vout(l,i,j,k,idust,2)**2.0)
                 
#endif
#if NDIM==3
                 wnorm= sqrt(vout(l,i,j,k,idust,1)**2.0+vout(l,i,j,k,idust,2)**2.0+vout(l,i,j,k,idust,3)**2.0)

#endif

           vmax = sqrt(qin(l,i,j,k,2)**2.0 + qin(l,i,j,k,3)**2.0  + qin(l,i,j,k,4)**2.0 )
           if (reduce_wdust) then   
                 do idim=1,ndim
                    if(wnorm>vmax) ew(idim)= vout(l,i,j,k,idust,idim)/wnorm        
                    if(wnorm>vmax) vout(l,i,j,k,idust,idim)=  ew(idim)*vmax
                    !print *, wnorm, vmax, ew(idim),cs, t_stop(idust),fgas(i,idim)
              end do

           
           end if
                  end do
           end do
        end do
     end do
  end do


end subroutine cmpvdust

subroutine ctoprimdust(uin,q,bf,gravin,dt,ngrid)
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
 
end subroutine ctoprimdust
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!!$subroutine set_vdust(ilevel)
!!$  use amr_commons
!!$  use hydro_commons
!!$  use units_commons
!!$  use cloud_module
!!$  use cooling_module,ONLY:kB,mH
!!$  use radiation_parameters
!!$
!!$  implicit none
!!$  integer::ilevel
!!$  integer::i,j,k,ivar,irad,ind,iskip,nx_loc,ind_cell1,idust
!!$  integer::ncache,igrid,ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
!!$  integer,dimension(1:3,1:2,1:8)::iii,jjj
!!$  real(dp)::scale,dx,dx_loc,d,u,v,w,eold,A,B,C,pressure
!!$
!!$  integer ,dimension(1:nvector),save::ind_grid,ind_cell
!!$  integer ,dimension(1:nvector,0:twondim),save::igridn
!!$  integer ,dimension(1:nvector,1:ndim),save::ind_left,ind_right
!!$  real(dp),dimension(1:nvector,1:ndim),save::dx_g,dx_d
!!$  real(dp)::usquare,emag,erad_loc,ekin,eps,sum_dust,enint
!!$  real(dp)::e_mag,e_kin,e_cons,e_prim,e_trunc,div,fact,e_r
!!$  real(dp)::Pgdivu,u_square,d_loc,Tp_loc,Tr_loc,cal_Teg
!!$  real(dp),dimension(1:nvector,1:ndim),save::Pleft,Pright
!!$  real(dp),dimension(1:nvector,1:ndim),save:: Aleft,Aright,Bleft,Bright,Cleft,Cright
!!$
!!$  real(dp),dimension(1:nvector,1:ndim)       :: gradP
!!$  real(dp),dimension(1:nvector,1:ndim)       :: rotB,BgradB
!!$  real(dp),dimension(1:nvector,1:ndim)       :: fgas,fmag
!!$  real(dp),dimension(1:nvector,1:ndim,1:ndust)       :: fdust
!!$  
!!$  real(dp),dimension(1:ndust)  :: t_stop
!!$  real(dp)  ::cs,pi,tstop_tot,t_stop_floor,dens_floor
!!$  real(dp), dimension(1:ndust) ::d_grain,l_grain
!!$  real(dp) :: dd,ee,cmp_Cv_eos,d0,r0
!!$  integer  :: ht
!!$  real(dp):: epsilon_0, dt_dust,wnorm,vmax
!!$  real(dp),dimension(1:ndim):: ew
!!$  real(dp),dimension(1:ndust):: dustMRN
!!$  epsilon_0 = dust_ratio(1)
!!$  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
!!$  if(numbtot(1,ilevel)==0)return
!!$  if(verbose)write(*,111)ilevel
!!$  nx_loc=icoarse_max-icoarse_min+1
!!$  scale=boxlen/dble(nx_loc)
!!$  dx=0.5d0**ilevel
!!$  dx_loc=dx*scale
!!$  t_stop_floor=0.0d0
!!$  sum_dust=0.0d0
!!$  Pleft=0.0; Pright=0.0
!!$  t_stop=0.0d0
!!$  sum_dust=0.0d0
!!$  pi =3.14159265358979323846_dp
!!$  dens_floor=0.0d0
!!$#if NDUST>0
!!$     do idust =1,ndust
!!$        dustMRN(idust) = dust_ratio(idust)/(1.0d0+dust_ratio(idust))
!!$     end do     
!!$     if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
!!$     do idust =1,ndust
!!$           sum_dust = sum_dust + dustMRN(idust)
!!$        end do   
!!$#endif   
!!$  r0=(alpha_dense_core*2.*6.67d-8*mass_c*scale_m*mu_gas*mH/(5.*kB*Tr_floor*(1.0d0-sum_dust)))/scale_l
!!$  d0 = 3.0d0*mass_c/(4.0d0*pi*r0**3.)
!!$  dens_floor=d0
!!$
!!$
!!$  vmax=vdust_max/scale_v
!!$  if(mrn.eqv..true.) then
!!$     call size_dust(l_grain)
!!$     do idust=1,ndust
!!$       l_grain(idust) = l_grain(idust)/scale_l
!!$       d_grain(idust)=grain_dens(idust)/scale_d
!!$    end do
!!$  else
!!$     do idust=1,ndust
!!$       d_grain(idust)=grain_dens(idust)/scale_d
!!$       l_grain(idust)=grain_size(idust)/scale_l
!!$    end do
!!$ endif
!!$  
!!$  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
!!$  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
!!$  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
!!$  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
!!$  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
!!$  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)
!!$  ncache=active(ilevel)%ngrid
!!$  do igrid=1,ncache,nvector
!!$   
!!$     ! Gather nvector grids
!!$     ngrid=MIN(nvector,ncache-igrid+1)
!!$     do i=1,ngrid
!!$        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
!!$     end do
!!$     
!!$     ! Gather neighboring grids
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
!!$     
!!$     ! Loop over cells
!!$     do ind=1,twotondim
!!$        
!!$        ! Compute central cell index
!!$        iskip=ncoarse+(ind-1)*ngridmax
!!$        do i=1,ngrid
!!$           ind_cell(i)=iskip+ind_grid(i)
!!$        end do
!!$        
!!$        ! Gather all neighboring velocities
!!$        do idim=1,ndim
!!$           id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
!!$           ih1=ncoarse+(id1-1)*ngridmax
!!$           do i=1,ngrid
!!$              if(igridn(i,ig1)>0)then
!!$                 dx_g(i,idim)=dx_loc              
!!$              if(energy_fix)then
!!$                 eold=uold(igridn(i,ig1)+ih1,nvar)
!!$              else
!!$                 ! Gather left thermal energy
!!$                 d=max(uold(igridn(i,ig1)+ih1,1),smallr)
!!$                 u=0.0d0; v=0.0d0; w=0.0d0
!!$                 if(ndim>0)u=uold(igridn(i,ig1)+ih1,2)/d
!!$                 if(ndim>1)v=uold(igridn(i,ig1)+ih1,3)/d
!!$                 if(ndim>2)w=uold(igridn(i,ig1)+ih1,4)/d
!!$                 A=0.5d0*(uold(igridn(i,ig1)+ih1,6)+uold(igridn(i,ig1)+ih1,nvar+1))
!!$                 B=0.5d0*(uold(igridn(i,ig1)+ih1,7)+uold(igridn(i,ig1)+ih1,nvar+2))
!!$                 C=0.5d0*(uold(igridn(i,ig1)+ih1,8)+uold(igridn(i,ig1)+ih1,nvar+3))
!!$                 eold=uold(igridn(i,ig1)+ih1,5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)
!!$#if NENER>0
!!$                 do irad=1,nener
!!$                    eold=eold-uold(igridn(i,ig1)+ih1,8+irad)
!!$                 end do
!!$#endif
!!$              endif
!!$              sum_dust=0.0d0
!!$              do idust = 1, Ndust
!!$                 sum_dust=sum_dust+uold(igridn(i,ig1)+ih1,firstindex_ndust+idust)/d
!!$              end do
!!$              call pressure_eos((1.0_dp-sum_dust)*d,eold,Pleft(i,idim))
!!$              Aleft(i,idim)=A
!!$              Bleft(i,idim)=B
!!$              Cleft(i,idim)=C
!!$              
!!$           else
!!$              dx_g(i,idim)=dx_loc*1.5d0         
!!$              if(energy_fix)then
!!$               eold=uold(ind_left(i,idim),nvar)
!!$            else
!!$              ! Gather left thermal energy
!!$               d=max(uold(ind_left(i,idim),1),smallr)
!!$               u=0.0; v=0.0; w=0.0
!!$               if(ndim>0)u=uold(ind_left(i,idim),2)/d
!!$               if(ndim>1)v=uold(ind_left(i,idim),3)/d
!!$               if(ndim>2)w=uold(ind_left(i,idim),4)/d
!!$               A=0.5d0*(uold(ind_left(i,idim),6)+uold(ind_left(i,idim),nvar+1))
!!$               B=0.5d0*(uold(ind_left(i,idim),7)+uold(ind_left(i,idim),nvar+2))
!!$               C=0.5d0*(uold(ind_left(i,idim),8)+uold(ind_left(i,idim),nvar+3))
!!$               eold=uold(ind_left(i,idim),5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)
!!$#if NENER>0
!!$               do irad=1,nener
!!$                  eold=eold-uold(ind_left(i,idim),8+irad)
!!$               end do
!!$#endif
!!$            endif
!!$            sum_dust=0.0d0
!!$            do idust = 1, Ndust
!!$               sum_dust=sum_dust+uold(ind_left(i,idim),firstindex_ndust+idust)/d
!!$            end do
!!$            call pressure_eos((1.0_dp-sum_dust)*d,eold,Pleft(i,idim))
!!$            Aleft(i,idim)=A
!!$            Bleft(i,idim)=B
!!$            Cleft(i,idim)=C
!!$         endif
!!$      end do
!!$   end do
!!$   do idim=1,ndim
!!$           id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
!!$           ih2=ncoarse+(id2-1)*ngridmax
!!$           do i=1,ngrid
!!$              if(igridn(i,ig2)>0)then
!!$              dx_d(i,idim)=dx_loc              
!!$              if(energy_fix)then
!!$               eold=uold(igridn(i,ig2)+ih2,nvar)
!!$              else
!!$              ! Gather right thermal energy
!!$              d=max(uold(igridn(i,ig2)+ih2,1),smallr)
!!$              u=0.0; v=0.0; w=0.0
!!$              if(ndim>0)u=uold(igridn(i,ig2)+ih2,2)/d
!!$              if(ndim>1)v=uold(igridn(i,ig2)+ih2,3)/d
!!$              if(ndim>2)w=uold(igridn(i,ig2)+ih2,4)/d
!!$              A=0.5d0*(uold(igridn(i,ig2)+ih2,6)+uold(igridn(i,ig2)+ih2,nvar+1))
!!$              B=0.5d0*(uold(igridn(i,ig2)+ih2,7)+uold(igridn(i,ig2)+ih2,nvar+2))
!!$              C=0.5d0*(uold(igridn(i,ig2)+ih2,8)+uold(igridn(i,ig2)+ih2,nvar+3))
!!$              eold=uold(igridn(i,ig2)+ih2,5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)
!!$#if NENER>0
!!$              do irad=1,nener
!!$                 eold=eold-uold(igridn(i,ig2)+ih2,8+irad)
!!$              end do
!!$#endif
!!$              endif    
!!$              sum_dust=0.0d0
!!$              do idust = 1, Ndust
!!$                 sum_dust=sum_dust+uold(igridn(i,ig2)+ih2,firstindex_ndust+idust)/d
!!$              end do
!!$              call pressure_eos((1.0_dp-sum_dust)*d,eold,Pright(i,idim))
!!$              Aright(i,idim)=A
!!$              Bright(i,idim)=B
!!$              Cright(i,idim)=C
!!$           else
!!$              dx_d(i,idim)=dx_loc*1.5d0
!!$              if(energy_fix)then
!!$              eold=uold(ind_right(i,idim),nvar)
!!$              else
!!$              ! Gather right thermal energy
!!$              d=max(uold(ind_right(i,idim),1),smallr)
!!$              u=0.0d0; v=0.0d0; w=0.0d0
!!$              if(ndim>0)u=uold(ind_right(i,idim),2)/d
!!$              if(ndim>1)v=uold(ind_right(i,idim),3)/d
!!$              if(ndim>2)w=uold(ind_right(i,idim),4)/d
!!$              A=0.5d0*(uold(ind_right(i,idim),6)+uold(ind_right(i,idim),nvar+1))
!!$              B=0.5d0*(uold(ind_right(i,idim),7)+uold(ind_right(i,idim),nvar+2))
!!$              C=0.5d0*(uold(ind_right(i,idim),8)+uold(ind_right(i,idim),nvar+3))
!!$              eold=uold(ind_right(i,idim),5)-0.5d0*d*(u**2+v**2+w**2)-0.5d0*(A**2+B**2+C**2)
!!$#if NENER>0
!!$              do irad=1,nener
!!$                 eold=eold-uold(ind_right(i,idim),8+irad)
!!$              end do
!!$#endif
!!$              endif                  
!!$              sum_dust=0.0d0
!!$              do idust = 1, Ndust
!!$                 sum_dust=sum_dust+uold(ind_right(i,idim),firstindex_ndust+idust)/d
!!$              end do
!!$              call pressure_eos((1.0_dp-sum_dust)*d,eold,Pright(i,idim))
!!$              Aright(i,idim)=A
!!$              Bright(i,idim)=B
!!$              Cright(i,idim)=C
!!$           endif
!!$        end do
!!$     end do
!!$     do idim=1,ndim
!!$           do i=1,ngrid
!!$              gradP(i,idim) = (Pright(i,idim)-Pleft(i,idim))/(dx_g(i,idim)+dx_d(i,idim))
!!$           end do
!!$        end do
!!$        do i=1,ngrid
!!$           d=max(uold(ind_cell(i),1),smallr)
!!$           sum_dust=0.0_dp
!!$           do idust = 1, ndust
!!$              sum_dust = sum_dust + uold(ind_cell(i),firstindex_ndust+idust)/d
!!$              fdust(i,idim,idust)=0.0d0!
!!$           end do
!!$           A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
!!$           B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
!!$           C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
!!$           rotB(i,1)=(Cright(i,2)-Cleft(i,2))/(dx_g(i,2)+dx_d(i,2))-(Bright(i,3)-Bleft(i,3))/(dx_g(i,3)+dx_d(i,3))
!!$           rotB(i,2)=(Aright(i,3)-Aleft(i,3))/(dx_g(i,3)+dx_d(i,3))-(Cright(i,1)-Cleft(i,1))/(dx_g(i,1)+dx_d(i,1))
!!$           rotB(i,3)=(Bright(i,1)-Bleft(i,1))/(dx_g(i,1)+dx_d(i,1))-(Aright(i,2)-Aleft(i,2))/(dx_g(i,2)+dx_d(i,2))
!!$           fmag(i,1)=(rotB(i,2)*C-rotB(i,3)*A)/d/(1.0-sum_dust)
!!$           fmag(i,2)=(rotB(i,3)*A-rotB(i,1)*C)/d/(1.0-sum_dust)
!!$           fmag(i,3)=(rotB(i,1)*B-rotB(i,2)*A)/d/(1.0-sum_dust)
!!$
!!$           !non drag gas forces  
!!$           do idim=1,ndim
!!$              fgas(i,idim)=-gradP(i,idim)/d/(1.0-sum_dust)!+fmag(i,idim)
!!$              !print *, fmag(i,idim),-gradP(i,idim)/d/(1.0-sum_dust)
!!$           end do
!!$        end do
!!$        do i=1,ngrid
!!$           if (energy_fix) then
!!$              d=max(uold(ind_cell(i),1),smallr)
!!$              enint=uold(ind_cell(i),nvar)
!!$           else
!!$              u=0.0d0; v=0.0d0; w=0.0d0
!!$
!!$              d=max(uold(ind_cell(i),1),smallr)
!!$              if(ndim>0)u=uold(ind_cell(i),2)/d
!!$              if(ndim>1)v=uold(ind_cell(i),3)/d
!!$              if(ndim>2)w=uold(ind_cell(i),4)/d
!!$              e_mag= 0.0_dp
!!$#ifdef SOLVERmhd                           
!!$              A=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
!!$              B=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
!!$              C=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
!!$              e_mag=0.5d0*(A**2+B**2+C**2)
!!$#endif
!!$              e_kin=0.5d0*d*(u**2+v**2+w**2)
!!$#if NENER>0
!!$              do irad=1,nener
!!$                 e_mag=e_mag+uold(ind_cell(i),8+irad)
!!$              end do
!!$#endif
!!$              if(dust_barr) then
!!$                 enint=uold(ind_cell(i),5)
!!$              else   
!!$                 enint=uold(ind_cell(i),5)-e_kin- e_mag
!!$              endif
!!$           endif
!!$           sum_dust=0.0_dp
!!$           do idust = 1, ndust
!!$              sum_dust = sum_dust + uold(ind_cell(i),firstindex_ndust+idust)/d
!!$           end do
!!$           call pressure_eos  ((1.0_dp-sum_dust)*d,enint,pressure)
!!$           call soundspeed_eos((1.0_dp-sum_dust)*d,enint, cs)
!!$           if(dust_barr)  cs = 1.0_dp
!!$           if(dust_barr) pressure = (1.0_dp-sum_dust)*d*cs*cs
!!$
!!$            sum_dust=0.0d0
!!$            do idust = 1, ndust
!!$               sum_dust=sum_dust+uold(ind_cell(i),firstindex_ndust+idust)/d
!!$            end do
!!$            tstop_tot=0.0d0
!!$            t_stop=0.0d0
!!$            do idust = 1,ndust
!!$               t_stop(idust) =  d_grain(idust)*l_grain(idust)*SQRT(pi*gamma/8.0_dp)/cs/(d-uold(ind_cell(i),firstindex_ndust+idust))
!!$               if(K_drag)  t_stop(idust) = uold(ind_cell(i),firstindex_ndust+idust)/K_dust(idust)
!!$               if(dust_barr) t_stop (idust)= 0.1_dp
!!$               tstop_tot= tstop_tot-t_stop(idust)*(uold(ind_cell(i),firstindex_ndust+idust)/d)
!!$            end do
!!$            do idust = 1,ndust
!!$               t_stop(idust) = t_stop(idust)+tstop_tot
!!$               do idim=1,ndim
!!$                  v_dust(ind_cell(i),idust,idim)= t_stop(idust)*(1.0d0-sum_dust)*(fdust(i,idim,idust)-fgas(i,idim))
!!$                 ! print *,v_dust(ind_cell(i),idust,idim), idim
!!$		end do	  
!!$              if(reduce_wdust) then
!!$#if NDIM==1
!!$                 wnorm= sqrt(v_dust(ind_cell(i),idust,1)**2.0)
!!$#endif
!!$#if NDIM==2
!!$                 wnorm= sqrt(v_dust(ind_cell(i),idust,1)**2.0+v_dust(ind_cell(i),idust,2)**2.0)
!!$                 
!!$#endif
!!$#if NDIM==3
!!$                 wnorm= sqrt(v_dust(ind_cell(i),idust,1)**2.0+v_dust(ind_cell(i),idust,2)**2.0+v_dust(ind_cell(i),idust,3)**2.0)
!!$#endif
!!$
!!$                 
!!$                 do idim=1,ndim
!!$                    if(wnorm>vmax) ew(idim)= v_dust(ind_cell(i),idust,idim)/wnorm        
!!$                    if(wnorm>vmax) v_dust(ind_cell(i),idust,idim)=  ew(idim)*vmax
!!$                    !print *, wnorm, vmax, ew(idim),cs, t_stop(idust),fgas(i,idim)
!!$              end do
!!$
!!$           
!!$            end if
!!$         end do
!!$         end do
!!$      end do
!!$enddo
!!$
!!$111 format('   Entering set_vdust for level ',i2)
!!$
!!$end subroutine set_vdust


