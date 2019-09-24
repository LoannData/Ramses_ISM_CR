!############################################################
!############################################################
!############################################################
!############################################################
subroutine rt_boundana(x,u,dx,ibound,ncell,ilevel)
  use amr_parameters!, ONLY: dp,ndim,nvector
  use rt_parameters!, ONLY: nrtvar, rt_boundary_var
  use amr_commons
  use rt_hydro_commons
  use hydro_commons
  ! CC 02/17
  use cloud_module
  implicit none
  integer ::ibound                          ! Index of boundary region
  integer ::ncell                           ! Number of active cells
  real(dp)::dx                              ! Cell size
  real(dp),dimension(1:nvector,1:nrtvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x   ! Cell center position.
  integer,dimension(1:nvector)::cell_index,cell_index2
  integer,dimension(1:nvector)::cell_levl,cell_levl2
  real(dp),dimension(1:nvector,1:ndim)::y
  real(dp),dimension(1:nvector,1:ndim)::xtcell1,xtcell2,xtcell3
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): N, U(i,2:ndim+1): Fx,Fy,Fz.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i,idi,ilevel
  real(dp)::diff,w1,w2

  dx = 0.5**ilevel
  
  ! CC 01/2017

  ! Periodicite et shear sur la boundary 3 (en bas)
!  if(ibound.eq.3) then
  if(ibound.eq.1) then
     do i=1,ncell
        y(i,1)=(1.0/boxlen)*mod(x(i,1)+Vshear*t,boxlen)
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1
        end if
        do idi=2,ndim
           y(i,idi)=(1.0/boxlen)*mod(x(i,idi),boxlen)
           if (y(i,idi).lt.0) then
              y(i,idi)=y(i,idi)+1
           end if
        end do
     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)
    
     ! interpolation : no asmr assumed here
     do i=1,ncell
        diff=y(i,1)-xtcell1(i,1)
        if(abs(diff)<dx/1000.) then
           diff = 0d0
        else
           diff=diff/abs(diff)
        end if 
        xtcell2(i,1)=mod(xtcell1(i,1)+dx*diff,1.d0)
        if(xtcell2(i,1)<0.) xtcell2(i,1)=xtcell2(i,1)+1.
        xtcell2(i,2)=xtcell1(i,2)
        xtcell2(i,3)=xtcell1(i,3)
     end do

     call get_cell_index4(cell_index2,cell_levl2,xtcell2,xtcell3,ilevel,ncell)

     do i=1,ncell
        if(diff.eq.0d0) then
           w1 = 1.d0
        else if(abs(xtcell2(i,1)-y(i,1)) .le. dx) then
           w1 =  abs(xtcell2(i,1)-y(i,1))  / dx
        else
           w1 =  ( 1.-abs(xtcell2(i,1)-y(i,1)) ) / dx
        end if
     
        w2=1.-w1

        do ivar=1,nrtvar
              u(i,ivar)=w1*rtuold(cell_index(i),ivar)+w2*rtuold(cell_index2(i),ivar)

              if(isnan(u(i,ivar))) then 
                 write(*,*) 'w1,w2,rtuold(cell_index(i),ivar),rtuold(cell_index2(i),ivar)'
                 write(*,*) w1,w2,rtuold(cell_index(i),ivar),rtuold(cell_index2(i),ivar)
              endif

        end do
     end do 
  end if

  ! Periodicite et shear sur la boundary 4 (en haut)
!  if(ibound.eq.4) then
  if(ibound.eq.2) then
     do i=1,ncell
        y(i,1)=(1.0/boxlen)*mod(x(i,1)-Vshear*t,boxlen)
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1
        end if
        do idi=2,ndim
           y(i,idi)=(1.0/boxlen)*mod(x(i,idi),boxlen)
           if (y(i,idi).lt.0) then
              y(i,idi)=y(i,idi)+1
           end if
        end do
     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)
     
     ! interpolation : no asmr assumed here
     do i=1,ncell
        diff=y(i,1)-xtcell1(i,1)
        if(abs(diff)<dx/1000.) then
            diff = 0d0
        else
            diff=diff/abs(diff)
        end if
        xtcell2(i,1)=mod(xtcell1(i,1)+dx*diff,1.d0)
        if(xtcell2(i,1)<0.) xtcell2(i,1)=xtcell2(i,1)+1.
        xtcell2(i,2)=xtcell1(i,2)
        xtcell2(i,3)=xtcell1(i,3)
     end do

     call get_cell_index4(cell_index2,cell_levl2,xtcell2,xtcell3,ilevel,ncell)

     do i=1,ncell
        if(diff.eq.0d0) then
           w1 = 1.d0
        else if(abs(xtcell2(i,1)-y(i,1)) .le. dx) then
           w1 =   abs(xtcell2(i,1)-y(i,1)) / dx
        else
           w1 =  ( 1.-abs(xtcell2(i,1)-y(i,1)) ) / dx
        end if
     
     w2=1.-w1

        do ivar=1,nrtvar
              u(i,ivar)=w1*rtuold(cell_index(i),ivar)+w2*rtuold(cell_index2(i),ivar)

              if(isnan(u(i,ivar))) then 
                 write(*,*) 'w1,w2,rtuold(cell_index(i),ivar),rtuold(cell_index2(i),ivar)'
                 write(*,*) w1,w2,rtuold(cell_index(i),ivar),rtuold(cell_index2(i),ivar)
              endif

        end do

     end do 

  end if

  ! Gradient nul sur les boundaries 5 et 6 (bas et haut)
  ! Les valeurs des cellules fantomes sont egales a celles de la cellule du
  ! domaine juste en dessous
!  if((ibound.eq.5).or.(ibound.eq.6)) then
  if((ibound.eq.3).or.(ibound.eq.4)) then
     do i=1,ncell
!        if(ibound.eq.6) then
        if(ibound.eq.4) then
           y(i,3)=1-0.5*(0.5)**ilevel
        else
           y(i,3)=0.5*(0.5)**ilevel
        end if
        y(i,2)=(1.0/boxlen)*mod(x(i,2),boxlen)
        if ((x(i,2).lt.boxlen) .and. (x(i,2).gt.0.)) then ! On cherche la cellule du domaine juste en dessous, modulo la périodicité  
           y(i,1)=(1.0/boxlen)*mod(x(i,1),boxlen)
        else if (x(i,2).lt.0.) then                       ! On decale cette fois avec le shear (vers la droite)
           y(i,1)=(1.0/boxlen)*mod(x(i,1)+Vshear*t,boxlen)
        else                                              ! On decale vers la gauche
           y(i,1)=(1.0/boxlen)*mod(x(i,1)-Vshear*t,boxlen)
        end if

        if (y(i,2).lt.0) then
           y(i,2)=y(i,2)+1
        end if
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        do ivar=1,nrtvar
              u(i,ivar)=rtuold(cell_index(i),ivar)
        end do
     end do
  end if

end subroutine rt_boundana
