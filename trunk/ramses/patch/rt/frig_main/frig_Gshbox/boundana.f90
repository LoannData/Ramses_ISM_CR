!############################################################
!############################################################
!############################################################
!############################################################
subroutine boundana(x,u,dx,ibound,ncell,ilevel)
  use amr_parameters!, ONLY: dp,ndim,nvector
  use hydro_parameters!, ONLY: nvar,boundary_var
  use amr_commons
  use hydro_commons
  use cloud_module
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  integer,dimension(1:nvector)::cell_index,cell_index2
  integer,dimension(1:nvector)::cell_levl,cell_levl2
  real(dp),dimension(1:nvector,1:ndim)::y
  real(dp),dimension(1:nvector,1:ndim)::xtcell1,xtcell2,xtcell3
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i,idi,ilevel

  real(dp)::dx_loc,emag,diff,w1,w2

  dx_loc = 0.5**ilevel * boxlen
  dx = 0.5**ilevel

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! CC 01/12/2016

  ! y sont les coordonnees de la cellule du domaine correspondant a la cellule
  ! fantome, obtenues par periodicite et par decalage vers la gauche (boundary
  ! #2) ou la droite (boundary #1) du au shear

  ! Periodicite et shear sur la boundary 3 (en bas)
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
    
     ! interpolation : no amr assumed here
     do i=1,ncell
        diff=y(i,1)-xtcell1(i,1)
        if(abs(diff)<dx/1000.) then
            diff=0d0
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
           w1 =  ( 1.-abs(xtcell2(i,1)-y(i,1))) / dx
        end if
     
        w2=1.-w1
       
        u(i,1) = w1*uold(cell_index(i),1) + w2*uold(cell_index2(i),1) 
        u(i,2) = w1*uold(cell_index(i),2) + w2*uold(cell_index2(i),2) - u(i,1)*Vshear
        do ivar=3,nvar+3
           u(i,ivar)= w1*uold(cell_index(i),ivar) + w2*uold(cell_index2(i),ivar)
        end do
        ! Correct kinetic energy
        u(i,ndim+2) = u(i,ndim+2) - w1*0.5*( (uold(cell_index(i),2))**2 + (uold(cell_index(i),3))**2 + (uold(cell_index(i),4))**2 )/uold(cell_index(i),1)
        u(i,ndim+2) = u(i,ndim+2) - w2*0.5*( (uold(cell_index2(i),2))**2 + (uold(cell_index2(i),3))**2 + (uold(cell_index2(i),4))**2 )/uold(cell_index2(i),1)
        u(i,ndim+2) = u(i,ndim+2) + 0.5*(u(i,2)**2+u(i,3)**2+u(i,4)**2) / u(i,1)

        ! Correct magnetic energy
        u(i,ndim+2) = u(i,ndim+2) - w1*0.125*( (uold(cell_index(i),6)+uold(cell_index(i),nvar+1))**2 + (uold(cell_index(i),7)+uold(cell_index(i),nvar+2))**2 + (uold(cell_index(i),8)+uold(cell_index(i),nvar+3))**2 )
        u(i,ndim+2) = u(i,ndim+2) - w2*0.125*( (uold(cell_index2(i),6)+uold(cell_index2(i),nvar+1))**2 + (uold(cell_index2(i),7)+uold(cell_index2(i),nvar+2))**2 + (uold(cell_index2(i),8)+uold(cell_index2(i),nvar+3))**2 )
        u(i,ndim+2) = u(i,ndim+2) + 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)

     end do 


     !now we must correct for div B. For that purpose we first grab the cells 
     !located in 0.5^levelmin and we propagate its By field
     do i=1,ncell
        y(i,2)=0.5*(0.5)**ilevel
        y(i,1)=(1.0/boxlen)*mod(x(i,1),boxlen)
        y(i,3)=(1.0/boxlen)*mod(x(i,3),boxlen)

        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1
        end if
        if (y(i,3).lt.0) then
           y(i,3)=y(i,3)+1
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        !calculate the magnetic energy
        emag = 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !remove it from total energy
        u(i,ndim+2)=u(i,ndim+2)-emag

        !impose the continuity of Bperp
        u(i,nvar+2)=uold(cell_index(i),7)

        !u(i,6)=uold(cell_index(i),6)
        !u(i,8)=uold(cell_index(i),8)
        !u(i,nvar+1)=uold(cell_index(i),nvar+1)
        !u(i,nvar+3)=uold(cell_index(i),nvar+3)

        !now impose div B within the cell
        u(i,7)=u(i,nvar+2)-u(i,6)+u(i,nvar+1)-u(i,8)+u(i,nvar+3)
     end do 


     !now depending whether the ghost cell is located at -0.5^levelmin / 2 or at -0.5^levelmin / 2 * 3 
     !(in this latter case) we must grab the neigbour of the one located at -0.5^levelmin / 2 to get div B =0
     do i=1,ncell

        y(i,1)=(1.0/boxlen)*mod(x(i,1)+Vshear*t,boxlen)
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1
        end if

        !shift the cells from dx_loc
        y(i,2)=(1.0/boxlen)*mod(x(i,2)+dx_loc,boxlen)
        if (y(i,2).lt.0) then
           y(i,2)=y(i,2)+1
        end if

        y(i,3)=(1.0/boxlen)*mod(x(i,3),boxlen)
        if (y(i,3).lt.0) then
           y(i,3)=y(i,3)+1
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        !check whether the correction to the normal field  must be applied
        if(x(i,2) .lt. -dx_loc) then
           !impose the continuity of Bperp by adding the difference of the tangential field
!           u(i,nvar+2)=u(i,nvar+2)-uold(cell_index(i),6)+uold(cell_index(i),nvar+1)-uold(cell_index(i),8)+uold(cell_index(i),nvar+3)
           !apply the correction to the other component
!           u(i,7)=u(i,7)-uold(cell_index(i),6)+uold(cell_index(i),nvar+1)-uold(cell_index(i),8)+uold(cell_index(i),nvar+3)

           u(i,nvar+2)=u(i,nvar+2)-u(i,6)+u(i,nvar+1)-u(i,8)+u(i,nvar+3)
           u(i,7)=u(i,7)-u(i,6)+u(i,nvar+1)-u(i,8)+u(i,nvar+3)
        endif
        !calculate the magnetic energy
        emag = 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !add it back
        u(i,ndim+2)=u(i,ndim+2)+emag
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

     ! interpolation : no amr assumed here
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
           w1=1.d0
        else if(abs(xtcell2(i,1)-y(i,1)) .le. dx) then
           w1 =  abs(xtcell2(i,1)-y(i,1))  /  dx
        else
           w1 =  ( 1.-abs(xtcell2(i,1)-y(i,1)) ) /  dx
        end if
     
        w2=1.-w1
 
        u(i,1) = w1*uold(cell_index(i),1) + w2*uold(cell_index2(i),1) 
        u(i,2) = w1*uold(cell_index(i),2) + w2*uold(cell_index2(i),2) + u(i,1)*Vshear
        do ivar=3,nvar+3
           u(i,ivar)= w1*uold(cell_index(i),ivar) + w2*uold(cell_index2(i),ivar)
        end do
        ! Correct kinetic energy
        u(i,ndim+2) = u(i,ndim+2) - w1*0.5*( (uold(cell_index(i),2))**2 + (uold(cell_index(i),3))**2 + (uold(cell_index(i),4))**2 )/uold(cell_index(i),1)
        u(i,ndim+2) = u(i,ndim+2) - w2*0.5*( (uold(cell_index2(i),2))**2 + (uold(cell_index2(i),3))**2 + (uold(cell_index2(i),4))**2 )/uold(cell_index2(i),1)
        u(i,ndim+2) = u(i,ndim+2) + 0.5*(u(i,2)**2+u(i,3)**2+u(i,4)**2) / u(i,1)

        ! Correct magnetic energy
        u(i,ndim+2) = u(i,ndim+2) - w1*0.125*( (uold(cell_index(i),6)+uold(cell_index(i),nvar+1))**2 + (uold(cell_index(i),7)+uold(cell_index(i),nvar+2))**2 + (uold(cell_index(i),8)+uold(cell_index(i),nvar+3))**2 )
        u(i,ndim+2) = u(i,ndim+2) - w2*0.125*( (uold(cell_index2(i),6)+uold(cell_index2(i),nvar+1))**2 + (uold(cell_index2(i),7)+uold(cell_index2(i),nvar+2))**2 + (uold(cell_index2(i),8)+uold(cell_index2(i),nvar+3))**2 )
        u(i,ndim+2) = u(i,ndim+2) + 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)

     end do 


     !now we must correct for div B. For that purpose we first grap the cells 
     !located in 0.5^levelmin and we propagate its By field
     do i=1,ncell
        y(i,2)=1.-0.5*(0.5)**ilevel
        y(i,1)=(1.0/boxlen)*mod(x(i,1),boxlen)
        y(i,3)=(1.0/boxlen)*mod(x(i,3),boxlen)

        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1
        end if
        if (y(i,3).lt.0) then
           y(i,3)=y(i,3)+1
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        !calculate the magnetic energy
        emag = 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !remove it from total energy
        u(i,ndim+2)=u(i,ndim+2)-emag

        !impose the continuity of Bperp
        u(i,7)=uold(cell_index(i),nvar+2)
        !now impose div B within the cell
        u(i,nvar+2)=u(i,7)+u(i,6)-u(i,nvar+1)+u(i,8)-u(i,nvar+3)
     end do 


     !now depending whether the ghost cell is located at -0.5^levelmin / 2 or at -0.5^levelmin / 2 * 3 
     !(in this latter case) we must grasp the neigbour of the one located at -0.5^levelmin / 2 to get div B =0
     do i=1,ncell

        y(i,1)=(1.0/boxlen)*mod(x(i,1)+Vshear*t,boxlen)
        if (y(i,1).lt.0) then
           y(i,1)=y(i,1)+1
        end if

        !shift the cells from dx_loc
        y(i,2)=(1.0/boxlen)*mod(x(i,2)-dx_loc,boxlen)
        if (y(i,2).lt.0) then
           y(i,2)=y(i,2)+1
        end if

        y(i,3)=(1.0/boxlen)*mod(x(i,3),boxlen)
        if (y(i,3).lt.0) then
           y(i,3)=y(i,3)+1
        end if

     end do

     call get_cell_index4(cell_index,cell_levl,y,xtcell1,ilevel,ncell)

     do i=1,ncell
        !check whether the correction to the normal field  must be applied
        if(x(i,2) .gt. boxlen+dx_loc) then
           !impose the continuity of Bperp by adding the difference of the tangential field
           u(i,7)=u(i,7)+uold(cell_index(i),6)-uold(cell_index(i),nvar+1)+uold(cell_index(i),8)-uold(cell_index(i),nvar+3)
           !apply the correction to the other component
           u(i,nvar+2)=u(i,nvar+2)+uold(cell_index(i),6)-uold(cell_index(i),nvar+1)+uold(cell_index(i),8)-uold(cell_index(i),nvar+3)

        !u(i,7)=u(i,7)+u(i,6)-u(i,nvar+1)+u(i,8)-u(i,nvar+3)
        !u(i,nvar+2)=u(i,nvar+2)+u(i,6)-u(i,nvar+1)+u(i,8)-u(i,nvar+3)

        endif
        !calculate the magnetic energy
        emag = 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !add it back
        u(i,ndim+2)=u(i,ndim+2)+emag
     end do 

  end if

  ! Gradient nul sur les boundaries 3 et 4 (bas et haut)
  ! Les valeurs des cellules fantomes sont egales a celles de la cellule du
  ! domaine juste en dessous

  if(ibound.eq.3) then
     do i=1,ncell
        y(i,3)=0.5*(0.5)**ilevel

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
        u(i,1)=uold(cell_index(i),1)
        do ivar=3,nvar+3
              u(i,ivar)=uold(cell_index(i),ivar)
        end do
        if ((x(i,2).lt.boxlen) .and. (x(i,2).gt.0.)) then
           u(i,2)=uold(cell_index(i),2)
        else if (x(i,2).lt.0.) then
           u(i,2)=uold(cell_index(i),2)-u(i,1)*Vshear
        else 
           u(i,2)=uold(cell_index(i),2)+u(i,1)*Vshear
        end if

        !impose outflowing conditions
        u(i,4) = min(u(i,4),0.)

        !modify kinetic energy
        u(i,ndim+2)=u(i,ndim+2)-0.5*(uold(cell_index(i),2)**2)/u(i,1)+0.5*(u(i,2)**2)/u(i,1)
        u(i,ndim+2)=u(i,ndim+2)-0.5*(uold(cell_index(i),3)**2)/u(i,1)+0.5*(u(i,3)**2)/u(i,1)
     end do

     do i=1,ncell
        !calculate the magnetic energy
        emag = 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !remove it from total energy
        u(i,ndim+2)=u(i,ndim+2)-emag

        !impose the continuity of Bperp
        u(i,nvar+3)=uold(cell_index(i),8)
        !now impose div B within the cell
        u(i,8)=u(i,nvar+3)-u(i,6)+u(i,nvar+1)-u(i,7)+u(i,nvar+2)

        !now depending whether the ghost cell is located at -0.5^levelmin / 2 or at -0.5^levelmin / 2 * 3 
        !(in this latter case) we must propagate the normal magnetic field through another cell

        if(x(i,3) .lt. -dx_loc) then
           !impose the continuity of Bperp by adding the difference of the tangential field
           u(i,nvar+3)=u(i,nvar+3)-uold(cell_index(i),6)+uold(cell_index(i),nvar+1)-uold(cell_index(i),7)+uold(cell_index(i),nvar+2)
           !apply the correction to the other component
           u(i,8)=u(i,8)-uold(cell_index(i),6)+uold(cell_index(i),nvar+1)-uold(cell_index(i),7)+uold(cell_index(i),nvar+2)
        endif
        !calculate the magnetic energy
        emag = 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !add it back
        u(i,ndim+2)=u(i,ndim+2)+emag

     end do 

  end if



  if(ibound.eq.4) then
     do i=1,ncell

        y(i,3)=1.-0.5*(0.5)**ilevel

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
        u(i,1)=uold(cell_index(i),1)
        do ivar=3,nvar+3
              u(i,ivar)=uold(cell_index(i),ivar)
        end do
        if ((x(i,2).lt.boxlen) .and. (x(i,2).gt.0.)) then
           u(i,2)=uold(cell_index(i),2)
        else if (x(i,2).lt.0.) then
           u(i,2)=uold(cell_index(i),2)-u(i,1)*Vshear
        else 
           u(i,2)=uold(cell_index(i),2)+u(i,1)*Vshear
        end if

        !impose outflowing conditions
        u(i,4) = max(u(i,4),0.)

        !modify kinetic energy
        u(i,ndim+2)=u(i,ndim+2)-0.5*(uold(cell_index(i),2)**2)/u(i,1)+0.5*(u(i,2)**2)/u(i,1)
        u(i,ndim+2)=u(i,ndim+2)-0.5*(uold(cell_index(i),3)**2)/u(i,1)+0.5*(u(i,3)**2)/u(i,1)
     end do

     do i=1,ncell
        !calculate the magnetic energy
        emag = 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !remove it from total energy
        u(i,ndim+2)=u(i,ndim+2)-emag

        !impose the continuity of Bperp
        u(i,8)=uold(cell_index(i),nvar+3)
        !now impose div B within the cell
        u(i,nvar+3)=u(i,8)+u(i,6)-u(i,nvar+1)+u(i,7)-u(i,nvar+2)

        !now depending whether the ghost cell is located at 1-0.5^levelmin / 2 or at 1-0.5^levelmin / 2 * 3 
        !(in this latter case) we must propagate the normal magnetic field through another cell

        if(x(i,3) .gt. boxlen-dx_loc) then
           !impose the continuity of Bperp by adding the difference of the tangential field
           u(i,8)=u(i,8)+uold(cell_index(i),6)-uold(cell_index(i),nvar+1)+uold(cell_index(i),7)-uold(cell_index(i),nvar+2)
           !apply the correction to the other component
           u(i,nvar+3)=u(i,nvar+3)+uold(cell_index(i),6)-uold(cell_index(i),nvar+1)+uold(cell_index(i),7)-uold(cell_index(i),nvar+2)
        endif
        !calculate the magnetic energy
        emag = 0.125*((u(i,6)+u(i,nvar+1))**2+(u(i,7)+u(i,nvar+2))**2+(u(i,8)+u(i,nvar+3))**2)
        !add it back
        u(i,ndim+2)=u(i,ndim+2)+emag

     end do 

  end if





end subroutine boundana

