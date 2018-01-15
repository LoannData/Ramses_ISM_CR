!###########################################################
!!! VAL (09/07/2014)
!!! VAL (update Nov 2014): I added dependence on T and the sticking coefficient for H2 formation
!
! This routine is a simple version of chem_step, but for only H2:
!
! Formation:    H + H + G(cst) -> H2 + G(cst)      
!               k1 = k1_0*sqrt(T/100K)*S(T)
!
!               k1_0 = 3*10^-17 [cm^3 s^-1]      (Kaufman+1999, Sternberg+2014)
!               S(T) = 1/(1 + (T/464K)^1.5) sticking coefficient (Le Bourlot+2012, Bron+2014, Meudon PDR code)
!
! Destruction:  H2 + hv -> 2 H 
!               k2 = kph = fshield( N(H2) ) * exp(-tau_{d,1000} ) * kph,0
!               kph,0 = 3.3*10^-11*factG0 [s^-1] (Glover & Mac Low2007, Draine & Bertoldi1996,  Sternberg+2014)
!###########################################################

SUBROUTINE simple_chem(ilevel)
  use amr_commons
  use hydro_commons
  !USE hydro_parameters
  !use chemo_parameters
  !use thermochemistry
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info,isink
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  !if(myid .EQ. 1) write(*,*) '***VAL: Entering SIMPLE_CHEM'
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     !if(verbose)write(*,110)ind_grid,ngrid,ilevel
     if(verbose)write(*,*) '***VAL: Ready to call simple_chemostep1' ! : ind_grid,ngrid,ilevel=',ind_grid,ngrid,ilevel
     call simple_chemostep1(ind_grid,ngrid,ilevel)
  end do

!110 format('   Ready to call chemostep1 with ind_grid,ngrid,ilevel=:',i2,i2,i2)
111 format('   Entering simple_chem for level',i2)

end subroutine simple_chem

!###########################################################                                          
!###########################################################                                          
!###########################################################                                          
!########################################################### 

subroutine simple_chemostep1(ind_grid,ngrid,ilevel) !
  use cooling_module,ONLY: kB,mH
  use amr_commons
  use hydro_commons
  use chemo_parameters
  use thermochemistry  , verbose_thermo=>verbose, mype_thermo=>mype
  implicit none
  integer                      :: ilevel,ngrid
  integer,dimension(1:nvector) :: ind_grid
  !
  ! variables as in cooling.f90
  integer                      :: i,j,ind,ind2,iskip,idim,nleaf
  real(dp)                     :: scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v, dt_ilev
  real(kind=8)                 :: nH=0.0_dp, nH2=0.0, ntot=0.0
  real(kind=8)                 :: k1_0=0.0_dp, k1=0.0_dp, k2=0.0, kph0, G0, fshield = 1.0 
  integer,dimension(1:nvector),save::ind_cell,ind_leaf

  real(kind=8),dimension(1:nvector),save     :: ekin,emag,T2,erad_loc
  integer                      :: neulS=8+nrad+n_extinct, neulP=5
  real(dp)                     :: testf = 1.d-8

  !!! FIRST ATTEMPT TO FORM H2 USING SIMPLE K
  ! As in cooling_fine:
  !
  ! G0 is the UV field defined by Habing and Draine
  G0 = 1.0_dp/1.7_dp
  G0 = G0*p_UV             ! p_UV: parameter of variation of UV
                           ! defined in the namelist
  k1_0 = 3.0d-17           ! cm3 s-1 Formation (Kaufman+1999) 3d-17/sqrt(100)
  kph0 = 3.3d-11*G0 !*testf        ! s-1 Photodissociation (Eq 29 Glover&MacLow2007)

  ! Loop over cells
  !if(myid .EQ. 1) write(*,*) '***VAL: Entering SIMPLE_CHEMOSTEP1, starting loop over cells'
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     ! Conversion factor from user units to cgs units
     call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
!     scale_nHmy=scale_d/(1.4*mp)   !!! VAL: ND chem
 
     dt_ilev = dtnew(ilevel)

!======================================================================
!VAL: I calculate pressure and then temperature (19/11/2014)

     ! Compute pressure                                                                                                                                                                                                                                                                                                      
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),neulP)
     end do
     do i=1,nleaf
        ekin(i)=0.0d0
     end do
     do idim=1,3
        do i=1,nleaf
           ekin(i)=ekin(i)+0.5_dp*uold(ind_leaf(i),idim+1)**2/uold(ind_leaf(i),1)
        end do
     end do
     do i=1,nleaf
        emag(i)=0.0d0
     enddo
     do idim=1,3
        do i=1,nleaf
           emag(i)=emag(i)+0.125_dp*(uold(ind_leaf(i),idim+neulP)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do
     do i=1,nleaf
        erad_loc(i)=0.0d0
     end do
     do j=1,nrad
        do i=1,nleaf
           erad_loc(i)=erad_loc(i)+uold(ind_leaf(i),8+j)
        enddo
     enddo
     ! Compute temperature 
!!!!! PH attention facteur (1-x)
     do i=1,nleaf
        
#if NSCHEM !=0                                                                                                                                                                                                                                                                                                               
        T2(i)=(gamma-1.0)*(T2(i)-ekin(i)-emag(i)-erad_loc(i))/(uold(ind_leaf(i),1) - uold(ind_leaf(i),neulS+1))
#else
        T2(i)=(gamma-1.0)*(T2(i)-ekin(i)-emag(i)-erad_loc(i))/(uold(ind_leaf(i),1))
#endif
        if(isnan(T2(i))) write(*,*) "WARNING: T2(i) is NaN: P,Ekin,Emag,Erad,dens,nH2", uold(ind_leaf(i),neulP), ekin(i), emag(i), erad_loc(i), uold(ind_leaf(i),1), uold(ind_leaf(i),neulS+1)
        if(T2(i) .LE. 0 ) write(*,*) "WARNING: T2(i) < 0: P,Ekin,Emag,Erad,dens,nH2", uold(ind_leaf(i),neulP), ekin(i), emag(i), erad_loc(i), uold(ind_leaf(i),1), uold(ind_leaf(i),neulS+1)
        if(isnan(T2(i)) .OR. T2(i) .LT. 0) stop
        
        T2(i) = T2(i)*scale_T2
        
     end do


!======================================================================

     !nH2 : just one value here
     do i=1,nleaf
        ! 1st: I read current value:
        ntot= uold(ind_leaf(i),1)             ! the ntot ???
        nH2 = uold(ind_leaf(i),neulS+1)       ! H2 density    

        if(myid .EQ. 1 .AND. nH2 .LT. 0) write(*,*) "WARNING nH2 LT 0:", ind, i, nH2, ntot

!        k2 = kph0 *uold(ind_leaf(i),neulS)*fshield    !kph0*<exp(-tau_(d,1000) )>*fshield
        if(isnan(uold(ind_leaf(i),neulS+2))) write(*,*) "WARNING (simple_chemostep1 BEFORE solve) uold( ,neulS+2) is nan, i, ind_leaf(i)",i, ind_leaf(i) 
        k2 = kph0 * uold(ind_leaf(i),neulS+2) ! kph_0*<fshield*exp(-tau_(d,1000) )> (eq 31, Glover&MacLow2006) 
        k1 = k1_0*sqrt(T2(i)/100.0_dp)/(1.0_dp + (T2(i)/464.0_dp)**1.5_dp)        ! k1 = k1_0*sqrt(T/100K)*S(T)

!        if(isnan(k2))  write(*,*) "WARNING (simple_chemostep1) k2 is NaN: uold=", uold(ind_leaf(i),neulS+2)
!        k2 = kph0
                                             
        ! 2nd: I solve the evolution of H2 for the next step
!        if(myid .EQ. 1) write(*,*) "CALLING SOLVE2_H2FORM"
        call solve2_H2form(nH2, ntot, k1, k2, dt_ilev)   !!! xH and xH2 must be in/out variables.
!        call solve_H2form(nH2, ntot, k1, k2, dt_ilev)   !!! xH and xH2 must be in/out variables.

        if(nH2/ntot .GT. 0.5 .OR. nH2/ntot .LT. 0.0) write(*,*) "WARNING: xH2,nH2,ntot", nH2/ntot, nH2,ntot
        if(isnan(nH2)) write(*,*)"WARNING(simple_chemostep1 AFTER): nH2 is nan, i, ind_leaf(i)",i, ind_leaf(i)

        ! 3rd: I actualize the value in the vector
        uold(ind_leaf(i),neulS+1) = nH2

     end do   ! leaf

  end do      ! cells

!  if(myid .EQ. 1) write(*,*) '***VAL: EXIT SIMPLE_CHEMOSTEP1, finished loop over cells'

end subroutine simple_chemostep1

!###########################################################                                          
!###########################################################                                          
!###########################################################                                          
!########################################################### 

subroutine solve_H2form(nH2, ntot, k1, k2, dt_ilev)
  implicit none

  real(kind=8),intent(INOUT)  :: nH2
  real(kind=8),intent(IN)     :: ntot, k1, k2, dt_ilev

  integer                     :: Nsteps, istep
  real(kind=8)                :: dtchem, dtstep, xH, xH2, xH2new, dt_tot, temps, dt   !nH2 = nH2/ntot
  real(kind=8)            :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  dt_tot = dt_ilev * scale_t ! * 3.08d18 / sqrt(kb_mm)
  dt = dt_tot/1000.0

!  xH  = nH  !/ntot
  xH2 = nH2/ntot
  xH = 1.0 - 2.0*xH2
  dtchem = 0.8/abs(k1 - k2*xH2)
  temps = 0.0

  ! EQUATION TO SOLVE:
  ! nH2(t+dt) = nH2(t) + (k1*nH**2 - k2*nH2)dt   

  
!!! to be written
! Idea???
  Nsteps = int(dt_ilev/dtchem) + 1   !DO IT BETTER. In the future use a do while, recalculate dt_chem at each timestep and compare as in cooling_fine. 

  do while ( temps < dt_tot)
!  do istep=1, Nsteps    

     xH2new = xH2 + (k1*xH*ntot - k2*xH2)*dt !*dtchem
     xH = 1.0 - 2.0*xH2new
     xH2 = xH2new
     if(xH2new .GT. 0.5) write(*,*) "WARNING: H2 fraction GT 0.5, xH, xH2", xH, xH2new
     if(xH2new .LT. 0.0) write(*,*) "WARNING: H2 fraction LT 0, xH, xH2", xH, xH2new
     if(xH .GT. 1.0) write(*,*) "WARNING: HI fraction GT 1, xH, xH2", xH, xH2new
     if(xH .LT. 0.0) write(*,*) "WARNING: HI fraction LT 0, xH, xH2", xH, xH2new

     temps = temps + dt

  end do

!  write(*,*) nH, nH2,dt_ilev, dt_tot, dt, dtchem, xH, xH2 
!  nH  = xH*ntot
  nH2 = xH2*ntot
end subroutine solve_H2form

!###########################################################                                          
!###########################################################                                          
!###########################################################                                          
!########################################################### 

subroutine solve2_H2form(nH2, ntot, k1, k2, dt_ilev)
  use amr_commons
  implicit none

  real(kind=8),intent(INOUT)  :: nH2
  real(kind=8),intent(IN)     :: ntot, k1, k2, dt_ilev

  integer                     :: Nsteps, istep
  real(kind=8)                :: dtchem, dtstep, xH, xH2, xH2new, dt_tot, temps, dt   !nH2 = nH2/ntot
  real(kind=8)                :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2,denom, k1ntot, dennew
  real(kind=8)                :: xH2min = 1.d-30, xH2max = 4.999d-1, xHmin = 1.d-30
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  if(isnan(nH2) .OR. isnan(ntot) .OR. isnan(k2) .OR. isnan(k1)) write(*,*) "WARNING: solve2_H2form invalid arguments nH2, ntot, k2, k1", nH2, ntot, k2, k1


  dennew = 1.0
  dt_tot = dt_ilev * scale_t ! * 3.08d18 / sqrt(kb_mm)
!  dt = dt_tot*1.d-5         ! Timestep seed
  dtchem = dt
!  denom = 1.0

!  xH  = nH  !/ntot
  if(myid .EQ. 1 .AND. nH2 .LT. 0.0) write(*,*) 'WARNING initial nH2 LT 0.0', nH2 
  xH2 = MAX(nH2/ntot,xH2min)
  if(myid .EQ. 1 .AND. xH2 .GT. 0.5) write(*,*) 'WARNING initial xH2 GT 0.5', xH2 
  xH2 = MIN(xH2, xH2max)
  xH = MAX(1.0 - 2.0*xH2, xHmin)
  
  k1ntot = k1*ntot
  if(isnan(k1ntot)) write(*,*) "WARNING: k1*ntot is nan, k1, ntot=", k1, ntot
  temps = 0.0

  !---------------------------------------------------------------
  ! EQUATION TO SOLVE:
  !( nH2(t+dt)- nH2(t))/dt = (k1*xH*n - k2*xH2)  
  !                        = k1*n*(1-2*xH2(t+dt)) - k2*xH2(t+dt)
  ! nH2(t+dt) =  (xH2(t) - k1*n*dt)/( 1 + (2*k1*n+k2)*dt)
  !---------------------------------------------------------------
!if(myid .EQ. 1) then
!  write(*,*) "=============================================="
!  write(*,*) "=====            TEST                    ====="
!  write(*,*) "       temps                   xH2            "
!endif

  !if(myid .EQ. 1) write(*,*) '***VAL: Entering SOLVE2_H2FORM'


  do while ( temps < dt_tot)



     denom= MAX((k1ntot - (2.0*k1ntot + k2)*xH2), 1./dt_tot)
     dtchem = MIN(ABS(xH2/denom), dt_tot)
     dtchem = 5d-2*MAX(1.0/(2*k1ntot + k2), dtchem)

!     dtchem = abs(xH2/(k1*ntot - (2.0*k1*ntot + k2)*xH2))      !PH estimation of time step
!     dtchem = 0.01*dtchem
!     if(myid .EQ. 1) write(*,*) '***VAL: SOLVE2_H2', denom, dtchem, dt_tot-temps

     dt = min(dtchem, dt_tot-temps)
     temps = temps + dt


!     if(myid.EQ.1) write(*,*)  temps, xH2
!     if(myid.EQ.1) write(*,*)  temps, dtchem, dt_tot-temps, dt

     dennew = (1.0 + (2.0*k1ntot + k2)*dt)
     if(isnan(dennew)) write(*,*) "WARNING: (IN) dennew is NaN", temps,k1ntot, k2, dt 
     xH2new = (xH2 + k1ntot*dt) / dennew              !(1.0 + (2.0*k1ntot + k2)*dt) !*dtchem
     !if(xH2new .GT. 0.5) xH2new = 0.5

     if(isnan(xH2new)) write(*,*) "WARNING: (IN) H2new fraction is NaN", temps, xH2, k1ntot, k2, dt  
     if(xH2new .GT. 0.5) write(*,*) "WARNING: (IN) xH2new GT 0.5,xH2", xH2new, MIN(MAX(xH2new,xH2min),0.5) 
     if(xH2new .LT. 0.0) write(*,*) "WARNING: (IN) xH2new LT 0,xH2", xH2new, MIN(MAX(xH2new,xH2min),0.5) 

     xH2new = MAX(xH2new,xH2min)
     xH2new = MIN(xH2new,xH2max)

     xH = 1.0 - 2.0*xH2new
     xH2 = xH2new




!     if(xH .GT. 1.0) write(*,*) "WARNING: (IN) HI fraction GT 1, xH, xH2", xH, xH2new
!     if(xH .LT. 0.0) write(*,*) "WARNING: (IN) HI fraction LT 0, xH, xH2", xH, xH2new


  end do

  if(isnan(xH2)) write(*,*) "WARNING: (OUT) H2 fraction is NaN"    
  if(xH2new .GT. 0.5 .OR. xH .LT. 0.0 ) write(*,*) "WARNING: (OUT) xH2 GT 0.5, xH, xH2", xH, xH2new
  if(xH2new .LT. 0.0 .OR. xH .GT. 1.0) write(*,*) "WARNING: (OUT) xH2 LT 0, xH, xH2", xH, xH2new
!  if(xH .GT. 1.0) write(*,*) "WARNING: HI (OUT) fraction GT 1, xH, xH2", xH, xH2new
!  if(xH .LT. 0.0) write(*,*) "WARNING: HI (OUT) fraction LT 0, xH, xH2", xH, xH2new



!  write(*,*) "=============================================="


!  write(*,*) nH, nH2,dt_ilev, dt_tot, dt, dtchem, xH, xH2 
!  nH  = xH*ntot

  if(isnan(xH2)) write(*,*) "WARNING: (END) xH2 is nan"
  xH2 = min(xH2, xH2max)
  xH2 = max(xH2, xH2min)
  nH2 = xH2*ntot

  if(nH2 .LT. 0.0) write(*,*) "WARNING END: ", ntot, xH2, nH2 

end subroutine solve2_H2form
