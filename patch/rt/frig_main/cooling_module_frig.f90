!! comes from frig version and then all the specific development done by Valeska
!! to take into account extinction have been moved there
!! PH 19/01/2017
!=======================================================================
subroutine solve_cooling_frig(nH,T2,zsolar,boost,dt,deltaT2,ncell)
!=======================================================================
  implicit none
  ! BRIDGE FUNCTION WITH SAME INTERFACE AS 
  ! Input/output variables to this function
  ! nH - hydrogen number density in PHYSICAL units
  ! T2 - temperature / mu in PHYSICAL units
  ! zsolar - Metallicity in solar units (Zphys / 0.02)
  ! boost - raditation boost - exp(-100*nH) if self_shielding=True
  ! dt - cooling timestep in seconds
  ! deltaT2 - temperature change in K/mu (??)
  ! ncell - number of elements in the vector
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar,boost
  ! Input/output variables to analytic function calc_temp 
  real(kind=8)::NN,TT, dt_tot_unicode
  ! Temporary variables
  integer::i
  real(kind=8)::TT_ini, mu
  ! Units
  real(kind=8) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  ! HARD-CODED mu TO MAKE TEMPERATURE AGREE WITH HENNEBELLE CODE
  mu = 1.4
  scale_T2 = scale_T2 * mu
  do i=1,ncell
     NN = nH(i) ! NOTE!! THE CODE BELOW ASSUMES scale_nH=1 !!
                ! SO WE LEAVE THIS AS IT IS TO KEEP UNITS CONSISTENCY
     TT = T2(i) / scale_T2
     TT_ini = TT
     dt_tot_unicode = dt / scale_t
     call calc_temp(NN,TT,dt_tot_unicode)
     deltaT2(i) = (TT - TT_ini) * scale_T2
  end do
end subroutine solve_cooling_frig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine  calc_temp(NN,TT,dt_tot_unicode)
    use amr_parameters
    use hydro_commons

    implicit none

    integer :: n,i,j,k,idim, iter, itermax,ii

    real(dp) :: dt, dt_tot, temps, dt_max, itermoy
    real(dp) :: rho,temp,dt_tot_unicode

    !alpha replaced by alpha_ct because of conflict with another alpha by PH 19/01/2017
    real(dp) :: mm,uma, kb, alpha_ct,mu,kb_mm
    real(dp) :: NN,TT, TTold, ref,ref2,dRefdT, eps, vardt,varrel, dTemp
    real(dp) :: rhoutot2
    real(dp) :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
    ! HARD-CODED mu TO MAKE TEMPERATURE AGREE WITH HENNEBELLE CODE
    mu = 1.4
    !
    ! Cette routine fonctionne en cgs
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    kb  =  1.38062d-16   ! erg/degre
    !  uma =  1.660531e-24  ! gramme
    !  mu  =  1.4
    !  mm = mu*uma
    !  kb_mm = kb / mm
    !  TT = TT  / kb  !/ kb_mm

    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_T2 = scale_T2 * mu

    if( TT .le. 0.) then
        TT = 50. / scale_T2
        return
    endif

    !if( TT*scale_T2 .gt. 50.) then
    !TT = 50. / scale_T2
    !return
    !endif

    vardt = 10.**(1./10.); varrel = 0.2

    dt_tot = dt_tot_unicode * scale_t ! * 3.08d18 / sqrt(kb_mm)
    TT     = TT * scale_T2

    !  nn = (rho/(gramme/cm3)) /mm

    itermax = 0 ; itermoy = 0.



    if (NN .le. smallr) then
        if( NN .le. 0)  write(*,*) 'prob dens',NN
        NN = smallr  !max(NN,smallr)
    endif


    ! alpha_ct = NN*kb_mm/(gamma-1.)
    alpha_ct = NN*kb/(gamma-1.)

    ! eps - a small offset of T to find gradient in T
    eps = 1d-5

    iter  = 0 ; temps = 0.
    do while ( temps < dt_tot)
        if (TT .lt.0) then
            write(*,*) 'prob Temp',TT, NN
            !         write(*,*) 'repair assuming isobariticity'
            NN = max(NN,smallr)
            TT = min(4000./NN,8000.)  !2.*4000. / NN
        endif


        TTold = TT

        ! Calculate cooling rate
        !NN is assumed to be in cc and TT in Kelvin
        if (TT < 10035.d0) then
            call cooling_low(TT,NN,ref)
            call cooling_low(TT*(1d0+eps),NN,ref2)
        else
            call cooling_high(TT,NN,ref)
            call cooling_high(TT*(1d0+eps),NN,ref2)
        end if
        
        ! dT = T*(1+eps)-T = eps*T
        dRefdT = (ref2-ref)/(TT*eps)

        ! TODO STG - COPY THIS FUNCTION UP TO HERE, USE ref, drefdT TO 
        !            REPLACE rt_cmp_metals SOMEHOW


        !       write(*,*) 'check',TTold, TT,NN,ref,dRefdT,iter


        if (iter == 0) then
            if (dRefDT .ne. 0.) then
                dt = abs(1.0E-1 * alpha_ct/dRefDT)
            else
                dt = 1.0E-1 * dt_tot
            endif
            dt_max = dt_tot - temps
            if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
        endif

        dTemp = ref/(alpha_ct/dt - dRefdT)

        eps = abs(dTemp/TT)
        if (eps > 0.2) dTemp = 0.2*TTold*dTemp/abs(dTemp)

        TT = TTold + dTemp
        if (TT < 0.) then
            write(*,*) 'Temperature negative !!!'
            write(*,*) 'TTold,TT   = ',TTold,TT
            write(*,*) 'rho   = ',rho
            TT = 100.  !*kelvin
        endif


        iter = iter + 1

        temps = temps + dt

        dt = vardt*varrel*dt/Max(vardt*eps, varrel)

        dt_max = dt_tot - temps
        if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
        !        write(*,987) temps, TT
        !987     format(E10.3,2x,E10.3)
        !        read (*,*)
    enddo


    !  if (TT .ge. 50.)  TT=50.

    !!now convert temperature in code units
    TT = TT / scale_T2

    return
end subroutine calc_temp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cooling_high(T,n,ref)
    use amr_parameters
    implicit none

    real(dp) :: T,P,N,x,ne                    !x est le taux d'ionisation
    real(dp) :: T2, ref2
    real(dp) :: froid,chaud,ref,nc, froidc
    real(dp) :: froidcII, froido, froidh
    real(dp) :: froidc_m,froidcII_m,froido_m
    real(dp) :: param, G0, epsilon,k1,k2,bet,froidrec
    real(dp) :: eps

    real(dp) :: logT, intcst, logT2
    real(dp) :: ion, neut
    real(dp) :: refion

    !taux de refroidissement base sur Dopita et Sutherland

    logT=log10(T)

    if (logT .LT. 4.0) then
        froid=0.1343*logT**3-1.3906*logT**2+5.1554*logT-31.967
    else if (logT .LT. 4.25) then
        froid=12.64*logT-75.56
    else if (logT .LT. 4.35) then
        froid=-0.3*logT-20.565
    else if (logT .LT. 4.9) then
        froid=1.745*logT-29.463
    else if (logT .LT. 5.4) then
        froid=-20.9125
    else if (logT .LT. 5.9) then
        froid=-1.795*logT-11.219
    else if (logT .LT. 6.2) then
        froid=-21.8095
    else if (logT .LT. 6.7) then
        froid=-1.261*logT-13.991
    else
        froid=-22.44
    endif

    froid=-1.0*10.0**(froid)

!    chaud = 1.E-25
    chaud = 0.
    ref= chaud*n + (n**2)*(froid)

end subroutine cooling_high

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine cooling_low(T,n,ref)

  use amr_parameters

  implicit none

  real(dp) :: T,P,N,x,ne,x_ana                 !x est le taux d'ionisation
  real(dp) :: froid,chaud,ref,nc, froidc
  real(dp) :: froidcII, froido, froidh
  real(dp) :: froidc_m,froidcII_m,froido_m
  real(dp) :: param, G0, epsilon,k1,k2,bet,froidrec


  !fonction de chauffage et de refroidissement calculee a partir des
  !refroidissements des differents elements



  !abondance de carbone 3.5 10-4, depletion 0.4

!!! on calcule l'ionisation
!!! on suppose que si x est superieure a 1.d-4 elle est domine par
!!! l'hydrogene et que la valeur inf est donne par le carbone
!!! et vaut 3.5 1.d-4 * depletion * densite

!!! Pour les electrons dus a l'hydrogene on prend la
!!! formule donnee par Wolfire et al. 2003 appendice C2.
!!! on prend un taux d'ionisation de 1.d-16 G0'=GO/1.7
!!! Z'd = 1 et phipah=0.5



  ne = 2.4d-3*((T/100d0)**0.25d0)/0.5d0 !formule C15 Wolfire et al. 2003

  ! Analytic ionisation in absence of photoionisation
  x_ana = ne / N   ! ionisation
  x_ana = min(x_ana,0.1d0)
  x_ana = max(x_ana,3.5d-4*0.4d0)
  x = x_ana ! (Different variables in case we need photoionised values later)

  !transition hyperfine a basse temperature: carbone et oxygene
  !chiffre pris dans la these de Karl Joulain

  !refroidissement par le carbone




  !      froidcII = ( 2.2d-23                     !excitation par H
  !     c          + 5.5d-20 * 2 / sqrt(T) * x ) !excitation par e
  !     c              * 3.5d-4 * 0.4d0 * exp(-92.d0 / T)


 ! NOTE - HERE WE USE THE NON-PHOTOIONISED RATES AS THIS MIGHT 
 !        BE TOO HIGH AT x=1
  froidcII =  92. * 1.38E-16 * 2. * (2.8E-7* ((T/100.)**(-0.5))*x_ana + 8.E-10*((T/100.)**(0.07))) &
       * 3.5E-4 * 0.4 * exp(-92./ T)
  !     c               3.d-4  * exp(-92. / T)


  !refroidissement par l'oxygene
  !abondance 8.6 10-4 depletion 0.8

  froido = 1.E-26 * sqrt(T) * (24. * exp(-228./ T) + 7. * exp(-326./ T) )


  !      froido =  230.d0*1.38d-16 * (
  !     c            1.4d-8*x + 9.2d-11 *(T /100.d0)**(0.67) )
  !     c     * exp(-230.d0 / T)

  !      froido = froido +  330.d0*1.38d-16 *(
  !     c            1.4d-8*x + 4.3d-11 *(T /100.d0)**(0.8) )
  !     c     * exp(-330.d0 / T)

  !      froido = froido +  98.d0*1.38d-16 * (
  !     c            5.d-9 *x + 1.1d-10* (T /100.d0)**(0.44) )
  !    c      * exp(-98.d0 / T)


  !       froido = 2.5d-27 * (T/100)**0.4 * exp(-228.d0 / T)


  !on tient compte de l'abondance du
  froido = froido * 4.5E-4


  !refroidissement par l'hydrogene
  ! formule de Spitzer 1978
  ! NOTE - RT function handles hydrogen cooling out of equilibrium
  froidh = 0d0
  if (.not. rt) then
     froidh = 7.3E-19 * x * exp(-118400./ T )
  endif

  !refroidissement par les raies metastables des metaux
  !chiffre pris dans Hollenbach and McKee 1989 (ApJ 342, 306)





  !carbone une fois ionise ,1 transitions 2P 4P
  ! le poids est 1
  ! 2P->4P :
  ! les expressions des coefficients d'excitation ont une dependance
  !en la temperature differente au dela de 10000K
  !abondance 3.5 d-4 depletion 0.4

         froidcII_m = 6.2d4 * 1.38d-16 * 1.d0 * &    !transition 2P->4P
        ( 2.3d-8* (T/10000.)**(-0.5) * x + 1.d-12 ) *exp(-6.2d4 / T) &
           * 3.5d-4 * 0.4




         if ( T .le. 1.d4 ) then
         froido_m = 2.3d4 * 1.38d-16 / 3.d0 * &
        ( 5.1d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.3d4/T)
  
         froido_m = froido_m + &
              4.9d4 * 1.38d-16 / 3.d0  * &
        ( 2.5d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-4.9d4/T)
  

         froido_m = froido_m + &
              2.6d4 * 1.38d-16 * 1.d0  * &
        ( 5.2d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.6d4/T)

         else

         froido_m = 2.3d4 * 1.38d-16 / 3.d0 * &
        ( 5.1d-9 * (T/10000.)**(0.17) * x + 1.d-12) *exp(-2.3d4/T)
  
         froido_m = froido_m + &
              4.9d4 * 1.38d-16 / 3.d0  * &
        ( 2.5d-9 * (T/10000.)**(0.13) * x + 1.d-12) *exp(-4.9d4/T)


         froido_m = froido_m + &
              2.6d4 * 1.38d-16 * 1.d0  * &
        ( 5.2d-9 * (T/10000.)**(0.15) * x + 1.d-12) *exp(-2.6d4/T)


         endif

  !! abondance de l'oxygene
         froido_m = froido_m *   4.5d-4



!!! on somme les refroidissements
  froid = froidcII  + froidh  + froido  + froido_m +  froidcII_m


  !      froid=froid*1.d-13    !conversion en MKS


  !refroidissement par le carbone neutre. On suppose l'equilibre
  ! de la reaction C + hv <=> C+ + e-
  ! les taux de reactions et de refroidissement sont pris dans
  !la these de Karl Joulain.

  ! abondance du carbone relative a n (MKS)


  !    C+ + e- => C
  !       k1 = 4.4d-12 * (T/300.)**(-0.61) !s^-1 cm^-3

  !       k1 = k1


  !    C => C+ + e-
  !       k2 = 2.2d-10


  ! on a : [C] = k1/k2 [C+] * [e-]
  ! on suppose que tout le carbone est sous forme C+
  ! et que [e-] = [C+]

  ! l'abondance relative de carbone
  !      nc = k1/k2 * (3.5d-4*0.4)**2 * n


  !      froidc =  1.0d-24 * ( 1.4d0 * exp( -23.d0 / T ) &
  !                     + 3.8d0 * exp( -62.d0 / T )   )

  !      froidc = froidc * nc !(nc est l'abondance relative du carbone)


  !       n=exp(log(10.d0)*logn) !ici pas besoin de log

  !       valeur utilisees au 23/08/98
  !       chaud=4.d0*exp(-24.5d0*log(10.d0))*1.d-7  !un peu empirique ....



!!!! on calcule le chauffage
!!! on prend en compte le chauffage sur les grains
!!! formules 1 et 2  de Wolfire et al. 1995

!!!! G0 est le flux UV par rapport au flux defini par Habing et
!!!! Draine

  G0 = 1./1.7

  param = G0 * sqrt(T)/(n*x)
  epsilon = 4.9E-2 / (1. + (param/1925.)**0.73)
  epsilon  = epsilon + 3.7E-2 * (T/1.E4)**0.7 / (1. + (param/5.E3) )


  chaud = 1.E-24 * epsilon


  ! pour un flux de rayonnement G0/1.7
  chaud = chaud * G0

  !refroidissement recombinaison sur les grains charges positivement
  bet = 0.74/(T**0.068)
  froidrec = 4.65E-30*(T**0.94)*(param**bet)*x



  !! chaud=1.d-32 !un peu empirique ....

  !      froidc=0.d0


  ref= chaud*n - (n**2)*(froid + froidrec) !!!+ froidc)

  return

end subroutine cooling_low

!###########################################################
!###########################################################
!###########################################################
!###########################################################
!! ph 11/2010 
!! estimate the column density through the box in 6 directions
!! from the edges of a grid using the tree structure of RAMSES. 
!! To get the full column density one must add the contribution 
!! of the cells within the oct. All information is local 
!! and do not require communication between cpus
!! VV 2012 
!! Modified version in order to take into account the contribution
!! of all the cells contained in the cubic shells. It considers 
!! directions based on spherical projections. 



subroutine column_density(ind_grid,ngrid,ilevel,column_dens, H2column_dens)
  use amr_commons  
  use hydro_commons
  use cooling_module
  implicit none

  integer,dimension(1:nvector)                          :: ind_grid
  integer                                               :: ilevel,ngrid
  real(dp),dimension(1:nvector,1:NdirExt_m,1:NdirExt_n) :: column_dens, H2column_dens

  integer,parameter                                     :: ndir=6       
  integer                                               :: i,ind,idim 
  integer                                               :: ix,iy,iz,idir,il,twice

  real(dp),dimension(1:nvector,1:3, 1:ndir)             :: xpart
  real(dp),dimension(1:nvector,1:3)                     :: xpart_int
  real(dp)                                              :: dx,dx_loc, xgaux
  integer                                               :: nx_loc,pos_son,ind_father

  real(dp),dimension(1:twotondim,1:3)                   :: xc
  real(dp),dimension(1:nvector,1:ndim)                  :: xx,xx_check
  real(dp),dimension(1:nvector)                         :: dist

  real(dp),dimension(1:nvector,1:ndir)                  :: col_dens,col_check 
  integer,dimension(1:nvector)                          :: cell_index,cell_levl
  real(dp),dimension(3,6)                               :: shift
  integer, dimension(1:ngrid,1:ndir)                    :: ind_lim1, ind_lim2

  
  if(numbtot(1,ilevel)==0)return
  
  !define various things needed to get the distance
  
  ! Mesh size at level ilevel in coarse cell units
  dx=0.5_dp**ilevel
  
  ! Rescaling factors
  !  nx_loc=(icoarse_max-icoarse_min+1)
  !  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  !  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  !  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  !  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  !  scale=dble(nx_loc)!/boxlen
  
  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5_dp)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5_dp)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5_dp)*dx
  end do

    
  !  ! this part checks that grids can find themselves and verify that coordinates are well recovered
  !  ! you must uncomment it if you need it
  !         do idim=1,ndim
  !           do i=1,ngrid
  !             xpart(i,idim)=xg(ind_grid(i),idim)
  !           end do
  !         end do
  !  
  !           call get_cell_index(cell_index,cell_levl,xpart,ilevel-1,ngrid)
  !  
  !            do i=1,ngrid
  !             if( son(cell_index(i)) .ne. ind_grid(i) ) then
  !              write(*,*) 'prob cell',cell_index(i),son(cell_index(i)),ind_grid(i),ilevel-1,cell_levl(i),i
  !             endif
  !  
  !                !get the father of the cell
  !                ind_father = mod(cell_index(i) -ncoarse,ngridmax)  !father(cell_index(i)) 
  !                !get the cell position in its oct
  !                pos_son = (cell_index(i)-ncoarse-ind_father)/ngridmax + 1
  !  
  !                !calculate the cell position
  !                write(*,*) 'dx',dx*2.
  !                do idim=1,ndim           
  !                  xx_check(i,idim)=xg(ind_father,idim)+xc(pos_son,idim)*2.
  !                enddo
  !                if(myid .eq. 1 .and. (i .eq. 1 .or. i .eq. 2) ) then 
  !                  write(*,*) 'xx_check',xx_check(i,1)-xg(ind_grid(i),1),xx_check(i,2)-xg(ind_grid(i),2),xx_check(i,3)-xg(ind_grid(i),3)
  !                endif
  !            end do 
  
  
  ! define the path direction 
  ! first index is for x,y,z while second is for direction
  ! x,-x,y,-y,z,-z in this order 
  shift(1,1) = 1. ; shift(2,1) = 0. ; shift(3,1) = 0. 
  shift(1,2) =-1. ; shift(2,2) = 0. ; shift(3,2) = 0. 
  shift(1,3) = 0. ; shift(2,3) = 1. ; shift(3,3) = 0. 
  shift(1,4) = 0. ; shift(2,4) =-1. ; shift(3,4) = 0. 
  shift(1,5) = 0. ; shift(2,5) = 0. ; shift(3,5) = 1. 
  shift(1,6) = 0. ; shift(2,6) = 0. ; shift(3,6) =-1.
  
  
  column_dens(:,:,:) = 0.
  H2column_dens(:,:,:) = 0.
  ind_lim1(:,:)= 0
  ind_lim2(:,:)= 0
  

  !-----    DETERMINATION OF BOUNDS  -----    
  !  Here we calculate the limits for the cubic shells considered at each 
  !  resolution level. We obtain external and internal limits.

  do idir=1,ndir        !  loop over the 6 directions
     do idim=1,ndim
        do i=1,ngrid
           xpart(i,idim,idir)=xg(ind_grid(i),idim)+shift(idim,idir)*dx
           col_dens(i,idir)=0.
           col_check(i,idir)=0.
         end do   !i
      end do      !idim
   end do         !idir

   

   do il = ilevel-1,1,-1            !  loop recursively over level
      dx_loc = 0.5_dp**(il)
      do idir=1,ndir                ! loop over the 6 directions
         do i=1,ngrid
            do idim=1,ndim
               if(shift(idim,idir).NE.0) then
                  xgaux = dx_loc*(INT(xg(ind_grid(i),idim)/dx_loc) + 0.5_dp)
                  ind_lim1(i,idir)= shift(idim,idir)*INT(ABS(xpart(i,idim,idir)-xgaux)/dx_loc )                  
               end if
            end do                  !idim
         end do                     !i
         
         !--- we do it twice because of oct structure           
         do twice=1,2               !do it twice
            do idim=1,ndim
               do i=1,ngrid
                  xpart(i,idim,idir)= xpart(i,idim,idir)+shift(idim,idir)*(dx_loc/2.)*twice
                  xpart_int(i,idim)= xpart(i,idim,idir)
               end do   !i
            end do      !idim
            
            call get_cell_index3(cell_index,cell_levl,xpart_int,il,ngrid)
            
            !--- NOT NECESSARY ---
            do i=1,ngrid ! loop over grid to calculate the column density
               if (cell_index(i) .ne. -1) then
                  if( cell_levl(i) .ne. il) then
                     write(*,*) 'problem in the calculation of column density'
                     write(*,*)  'cell_levl(i),il,ilevel',cell_levl(i),il,ilevel
                     stop
                  endif

                  col_dens(i,idir) = col_dens(i,idir) + dx_loc*uold(cell_index(i),1)
                  col_check(i,idir) = col_check(i,idir) + dx_loc
               end if
            end do !end loop over grid
            !----------- not necessary ---              
            
         end do        !end of do it twice
         
         
         !move the particle at the edge of the cell at level il
         do idim=1,ndim
            do i=1,ngrid

               if(shift(idim,idir) .NE. 0) then
                  xgaux = dx_loc*(INT(xg(ind_grid(i),idim)/dx_loc) + 0.5_dp)
                  ind_lim2(i,idir)= shift(idim,idir)*INT(ABS(xpart(i,idim,idir)-xgaux)/dx_loc )
               end if

               xpart(i,idim,idir)= xpart(i,idim,idir)+shift(idim,idir)*dx_loc/2.
               xpart_int(i,idim)= xpart(i,idim,idir)
                              
            end do
         end do
         
         
         !now before starting the next iteration one must check whether the particle is 
         !at the interface of a cell at level il-1 or whether one is in the center of such a cell
         !in the first case it is fine in the other case, we must jump from another
         !dx_loc to reach the next boundary of the next cell at level il-1
         
         !get cells at level il-1
         if (il .eq. 1) exit
         call get_cell_index3(cell_index,cell_levl,xpart_int,il-1,ngrid)
         
         do i=1,ngrid
            if (cell_index(i) .ne. -1)  then
               !get the father of the cell
               ind_father = mod(cell_index(i) -ncoarse,ngridmax)  !father(cell_index(i)) 
               !get the cell position in its oct
               pos_son = (cell_index(i)-ncoarse-ind_father)/ngridmax + 1
               
               !calculate the cell position
               !note that in principle pos_son is enough to get the requested information
               !here we choose to calculate the coordinate instead.            
               do idim=1,ndim           
                  xx_check(i,idim)=xg(ind_father,idim)+xc(pos_son,idim)*(0.5**(il-1-ilevel))
               enddo
            end if
         end do                             !ok
         
         
         !now calculate the distance between these cells and the particle
         !along the direction of interest
         dist=0
         do idim=1,ndim           
            do i=1,ngrid
               if (cell_index(i) .ne. -1)  dist(i) = dist(i) + (( xx_check(i,idim) - xpart(i,idim,idir) )*shift(idim,idir))**2
            end do
         end do
         do i=1,ngrid
            dist(i)=sqrt(dist(i))
         end do
         
         
         !finally if we are not at an interface, we move the particle to the next interface
         do i=1,ngrid
            if( dist(i) .lt. dx_loc/4.) then
               do idim=1,ndim
                  if (cell_index(i) .ne. -1) then
                     xpart(i,idim,idir)=xpart(i,idim,idir)+shift(idim,idir)*dx_loc
                     
                     !----------------------------------------------------------
                     ! we set the outer limit for the cubic shell
                     if(shift(idim,idir) .NE. 0) then
                        xgaux = dx_loc*(INT(xg(ind_grid(i),idim)/dx_loc) + 0.5_dp)
                        ind_lim2(i,idir)= shift(idim,idir)*INT(ABS(xpart(i,idim,idir)-xgaux)/dx_loc )
                     end if
                     !----------------------------------------------------------
                  endif
               end do
               
               !--- Check ---
               if (cell_index(i) .ne. -1) then
                  col_dens (i,idir) = col_dens (i,idir) + dx_loc*uold(cell_index(i),1)
                  col_check(i,idir) = col_check(i,idir) + dx_loc
               endif
               !--- check ---
               
            end if                          ! dist
         end do                             ! i
      end do                                ! end of loop over direction idir
      
      !-----  ADD CONTRIBUTIONS TO COLUMN DENSITY  -----------------
      ! once we have calculated the limits we call the contribution routine,
      ! that sums up to the column density

      call  contribution(ind_grid, ngrid, ilevel, il, ind_lim1, ind_lim2, column_dens,H2column_dens)
      
      
   end do                                   !end of loop recursively over level
   
   
   !now check that the whole grid has been properly covered
   do i=1,ngrid !end loop over grid
      !       if(myid .eq. 1 .and. i .eq. 1) write(*,*) 'col tot x',col_check(i,1)+col_check(i,2)+dx*2. 
      if(col_check(i,1)+col_check(i,2)+dx*2. .ne. 1.) then 
         write(*,*) 'col tot x',col_check(i,1)+col_check(i,2)+dx*2. 
         write(*,*) 'col tot x',col_check(i,1),col_check(i,2),dx*2. 
         stop
      endif
      if(col_check(i,3)+col_check(i,4)+dx*2. .ne. 1.) then 
         write(*,*) 'col tot y',col_check(i,3)+col_check(i,4)+dx*2. 
         write(*,*) 'col tot y',col_check(i,3),col_check(i,4),dx*2. 
         stop
      endif
      if(col_check(i,5)+col_check(i,6)+dx*2. .ne. 1.) then 
         write(*,*) 'col tot z',col_check(i,5)+col_check(i,6)+dx*2. 
         write(*,*) 'col tot z',col_check(i,5),col_check(i,6),dx*2. 
         stop
      endif
   enddo !end loop over grid
   
   
 end subroutine column_density
!#########################################################
!#########################################################
!#########################################################
!#########################################################
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! routine provided by Romain on July 2010 and included by PH on November 2010

subroutine get_cell_index3(cell_index,cell_levl,xpart,ilevel,np)
  use amr_commons
  implicit none
  integer                                :: np,ilevel
  integer,dimension(1:nvector)           :: cell_index,cell_levl
  real(dp),dimension(1:nvector,1:3)      :: xpart
  ! This function returns the index of the cell, at maximum level
  ! ilevel, in which the input particle sits
  real(dp)                               :: xx,yy,zz
  integer                                :: i,j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0
  
  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     call clean_stop
  end if
  
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  do i=1,np
     xx = xpart(i,1) + (nx-1)/2.0
     yy = xpart(i,2) + (ny-1)/2.0
     zz = xpart(i,3) + (nz-1)/2.0
     
     if( ((xx .le. 0) .or. (xx .ge. 1.)) .or. ((yy .le. 0) .or. (yy .ge. 1.)) .or. ((zz .le. 0) .or. (zz .ge. 1.)) ) then 
        cell_index(i)=-1.
     else 
        igrid=igrid0
        do j=1,ilevel
           ii=1; jj=1; kk=1
           if(xx<xg(igrid,1))ii=0
           if(yy<xg(igrid,2))jj=0
           if(zz<xg(igrid,3))kk=0
           ind=1+ii+2*jj+4*kk
           iskip=ncoarse+(ind-1)*ngridmax
           ind_cell=iskip+igrid
           igrid=son(ind_cell)
           if(igrid==0.or.j==ilevel)  exit
        end do
        cell_index(i)=ind_cell
        cell_levl(i)=j
     endif
  end do
  
  
end subroutine get_cell_index3


!#########################################################
!#########################################################
!#########################################################
!#########################################################
!! 2012 VV 
!! Modified version of get_cell_index:
!! This subroutine gives the cell index of just one cell
!! This function returns the index of one cell, at maximum level
!! ilevel, in which the input particle sits

subroutine get_cell_index2(cell_index2,cell_levl2,xpart2,ilevel)
  use amr_commons
  implicit none

  integer                  ::cell_index2,cell_levl2,ilevel
  real(dp),dimension(1:3)  ::xpart2
  real(dp)                 ::xx,yy,zz
  integer                  ::j,ii,jj,kk,ind,iskip,igrid,ind_cell,igrid0
  
  if ((nx.eq.1).and.(ny.eq.1).and.(nz.eq.1)) then
  else if ((nx.eq.3).and.(ny.eq.3).and.(nz.eq.3)) then
  else
     write(*,*)"nx=ny=nz != 1,3 is not supported."
     call clean_stop
  end if
  
  igrid0=son(1+icoarse_min+jcoarse_min*nx+kcoarse_min*nx*ny)
  
  xx = xpart2(1) + (nx-1)/2.0
  yy = xpart2(2) + (ny-1)/2.0
  zz = xpart2(3) + (nz-1)/2.0
  
  if( ((xx .le. 0) .or. (xx .ge. 1.)) .or. ((yy .le. 0) .or. (yy .ge. 1.)) .or. ((zz .le. 0) .or. (zz .ge. 1.)) ) then
     cell_index2=-1.
  else
     igrid=igrid0
     do j=1,ilevel
        ii=1; jj=1; kk=1
        if(xx<xg(igrid,1))ii=0
        if(yy<xg(igrid,2))jj=0
        if(zz<xg(igrid,3))kk=0
        ind=1+ii+2*jj+4*kk
        iskip=ncoarse+(ind-1)*ngridmax
        ind_cell=iskip+igrid
        igrid=son(ind_cell)
        if(igrid==0.or.j==ilevel)  exit
     end do
     cell_index2=ind_cell
     cell_levl2=j
  endif
  
end subroutine get_cell_index2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 02.04.2012 VV
!! Contribution to directions. Calculation of "column_dens" used in
!! subroutine column_density
!! It is here where we actually add the contributions of the cells in the shell
!! to the column densities

subroutine contribution(ind_grid, ngrid, ilevel, il, ind_lim1, ind_lim2, column_dens, H2column_dens)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none

  integer,parameter                                                   :: ndir=6  
  integer, intent (in)                                                :: ilevel, il,ngrid
  integer,dimension(1:nvector)                                        :: ind_grid
  integer, dimension(1:ngrid,1:ndir), intent (in)                     :: ind_lim1, ind_lim2
  real(dp),dimension(1:nvector,1:NdirExt_m,1:NdirExt_n),intent(inout) :: column_dens, H2column_dens

  !!take care "nent" should be added here. It is fine if it is 0 of course 
  integer                                               :: neulS=8+nrad+nextinct

  integer                                               :: i, inx,iny,inz, deltam
  integer                                               :: m, n, mloop, nloop, nl, mn
  real(dp)                                              :: dx_loc, dx_cross_ext, halfdx, quartdx
  integer, dimension(1:NdirExt_n)                       :: deltan1, deltan2
  real(dp),dimension(1:3)                               :: xgcoarse, xoct
  integer, dimension(1:3)                               :: indaux, indaux2
  integer                                               :: ind_oct, ind_oct1, ind_oct2, reg
  real(dp),dimension(1:3)                               :: xpart2, xzero
  real(dp)                                              :: xg_il, yg_il, zg_il     !xg at level il
  integer                                               :: cell_ind2, cell_levl2
  integer,dimension(1:6,1:3)                            :: l_inf, l_sup

  integer                                               :: ix, iy,iz

  !---  m, n loop limits  ------------------------------------------
  ! here we define the limits for the loops around the direction to the cell center
  ! For the vertical directions (m =0, NdirExt_m), the azimuthal angle covers 2*pi,
  ! for other directions we use +/- 1/8 of the total directions
  !----------------------------------------------------------------
  do m = 1, NdirExt_m
     ! for the vertical directions 
     if(m .EQ. 1 .OR. m .EQ. NdirExt_m) then
        deltan1(m) = INT(NdirExt_n/2.0_dp)
        deltan2(m) = deltan1(m) -1
     else
        deltan1(m) = INT(NdirExt_n/8)
        deltan2(m) = deltan1(m)
     end if
  end do
  deltam = INT((NdirExt_m-1)/4.)

  dx_loc = 0.5_dp**(il)
  halfdx  = dx_loc/2_dp
  quartdx = dx_loc/4.0_dp


  !-------------------------------------------------!
  !       CONTRIBUTION TO DIRECTIONS                !
  !-------------------------------------------------!

  !THETA AND PHI   
  
  !LOOP over cells
  do i=1,ngrid
     
     xzero(1) = xg(ind_grid(i),1)       ! xzero is the grid center 
     xzero(2) = xg(ind_grid(i),2)  
     xzero(3) = xg(ind_grid(i),3)  
     
     !! we have to determine in which configuration the target cell and the treated cell are
     !! in order to use the precalculated values for the geometrical corrections.
     !! the configuration depends on the difference of levels and their relative positions in the octs

     !-----  OBTAINING ind within the oct   -----         
     if(il .EQ. ilevel-1) then
        ind_oct = 73
        
     else if(il .EQ. ilevel -2) then
        do inx =1, 3
           xoct(inx)     = (INT(xzero(inx)/dx_loc) +0.5_dp)*dx_loc
           indaux(inx)   = INT((xzero(inx)-xoct(inx))/halfdx +0.5_dp)
        end do
        ind_oct = indaux(1) + 1 + 2*indaux(2) + 4*indaux(3)
        ind_oct = ind_oct + 64
        
     else         !if(il .GE. ilevel -3) then
        
        do inx =1, 3
           xgcoarse(inx) = (INT(xzero(inx)/halfdx) +0.5_dp)*halfdx
           xoct(inx)     = (INT(xzero(inx)/dx_loc) +0.5_dp)*dx_loc
           indaux2(inx)  = INT((xzero(inx)-xgcoarse(inx))/quartdx +0.5_dp)
           indaux(inx)   = INT((xgcoarse(inx)-xoct(inx))/halfdx +0.5_dp)
        end do
        ind_oct1 = indaux(1) + 1 + 2*indaux(2) + 4*indaux(3)    !coarse level
        ind_oct2 = indaux2(1) + 1 + 2*indaux2(2) + 4*indaux2(3)    !fine level
        ind_oct = (ind_oct1-1)*8 + ind_oct2
        
     end if


     !-----   SPLITTED LOOPS to avoid the "if" for skipping treated cells   -----
     ! we cut the cubic shell in 6 regions
     ! here we define the limits for the splitted regions and we integrate the column density.
     ! 6 REGIONS :

     l_inf(1,3) = ind_lim2(i,6)  
     l_sup(1,3) = ind_lim1(i,6)-1
     l_inf(1,2) = ind_lim2(i,4)
     l_sup(1,2) = ind_lim2(i,3)
     l_inf(1,1) = ind_lim2(i,2)
     l_sup(1,1) = ind_lim2(i,1)

     l_inf(2,3) = ind_lim1(i,5)+1
     l_sup(2,3) = ind_lim2(i,5)
     l_inf(2,2) = ind_lim2(i,4)
     l_sup(2,2) = ind_lim2(i,3)
     l_inf(2,1) = ind_lim2(i,2)
     l_sup(2,1) = ind_lim2(i,1)

     l_inf(3,3) = ind_lim1(i,6)
     l_sup(3,3) = ind_lim1(i,5)
     l_inf(3,2) = ind_lim2(i,4)
     l_sup(3,2) = ind_lim1(i,4)-1
     l_inf(3,1) = ind_lim2(i,2)
     l_sup(3,1) = ind_lim2(i,1)

     l_inf(4,3) = ind_lim1(i,6)
     l_sup(4,3) = ind_lim1(i,5)
     l_inf(4,2) = ind_lim1(i,3)+1
     l_sup(4,2) = ind_lim2(i,3)
     l_inf(4,1) = ind_lim2(i,2)
     l_sup(4,1) = ind_lim2(i,1)

     l_inf(5,3) = ind_lim1(i,6)
     l_sup(5,3) = ind_lim1(i,5)
     l_inf(5,2) = ind_lim1(i,4)
     l_sup(5,2) = ind_lim1(i,3)
     l_inf(5,1) = ind_lim2(i,2)
     l_sup(5,1) = ind_lim1(i,2)-1

     l_inf(6,3) = ind_lim1(i,6)
     l_sup(6,3) = ind_lim1(i,5)
     l_inf(6,2) = ind_lim1(i,4)
     l_sup(6,2) = ind_lim1(i,3)
     l_inf(6,1) = ind_lim1(i,1)+1
     l_sup(6,1) = ind_lim2(i,1)
     

     xg_il = dx_loc*(INT( xg(ind_grid(i),1)/dx_loc) + 0.5_dp)
     yg_il = dx_loc*(INT( xg(ind_grid(i),2)/dx_loc) + 0.5_dp)
     zg_il = dx_loc*(INT( xg(ind_grid(i),3)/dx_loc) + 0.5_dp)

     do reg=1,6

        do inz=l_inf(reg,3), l_sup(reg,3)
           iz = inz + 6                   ! +6 : respect to the center of the cubic shell 
           xpart2(3)= zg_il + dx_loc*inz 
!           if(xpart2(3) .NE. dx_loc*(INT( xg(ind_grid(i),3)/dx_loc) +0.5_dp + inz)) write(*,*) 'WARNING : zpart_il=',xpart2(3), 'zpart', dx_loc*(INT( xg(ind_grid(i),3)/dx_loc) +0.5_dp + inz) 

           do iny=l_inf(reg,2), l_sup(reg,2)    
              iy = iny + 6
              xpart2(2)= yg_il + dx_loc*iny
!              if(xpart2(2) .ne. dx_loc*(INT( xg(ind_grid(i),2)/dx_loc) +0.5_dp + iny)) write(*,*) 'WARNING : ypart_il=',xpart2(2), 'ypart', dx_loc*(INT( xg(ind_grid(i),2)/dx_loc) +0.5_dp + iny)

              do inx=l_inf(reg,1), l_sup(reg,1) 
                 ix = inx +6
                 xpart2(1)= xg_il + dx_loc*inx
!                 if(xpart2(1) .ne. dx_loc*(INT(xg(ind_grid(i),1)/dx_loc) + 0.5_dp + inx)) write(*,*) 'WARNING : xpart_il=',xpart2(1), 'xpart', dx_loc*(INT(xg(ind_grid(i),1)/dx_loc) + 0.5_dp + inx)

                 !+++  04/02/2013  ++++++++++++++++++++
                 ! here we obtain the direction to the center of the cell in the shell
                 
                 
                 mn = dirMN_ext(ind_oct, ix, iy, iz)
                 
!                 if(mn .GT. 12) write(*,*) 'WARNING: mn', mn

                 if(Mdx_ext_logical(ind_oct, ix, iy, iz, mn)) then
                    
                    call get_cell_index2(cell_ind2,cell_levl2,xpart2,il)
                    
                    if (cell_ind2 .NE. -1) then  
                       m = dirM_ext(ind_oct, ix, iy, iz)
                       n = dirN_ext(ind_oct, ix, iy, iz)
                       
                       !! and then we iterate around this direction.
                       do mloop = max(1, m-deltam), min(NdirExt_m, m+deltam)
                          do nloop = -deltan1(mloop) -1, deltan2(mloop) -1
                             nl = 1+ mod(n+ nloop+ NdirExt_n, NdirExt_n)       !cyclic value
                             mn = (mloop -1)*NdirExt_n + nl                                 
                             dx_cross_ext = dx_loc*Mdx_ext(ind_oct, ix, iy, iz, mn)   !!! ajouter true/false
                             column_dens(i,mloop,nl) = column_dens(i,mloop,nl) + dx_cross_ext*uold(cell_ind2,1)   
!                             column_dens(i,m,n) = column_dens(i,m,n) + dx_cross_ext*uold(cell_ind2,1)   
#if NSCHEM != 0
                             !if(myid .EQ. 1) write(*,*) "***VAL: Calculating H2column_dens, neulS+1=", neulS+1, "nH2=", uold(cell_ind2,neulS+1)
                             H2column_dens(i,mloop,nl) = H2column_dens(i,mloop,nl) + dx_cross_ext*uold(cell_ind2,neulS+1)
                             if(isnan(H2column_dens(i,mloop,nl))) write(*,*) "WARNING: CONT",uold(cell_ind2,neulS+1), Mdx_ext(ind_oct, ix, iy, iz, mn), dx_loc, mloop, nloop, nl, mn, m, n, ind_oct
#endif
                          end do
                       end do
                    end if                                               ! cell_ind2 .ne. -1
                 end if
                 !++++++++++++++++++++++++++++++++++++++
              end do                 !x
           end do                    !y
        end do                       !z
     end do                          !reg
     
  end do                             !i
  
end subroutine contribution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! This routine has been modified in order to include the effect of the extinction
!! due to the presence of matter

!PH modifies the name to avoid conflict with calc_temp from previous SAm/frig version
!which does not take into account extinction 
!19/01/2017
subroutine  calc_temp_extinc(NN,TT,dt_tot_unicode,vcolumn_dens,coeff_chi) 
  use amr_parameters
  use hydro_commons
  
  implicit none
  
  real(dp)                                                   :: NN,TT, dt_tot_unicode,coeff_chi
  real(dp), dimension(1:NdirExt_m,1:NdirExt_n),intent(inout) :: vcolumn_dens
  
  integer             :: n,i,j,k,idim, iter, itermax               
  real(dp)            :: dt, dt_tot, temps, dt_max, itermoy,extinct
  real(dp)            :: rho,temp
  real(dp)            :: mm,uma, kb, alpha1,mu,kb_mm
  real(dp)            :: TTold, ref,dRefdT, eps, vardt,varrel, dTemp
  real(dp)            :: rhoutot2
  real(dp)            :: scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2
  
  
  ! Cette routine fonctionne en cgs
  ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
  kb  =  1.38062d-16   ! erg/degre
  !  uma =  1.660531e-24  ! gramme
  !  mu  =  1.4
  !  mm = mu*uma 
  !  kb_mm = kb / mm 
  !  TT = TT  / kb  !/ kb_mm 
  
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  
  if( TT .le. 0.) then 
     TT = 50. / scale_T2
     return
  endif
  
  !use for the shell problem to keep high the temperature of the coronal gas 
  if( TT * scale_T2 .ge. 20000. ) then 
     return
  endif
  
  
  !if( TT*scale_T2 .gt. 50.) then
  !TT = 50. / scale_T2
  !return
  !endif
  
  
  vardt = 10.**(1./10.); varrel = 0.2
  
  
  
  dt_tot = dt_tot_unicode * scale_t ! * 3.08d18 / sqrt(kb_mm)
  TT     = TT * scale_T2
  
  
  !  nn = (rho/(gramme/cm3)) /mm
  
  itermax = 0 ; itermoy = 0.
  
  
  
  if (NN .le. smallr) then 
     if( NN .le. 0)  write(*,*) 'prob dens',NN
     NN = smallr  !max(NN,smallr)
  endif
  
  
  !     alpha1 = NN*kb_mm/(gamma-1.)
  alpha1 = NN*kb/(gamma-1.)
  
  iter = 0 ; temps = 0.
  do while ( temps < dt_tot)
     
     
     if (TT .lt.0) then
        write(*,*) 'prob Temp',TT, NN
        !         write(*,*) 'repair assuming isobariticity'
        NN = max(NN,smallr)
        TT = min(4000./NN,8000.)  !2.*4000. / NN 
     endif
     
     
     TTold = TT       
     
     !NN is assumed to be in cc and TT in Kelvin
     
     
     !! here we pass the value of the column density
     
     call chaud_froid_2(TT,NN,ref,dRefdT,vcolumn_dens,scale_l,coeff_chi)             
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     
     if (iter == 0) then
        if (dRefDT .ne. 0.) then 
           dt = abs(1.0E-1 * alpha1/dRefDT) 
        else 
           dt = 1.0E-1 * dt_tot 
        endif
        dt_max = dt_tot - temps
        if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
     endif
     
     dTemp = ref/(alpha1/dt - dRefdT) 
     
     eps = abs(dTemp/TT)
     if (eps > 0.2) dTemp = 0.2*TTold*dTemp/abs(dTemp)
     
     TT = TTold + dTemp
     if (TT < 0.) then
        write(*,*) 'Temperature negative !!!'
        write(*,*) 'TTold,TT   = ',TTold,TT
        write(*,*) 'rho   = ',rho 
        TT = 100.  !*kelvin
     endif
     
     
     iter = iter + 1
     
     temps = temps + dt
     
     dt = vardt*varrel*dt/Max(vardt*eps, varrel)
     
     dt_max = dt_tot - temps
     if (dt > 0.7*dt_max) dt = dt_max*(1.+1.0E-12)
     !        write(*,987) temps, TT
     !987     format(E10.3,2x,E10.3)
     !        read (*,*)
  enddo
  
  
  !  if (TT .ge. 50.)  TT=50.
  
  !!now convert temperature in code units
  TT = TT / scale_T2
  
  
  return
end subroutine calc_temp_extinc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2012 VV
!! This routine gives the value of the distance crossed inside a cell
!! centered at xcel, of side dx0 in the direction m,n 
!! defined with respect to x0

subroutine get_dx(x0,xcel,m,n,dx0,dx_cross)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none

  real(dp),dimension(1:3), intent (in)  :: x0, xcel
  integer, intent (in)                  :: m,n
  real(dp), intent (in)                 :: dx0
  real(dp), intent (out)                :: dx_cross

  real(dp), dimension(3,2)              :: xplane
  real(dp), dimension(2,3)              :: sol                  !2 solutions et 3 coordonnees par solution
  real(dp), dimension(3)                :: xx

  real(dp)                              :: dx_2,xal1, xal2, solaux
  real(dp), dimension(3)                :: delx
  integer                               :: nfound, i, j, indnf, mod1,mod2


  nfound  = 0
  sol(:,:)= 0.0_dp
  dx_2    = dx0/2.0_dp

  do i=1, 3
     delx(i)     = xcel(i) - x0(i)     ! position of the crossed cell relative to x0(the target)
     xplane(i,1) = delx(i) + dx_2      ! faces of crossed cell 
     xplane(i,2) = delx(i) - dx_2
  end do
    
  ! we write the equations in a cyclic manner in order to have
  ! the same structure. The equation system comes from:
  ! ycos(theta) - zsin(theta)sin(phi) = 0
  ! xcos(theta) - zsin(theta)cos(phi) = 0
  ! xsin(theta)sin(phi) - ysin(theta)cos(phi) = 0
 
  iloop: do i=1, 3
     mod1 = mod13(i)              !just the value of i+1 in modulo 3 (cyclic values)
     mod2 = mod23(i)
     xal1 = xalpha(m,n,i,1)
     xal2 = xalpha(m,n,i,2)

     do j=1, 2                   !Loop for both planes (+/-dx/2)
        xx(i)    = xplane(i,j)   !we set the position in the plane of one of the cell faces 
        xx(mod1) = xal1*xx(i)    !we look for the intersection with the other planes.
        xx(mod2) = xal2*xx(i)
        
        !! we look for any of the 6 planes if the line that passes through the 
        !! point x0 in the direction defined by m,n intersects the cell in the shell.   
        if( (abs( xx(mod1) - delx(mod1) ) .LE. dx_2) .AND. (abs(xx(mod2) - delx(mod2) ) .LE. dx_2) ) then 
           nfound = nfound + 1
           do indnf=1, 3
              sol(nfound, indnf) = xx(indnf)
           end do
        end if
        if(nfound .EQ. 2) EXIT iloop   ! if we have found 2 points we stop
 
     end do
  end do iloop
   
  dx_cross = 0.0_dp

  !! if we find two intersection points we calculate the length of the 
  !! segment crossed through the cell in the shell
  if(nfound .EQ. 2) then
     do i=1, 3
        solaux = (sol(1,i)-sol(2,i))
        dx_cross = dx_cross + solaux*solaux
     end do
  end if
  dx_cross = sqrt(dx_cross)

end subroutine get_dx


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2012 VV
!! This routine gives the closest direction m,n to the direction from x0 to x1

subroutine get_mn(x0,x1,m,n)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none

  real(kind= dp),dimension(1:3), intent (in)  :: x0, x1
  integer, intent (out)                       :: m,n
  real(kind= dp)                              :: rr, r, phi, cos_theta, cos_phi, sin_phi

  !! we define the values of the spherical coordinates

  rr= (x1(1)-x0(1))**2 + (x1(2)-x0(2))**2  
  r= rr+(x1(3)-x0(3))**2
  rr= sqrt(rr)
  r= sqrt(r)
  cos_theta= (x1(3)-x0(3))/r

  ! the calculation of m is straightforward
  m = min(INT((1.0-cos_theta)*NdirExt_m/2.0) + 1,NdirExt_m)
  
  ! for the calculation of n we have to analyze each case
  if(rr .EQ. 0) then 
     n = 1                       ! the vertical directions are degenerated in phi, then we use 1.
   
  else                           ! otherwise it depends on the values of sin and cos
     cos_phi= (x1(1)-x0(1))/rr
     sin_phi= (x1(2)-x0(2))/rr
!     if(abs(cos_phi) .GT. 1 .OR. abs(sin_phi) .GT. 1) write(*,*) 'WARNING: forbidden values for cos, sin =', cos_phi, sin_phi
     
     if(sin_phi .GE. 0) then
        phi = acos(cos_phi)
     else
        phi = 2.*pi_g - acos(cos_phi)
     end if
     

     n =  mod(INT(phi*NdirExt_n/(2.0*pi_g)+0.5),NdirExt_n) + 1

  end if
  
end subroutine get_mn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2012 VV
!! This routine precalculates the correction factors and the value of pi (pi_g).
!! It is called from "adaptative_loop.f90" if in the namelist radiative=.true.

subroutine init_radiative
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none

  ! geometrical corrections
  real(dp)                            :: phi, cos_theta, sin_theta, cos_phi, sin_phi, cos_theta2
  real(dp)                            :: dx_factor
  real(dp), dimension(1:3)            :: xzer, xpar
  integer                             :: m, n, ix, iy, iz, mn
  integer                             :: ind, ind1, ind2, ii, i
  real(dp),dimension(1:twotondim,1:3) :: xc

  allocate(xalpha(1:NdirExt_m, 1:NdirExt_n, 1:3, 1:2) )
  allocate(Mdirection(1:twotondim,1:twotondim),Ndirection(1:twotondim,1:twotondim))
  allocate(Mdx_cross_int(1:NdirExt_m, 1:NdirExt_n), Mdx_cross_loc(1:twotondim,1:twotondim, 1:NdirExt_m, 1:NdirExt_n))
  allocate(Mdx_ext(1:73, 1:11, 1:11, 1:11, 1:(NdirExt_m*NdirExt_n) ))
  allocate(Mdx_ext_logical(1:73, 1:11, 1:11, 1:11, 1:(NdirExt_m*NdirExt_n) ))
  allocate(dirM_ext(1:73, 1:11, 1:11, 1:11), dirN_ext(1:73, 1:11, 1:11, 1:11), dirMN_ext(1:73, 1:11, 1:11, 1:11) )



  Mdx_ext(:,:,:,:,:)=0
  Mdx_ext_logical(:,:,:,:,:)= .false.  
  dirM_ext(:,:,:,:) = 0
  dirN_ext(:,:,:,:) = 0
  dirMN_ext(:,:,:,:) = 0

  pi_g = 2.0*acos(0.0)    !global value of pi for cooling_fine.f90

  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5_dp)
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5_dp)
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5_dp)
  end do


  !----- xalpha for  GET_DX   -----
  do n=1, NdirExt_n
     phi       = 2.0_dp*pi_g*(n-1)/NdirExt_n
     cos_phi   =  dsign(max(abs(cos(phi)),1D-10),cos(phi))
     sin_phi   =  dsign(max(abs(sin(phi)),1D-10),sin(phi))
     
     do m=1, NdirExt_m
        cos_theta = (1.0_dp-2.0_dp*m)/NdirExt_m +1.0_dp
        cos_theta = dsign(max(abs(cos_theta),1D-10),cos_theta)
        sin_theta = sqrt(1.0_dp -cos_theta**2)

        xalpha(m,n,1,1) = sin_phi/cos_phi
        xalpha(m,n,1,2) = cos_theta/(cos_phi*sin_theta)
        xalpha(m,n,2,1) = cos_theta/(sin_phi*sin_theta)
        xalpha(m,n,2,2) = cos_phi/sin_phi
        xalpha(m,n,3,1) = cos_phi*sin_theta/cos_theta
        xalpha(m,n,3,2) = sin_phi*sin_theta/cos_theta

     end do
  end do

  ! Precalculation of function modulo 3
  do i=0, 2
     mod13(i+1)    = 1 + mod(i+1,3)
     mod23(i+1)    = 1 + mod(i+2,3)
  end do



  !-----   CALCULATION OF GEOMETRICAL CORRECTIONS   -----
  ! the correction factors are calculated for a cell of side 1
  ! the value of dx_crossed at each level can be calculated by
  ! multiplying the correction factor by dx
  
  ! Internal and Local corrections----------------------------
  do ind=1,twotondim
     xzer(1) = xc(ind,1)
     xzer(2) = xc(ind,2)
     xzer(3) = xc(ind,3)
     do n = 1,NdirExt_n
        do m = 1,NdirExt_m
           call get_dx(xzer,xzer,m,n,1.0_dp,dx_factor)        
           Mdx_cross_int(m,n) = dx_factor/2.0_dp   
        end do
     end do
     do ii=1,twotondim
        if(ii .NE. ind) then
           xpar(1) = xc(ii,1)
           xpar(2) = xc(ii,2)
           xpar(3) = xc(ii,3)
           call get_mn(xzer,xpar, m,n)
           Mdirection(ind,ii) = m
           Ndirection(ind,ii) = n
           do n=1,NdirExt_n
              do m=1,NdirExt_m  
                 call get_dx(xzer,xpar,m,n,1.0_dp,dx_factor)
                 Mdx_cross_loc(ind,ii,m,n) = dx_factor
              end do
           end do
        end if
     end do
  end do
  
  ! External corrections
  !xzer is the cell position where we want to calculate the contribution
  !to the screening caused by the surrounding cells
  !ind  1 to 64 : position of a cell in a ilevel -3 configuration
  !ind 65 to 72 :                         ilevel -2
  !ind     = 73 :                         ilevel -1

  do ind1=1, 9
     ind = 64 + ind1
     ! if(myid .EQ. 1)write(*,*) ind1, ind
     if(ind1 .EQ. 9) then
        xzer(:) = 0.0_dp      !position of the center of the oct
     else
        xzer(1) = 0.5_dp*xc(ind1,1)  !position of a cell in the oct
        xzer(2) = 0.5_dp*xc(ind1,2)  !xc=0.5 et xzer=0.25
        xzer(3) = 0.5_dp*xc(ind1,3)
     end if
     !if(myid .EQ. 1) write(*,*) ind, xzer(1),xzer(2),xzer(3)

     !Here we take into account the surrounding cells. We consider
     !+-5 cells from xzer in each direction. -6 is to avoid negative index
     do iz=1, 11
        xpar(3) = REAL(iz) - 6.0_dp        !-5,...,+5
        do iy=1, 11
           xpar(2) = REAL(iy) - 6.0_dp
           do ix=1, 11
              xpar(1) = REAL(ix) - 6.0_dp

              !+++ 04/02/2013
              call get_mn(xzer,xpar, m,n)
              dirM_ext(ind,ix,iy,iz) = m
              dirN_ext(ind,ix,iy,iz) = n
              dirMN_ext(ind,ix,iy,iz) = (m -1)*NdirExt_n + n

!              if(myid .EQ.1 .AND. n .GT.3 .AND. xpar(3) .GT. 0) write(*,*) 'WARNING: x=', xpar(1), xpar(2), xpar(3), 'i=', ix, iy, iz
!              if(dirMN_ext(ind,ix,iy,iz) .GT. 12) write(*,*) 'WARNING : ind=',ind, 'm=',m, 'n=',n, 'mn=', (m -1)*NdirExt_n + n

              do n=1, NdirExt_n
                 do m=1, NdirExt_m
                    call get_dx(xzer, xpar, m, n, 1.0_dp,dx_factor)
                    mn = (m-1)*NdirExt_n + n
                    Mdx_ext(ind,ix,iy,iz,mn) = dx_factor
                    if(dx_factor .GT. 0) Mdx_ext_logical(ind,ix,iy,iz,mn) = .true.
                    ! if (myid .EQ. 1)write(*,*)ind, ix,iy,iz,m,n, dx_factor
                 end do
              end do
           end do
        end do
     end do
  end do

  do ind1=1,8
     do ind2=1,8
        ! if(myid .EQ. 1) write(*,*) 'ind1, ind2',  ind1, ind2

        ind = (ind1-1)*8 + ind2    !ind 1 to 64

        !position of a cell in a grid of level ilevel-3
        xzer(1) = 0.5_dp*xc(ind1,1) + 0.25_dp*xc(ind2,1)
        xzer(2) = 0.5_dp*xc(ind1,2) + 0.25_dp*xc(ind2,2)
        xzer(3) = 0.5_dp*xc(ind1,3) + 0.25_dp*xc(ind2,3)
        do iz=1, 11
           xpar(3) = REAL(iz) - 6.0_dp
           do iy=1, 11
              xpar(2) = REAL(iy) - 6.0_dp
              do ix=1, 11
                 xpar(1) = REAL(ix) - 6.0_dp

                 !+++ 04/02/2013: 
                 call get_mn(xzer,xpar, m,n)
                 dirM_ext(ind,ix,iy,iz) = m
                 dirN_ext(ind,ix,iy,iz) = n
                 dirMN_ext(ind,ix,iy,iz) = (m -1)*NdirExt_n + n
 
!                 if(dirMN_ext(ind,ix,iy,iz) .GT. 12) write(*,*) 'WARNING : ind=',ind,'m=',m, 'n=',n, 'mn=', (m -1)*NdirExt_n + n                 

                 do n=1, NdirExt_n
                    do m=1, NdirExt_m
                       call get_dx(xzer, xpar, m, n, 1.0_dp,dx_factor)
                       mn = (m-1)*NdirExt_n + n
                       Mdx_ext(ind,ix,iy,iz,mn) = dx_factor
                       if(dx_factor .GT. 0) Mdx_ext_logical(ind,ix,iy,iz,mn) = .true.
                       ! if (myid .EQ. 1)write(*,*)ind, ix,iy,iz,m,n, dx_factor
                    end do   !m
                 end do      !n
              end do         !ix
           end do            !iy
        end do               !iz
     end do                  !ind2
  end do                     !ind1

end subroutine init_radiative
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!=====================================================================================================
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 2012 VV Modified to include the extinction
subroutine chaud_froid_2(T,n,ref,dRefDT,vcolumn_dens,sc_l,coeff_chi)    
  use hydro_commons,only:NdirExt_m,NdirExt_n
  use amr_parameters
  
  implicit none
  
  real(dp)                                                   :: T, n, ref, dRefDT, sc_l,coeff_chi    
  real(dp), dimension(1:NdirExt_m,1:NdirExt_n),intent(inout) :: vcolumn_dens

  integer                :: nnn,index_m,index_n          
  real(dp)               :: P,x,ne,extinct,coef             !x est le taux d'ionisation
  real(dp)               :: T2, ref2
  real(dp)               :: froid,chaud,nc, froidc    
  real(dp)               :: froidcII, froido, froidh
  real(dp)               :: froidc_m,froidcII_m,froido_m
  real(dp)               :: param, G0, epsilon,k1,k2,bet,froidrec
  real(dp)               :: eps

  !fonction de chauffage et de refroidissement calculee a partir des
  !refroidissements des differents elements



  !abondance de carbone 3.5 10-4, depletion 0.4

!!! on calcule l'ionisation 
!!! on suppose que si x est superieure a 1.d-4 elle est domine par 
!!! l'hydrogene et que la valeur inf est donne par le carbone 
!!! et vaut 3.5 1.d-4 * depletion * densite

!!! Pour les electrons dus a l'hydrogene on prend la 
!!! formule donnee par Wolfire et al. 2003 appendice C2.
!!! on prend un taux d'ionisation de 1.d-16 G0'=GO/1.7
!!! Z'd = 1 et phipah=0.5



  ne = 2.4E-3*((T/100.)**0.25_dp)/0.5_dp !formule C15 Wolfire et al. 2003 

  x = ne / N   ! ionisation

  x = min(x,0.1)

  x = max(x,3.5E-4*0.4_dp)


  !transition hyperfine a basse temperature: carbone et oxygene
  !chiffre pris dans la these de Karl Joulain 

  !refroidissement par le carbone

  !      froidcII = ( 2.2d-23                     !excitation par H
  !     c          + 5.5d-20 * 2 / sqrt(T) * x ) !excitation par e 
  !     c              * 3.5d-4 * 0.4d0 * exp(-92.d0 / T)


  froidcII =  92. * 1.38E-16 * 2. * (2.8E-7* ((T/100.)**(-0.5_dp))*x + 8.E-10*((T/100.)**(0.07_dp))) &
       * 3.5E-4 * 0.4_dp * exp(-92.0_dp/ T)
  !     c               3.d-4  * exp(-92. / T)


  !refroidissement par l'oxygene 
  !abondance 8.6 10-4 depletion 0.8

  froido = 1.E-26 * sqrt(T) * (24. * exp(-228./ T) + 7. * exp(-326./ T) )


  !      froido =  230.d0*1.38d-16 * (
  !     c            1.4d-8*x + 9.2d-11 *(T /100.d0)**(0.67) ) 
  !     c     * exp(-230.d0 / T)  

  !      froido = froido +  330.d0*1.38d-16 *( 
  !     c            1.4d-8*x + 4.3d-11 *(T /100.d0)**(0.8) ) 
  !     c     * exp(-330.d0 / T) 

  !      froido = froido +  98.d0*1.38d-16 * (
  !     c            5.d-9 *x + 1.1d-10* (T /100.d0)**(0.44) ) 
  !    c      * exp(-98.d0 / T) 


  !       froido = 2.5d-27 * (T/100)**0.4 * exp(-228.d0 / T)  


  !on tient compte de l'abondance du 
  froido = froido * 4.5E-4 


  !refroidissement par l'hydrogene 
  ! formule de Spitzer 1978
  froidh = 7.3E-19 * x * exp(-118400./ T )


  !refroidissement par les raies metastables des metaux
  !chiffre pris dans Hollenbach and McKee 1989 (ApJ 342, 306)


  !carbone une fois ionise ,1 transitions 2P 4P
  ! le poids est 1
  ! 2P->4P : 
  ! les expressions des coefficients d'excitation ont une dependance
  !en la temperature differente au dela de 10000K 
  !abondance 3.5 d-4 depletion 0.4

  !       froidcII_m = 6.2d4 * 1.38d-16 * 1.d0 *     !transition 2P->4P
  !     c ( 2.3d-8* (T/10000.)**(-0.5) * x + 1.d-12 ) *exp(-6.2d4 / T) 
  !     c    * 3.5d-4 * 0.4 




  !       if ( T .le. 1.d4 ) then 
  !       froido_m = 2.3d4 * 1.38d-16 / 3.d0 *   
  !     c ( 5.1d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.3d4/T) 
  !      
  !       froido_m = froido_m + 
  !     c       4.9d4 * 1.38d-16 / 3.d0  *   
  !     c ( 2.5d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-4.9d4/T) 
  !

  !       froido_m = froido_m + 
  !     c       2.6d4 * 1.38d-16 * 1.d0  *   
  !     c ( 5.2d-9 * (T/10000.)**(0.57) * x + 1.d-12) *exp(-2.6d4/T) 

  !       else 

  !       froido_m = 2.3d4 * 1.38d-16 / 3.d0 *    
  !     c ( 5.1d-9 * (T/10000.)**(0.17) * x + 1.d-12) *exp(-2.3d4/T) 
  !      
  !       froido_m = froido_m + 
  !     c       4.9d4 * 1.38d-16 / 3.d0  *     
  !     c ( 2.5d-9 * (T/10000.)**(0.13) * x + 1.d-12) *exp(-4.9d4/T) 


  !       froido_m = froido_m + 
  !     c       2.6d4 * 1.38d-16 * 1.d0  *     
  !     c ( 5.2d-9 * (T/10000.)**(0.15) * x + 1.d-12) *exp(-2.6d4/T) 


  !       endif

  !! abondance de l'oxygene
  !       froido_m = froido_m *   4.5d-4 



!!! on somme les refroidissements
  froid = froidcII  + froidh  + froido  !+ froido_m !+  froidcII_m 


  !      froid=froid*1.d-13    !conversion en MKS


  !refroidissement par le carbone neutre. On suppose l'equilibre
  ! de la reaction C + hv <=> C+ + e-
  ! les taux de reactions et de refroidissement sont pris dans
  !la these de Karl Joulain. 

  ! abondance du carbone relative a n (MKS)     


  !    C+ + e- => C 
  !       k1 = 4.4d-12 * (T/300.)**(-0.61) !s^-1 cm^-3

  !       k1 = k1 


  !    C => C+ + e-
  !       k2 = 2.2d-10 


  ! on a : [C] = k1/k2 [C+] * [e-]
  ! on suppose que tout le carbone est sous forme C+
  ! et que [e-] = [C+]

  ! l'abondance relative de carbone
  !      nc = k1/k2 * (3.5d-4*0.4)**2 * n


  !      froidc =  1.0d-24 * ( 1.4d0 * exp( -23.d0 / T )  
  !     c                + 3.8d0 * exp( -62.d0 / T )   )

  !      froidc = froidc * nc !(nc est l'abondance relative du carbone)


  !       n=exp(log(10.d0)*logn) !ici pas besoin de log

  !       valeur utilisees au 23/08/98
  !       chaud=4.d0*exp(-24.5d0*log(10.d0))*1.d-7  !un peu empirique ....



!!!! on calcule le chauffage
!!! on prend en compte le chauffage sur les grains 
!!! formules 1 et 2  de Wolfire et al. 1995 


  !! G0 is the UV field defined by Habing and Draine
  G0 = 1.0_dp/1.7_dp
  G0 = G0*p_UV             ! p_UV: parameter of variation of UV
                           ! defined in the namelist
  coeff_chi  = 1.0D0

  !---  we calculate the extinction using the column density   ----
  if(extinction) then
     extinct=0.0_dp

     !! Loop in theta and phi 
     do index_m=1,NdirExt_m
        do index_n=1,NdirExt_n

           ! now take the exponential and sum over directions 
           extinct = extinct + exp(-vcolumn_dens(index_m,index_n)*coef) 
        end do
     end do
     coeff_chi  = extinct/(NdirExt_m*NdirExt_n)
     G0 = G0*coeff_chi  

  end if
  !G0 has been modified and it considers the extinction
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  

  param = G0 * sqrt(T)/(n*x)
  epsilon = 4.9E-2 / (1. + (param/1925.)**0.73_dp)
  epsilon  = epsilon + 3.7E-2 * (T/1.E4)**0.7_dp / (1. + (param/5.E3) )

  chaud = 1.E-24 * epsilon 

  ! pour un flux de rayonnement G0/1.7 
  chaud = chaud * G0  
  chaud = chaud + 1.0E-27_dp  !heating including the cosmic-ray heating rate
                            !for dark cloud cores
                            !REF: (Eq 3) Goldsmith 2001


  !refroidissement recombinaison sur les grains charges positivement
  bet = 0.74/(T**0.068)
  froidrec = 4.65E-30*(T**0.94)*(param**bet)*x

  !! chaud=1.d-32 !un peu empirique ....

  !      froidc=0.d0


  ref= chaud*n - (n**2)*(froid + froidrec) !!!+ froidc)

!------------------------------------------------------

  eps = 1.0E-5
  T2 = T*(1.+eps)

  ne = 2.4E-3*((T2/100.)**0.25)/0.5_dp    ! (eqn C15) Wolfire et al. 2003 
  x = ne / N                              ! ionisation
  x = min(x,0.1)
  x = max(x,3.5E-4*0.4)


  froidcII =  92. * 1.38E-16 * 2. * (2.8E-7* ((T2/100.)**(-0.5_dp))*x + 8.E-10*((T2/100.0_dp)**(0.07_dp))) &
       * 3.5E-4 * 0.4_dp * exp(-92.0_dp/ T2)

  froido = 1.E-26 * sqrt(T2) * (24. * exp(-228./ T2) + 7. * exp(-326./ T2) )
  froido = froido * 4.5E-4 

  froidh = 7.3E-19 * x * exp(-118400./ T2 )

  froid = froidcII  + froidh  + froido  !+ froido_m !+  froidcII_m 


  param = G0 * sqrt(T2)/(n*x)
  epsilon = 4.9E-2 / (1. + (param/1925.)**0.73)
  epsilon  = epsilon + 3.7E-2 * (T2/1.E4)**0.7 / (1. + (param/5.E3) )


  chaud = 1.E-24 * epsilon 

  chaud = chaud * G0  
  chaud = chaud + 1.0E-27_dp !heating including the cosmic-ray heating rate
  
  bet = 0.74/(T2**0.068_dp)
  froidrec = 4.65E-30*(T2**0.94_dp)*(param**bet)*x

  ref2= chaud*n - (n**2)*(froid + froidrec) !!!+ froidc)

  dRefDT = (ref2-ref)/T/eps

  return 

end subroutine chaud_froid_2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!calculate the extinction
!! 21/11/2017
!!PH adapts the cooling_fine adapted by Benoit from Valeska's original routine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine extinction_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute extinction for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,info
  integer,dimension(1:nvector),save::ind_grid
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  ! files
  character(LEN=5)                    :: nsort, nocpu
  character(LEN = 80)                 :: filenamex,filenamey,filenamez
  integer::uleidx,uleidy,uleidz

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  !  if(myid .EQ. 1) write(*,*) 'STATUS: starting cooling_fine, ilevel=', ilevel
  ! Valeska
  !--------------------------------------------------------------------
  ! The write option permits us to construct a column density map projected 
  ! in x, y, and z by calculating the column densities along positive and 
  ! negative directions with respect to a slice that passes through the 
  ! respective mid planes.
  !--------------------------------------------------------------------
  if(writing) then
     call title(ifout-1, nsort)
     call title(myid, nocpu)
     filenamex= TRIM(nsort)//'_test_densX_'//TRIM(nocpu)//'.dat'        !ex.:00001_test_densX_00010.dat
     filenamey= TRIM(nsort)//'_test_densY_'//TRIM(nocpu)//'.dat'
     filenamez= TRIM(nsort)//'_test_densZ_'//TRIM(nocpu)//'.dat'
     
     uleidx = myid + 100                                                !integer
     uleidy = myid + 200
     uleidz = myid + 300
     
     open(unit=uleidx, file=filenamex, form='formatted', status='unknown',position='append')
     open(unit=uleidy, file=filenamey, form='formatted', status='unknown',position='append')
     open(unit=uleidz, file=filenamez, form='formatted', status='unknown',position='append')
  end if
  !--------------------------------------------------------------------
  ! Valeska

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call extinctionfine1(ind_grid,ngrid,ilevel)
  end do

  
111 format('   Entering extinction_fine for level',i2)

end subroutine extinction_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine extinctionfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module

  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc,ix,iy,iz,info
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant
  integer,dimension(1:nvector),save::ind_cell,ind_leaf,ind_leaf_loc
  real(kind=8),dimension(1:nvector),save::nH,T2,T2_new,delta_T2,ekk,err,emag
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  integer::irad

  real(dp) :: barotrop1D,mincolumn_dens
  real(dp)                                   :: x0, y0, z0,coeff_chi,coef
  double precision                           :: v_extinction,extinct
  integer::uleidx,uleidy,uleidz,uleidh,igrid,ii,indc2,iskip2,ind_ll

  real(dp),dimension(1:twotondim,1:3)        :: xc                                               !xc: coordinates of center/center grid

  !-------------- SPHERICAL DIRECTIONS ------------------------------------------------!
!  real(dp),dimension(1:nvector,1:ndir)                 :: col_dens                     !
  real(dp),dimension(1:nvector,1:NdirExt_m,1:NdirExt_n):: column_dens,column_dens_loc  !
  real(dp),dimension(1:nvector,1:NdirExt_m,1:NdirExt_n):: H2column_dens,H2column_dens_loc  ! H2 column density
!  real(dp),dimension(ndir)                             :: vcol_dens                    !
  real(dp),dimension(1:NdirExt_m,1:NdirExt_n)          :: vcolumn_dens
  real(dp),dimension(1:3)                              :: xpos
  real(dp)                                             :: dx_cross_int, dx_cross_loc
  integer                                              :: index_m,index_n, mmmm,nnnn
  integer                                              :: m, n, mloop, nloop, nl   
  integer, dimension(1:NdirExt_n)                      :: deltan1, deltan2
  integer                                              :: deltam
  !-----   simple_chem   --------------------------------------------------------------!
  real(dp)                                             :: xshield, fshield, kph, m_kph, b5   ! cgs                                                 
  real(dp)                                             :: small_exp=1d-10, e_valexp, small_arg=-50.0_dp
  real(dp)                                             :: cst1 = -8.5d-4, one=1.0, cst2, coeff, valexp
  real(dp)                                             :: mH2 = 3.347449d-24           ! g                                                         
  real(dp)                                             :: kbcgs  =  1.380658d-16       ! erg/degre                                                 
  real(dp)                                             :: bturb = 2d5    
     ! 2km/s (Maryvonne Gerin, Jacques Le Bourlot, Pierre Lesaffre) 1:Gnedin+2009 !7.1d5(Krumholz2012)
  real(kind=8),dimension(1:nvector),save::nH2
  !!take care "nent" should be added here. It is fine if it is 0 of course 
  integer                                               :: neulS=8+nrad+nextinct
  !------------------------------------------------------------------------------------!


  !Valeska
!  vcol_dens(:) = 0.
  column_dens(:,:,:) = 0.
  H2column_dens(:,:,:) = 0.
  vcolumn_dens(:,:) = 0.

  
  if(writing) then
     !---position of reference to calculate the column density maps---
     ! we add 1.0D-09 in order to avoid the exact center (there are not cells centered in 0.5L)
     x0 = 0.5D0 + 1.0D-09
     y0 = 0.5D0 + 1.0D-09
     z0 = 0.5D0 + 1.0D-09
     
     !---Units uleidx, uleidy, uleidz---
     uleidx = myid + 100
     uleidy = myid + 200
     uleidz = myid + 300
  end if
  
294 FORMAT(I10,4ES14.5) 
295 FORMAT(5ES14.5)   
296 FORMAT(I10,5ES14.5)                
  
  if(numbtot(1,ilevel)==0)return
  
  !get the column density within the box from the grid faces

  !-----   EXTERNAL CONTRIBUTION   -----
  call column_density(ind_grid,ngrid,ilevel,column_dens,H2column_dens) 
  !Valeska

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

  !scaling coefficent between column density and Av 
  coef = 2.d-21 * scale_l * boxlen       !cm^2; Eq.34 Glover & Mac Low 2007


  !--- Set position of cell centers relative to grid center ---
  do ind=1,twotondim
     iz=(ind-1)/4                               !0 for ind=1,2,3,4; 1 for ind=5,6,7,8
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5_dp)*dx   
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5_dp)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5_dp)*dx
  end do

!******************************************************
!--- GEOMETRICAL CORRECTIONS:  Internal and Local -----
!  do ind=1,twotondim
!     xpos(1) = xc(ind,1)
!     xpos(2) = xc(ind,2)
!     xpos(3) = xc(ind,3)
!     
!     do index_m = 1,NdirExt_m
!        do index_n = 1,NdirExt_n
!           call get_dx(xpos,xpos,index_m,index_n,dx,dx_cross_int)        
!           Mdx_cross_int2(index_m,index_n) = dx_cross_int/2.0_dp
!        end do
!     end do
!     
!     do ii=1,twotondim
!        if(ii .NE. ind) then
!           xcell(1) = xc(ii,1)
!           xcell(2) = xc(ii,2)
!           xcell(3) = xc(ii,3)
!           call get_mn(xpos,xcell, m,n)
!           Mdirection2(ind,ii) = m
!           Ndirection2(ind,ii) = n
!           
!           do index_m=1,NdirExt_m  
!              do index_n=1,NdirExt_n
!                 call get_dx(xpos,xcell,index_m,index_n,dx,dx_cross_loc)
!                 Mdx_cross_loc2(ind,ii,index_m,index_n) = dx_cross_loc
!              end do
!           end do
!        end if
!     end do
!
!  end do
!----------------------------------------------------
!******************************************************


  !---  PRECALCULATION of  m, n loop limits  ------------------------------------------
  ! each cell can contribute to the column density in several directions, then
  ! we define here the limits for the loop around the direction to the cell center.
  ! For the vertical directions (m =1, NdirExt_m), the azimuthal angle has to cover 2*pi,
  ! for other directions we use +/- 1/8 of the total directions 
  !-------------------------------------------------------------------------------------
     do m = 1, NdirExt_m
        if(m .EQ. 1 .OR. m .EQ. NdirExt_m) then
           deltan1(m) = INT(NdirExt_n/2.0_dp)
           deltan2(m) = deltan1(m) -1
        else
           deltan1(m) = INT(NdirExt_n/8)
           deltan2(m) = deltan1(m)
        end if
     end do
     deltam = INT((NdirExt_m-1)/4.)


  ! Loop over cells
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
           ind_leaf_loc(nleaf)=i       !index within local group of cells
        end if
     end do
     if(nleaf.eq.0)cycle

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do


#if NEXTINCT>0
     do i=1,nleaf
        nH2(i)=MAX(uold(ind_leaf(i),neulS+1),smallr)
     end do
#endif


#if NEXTINCT>1
     do i=1,nleaf                                  !loop over leaf cells 
        ind_ll=ind_leaf_loc(i)
        
        igrid=ind_grid(ind_leaf_loc(i))            !index father  
        
        xpos(1) = xg(igrid,1) + xc(ind,1)          !grid position + leaf position relative to grid center
        xpos(2) = xg(igrid,2) + xc(ind,2)
        xpos(3) = xg(igrid,3) + xc(ind,3)
               
                   
           !------     +  INTERNAL CONTRIBUTION        ------
           ! Here we sum up the contribution due to the cell itself. Its 1/2 in each direction.
           
           do index_n=1,NdirExt_n      !loop over directions to compute the screaning of half the cells on itself 
              do index_m=1,NdirExt_m 
                 
                 column_dens_loc(ind_ll,index_m,index_n) = column_dens(ind_ll,index_m,index_n) + dx*Mdx_cross_int(index_m,index_n)*nH(i)
#if NEXTINCT>1
                 H2column_dens_loc(ind_ll,index_m,index_n) = H2column_dens(ind_ll,index_m,index_n) + dx*Mdx_cross_int(index_m,index_n)*nH2(i)
                 !if(isnan(H2column_dens_loc(ind_ll,mloop,nl))) write(*,*) "WARNING: INT",H2column_dens(ind_ll,index_m,index_n), nH2(i), Mdx_cross_int(index_m,index_n), dx, index_m, index_n
#endif
                 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~         
                 ! column_dens_loc(ind_ll,index_m,index_n) = dx*Mdx_cross_int(index_m,index_n)*nH(i) !if just internal contribution is needed
                 ! column_dens_loc(ind_ll,index_m,index_n) = column_dens(ind_ll,index_m,index_n)  !if just the external component is needed
                 !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              end do
           end do
           
           !-------     + LOCAL CONTRIBUTION     -----
           ! and here we add the contributon due to its siblings (in the same oct).
           
           do ii=1,twotondim    !1-8 the cell in the oct! 
              if(ii .NE. ind) then
                 
                 iskip2  = ncoarse+(ii-1)*ngridmax 
                 indc2   = iskip2+ ind_grid(ind_ll) !indice for the cell crossed along the direction ind_dir
                 
                 !--------------------------------------------------  
                 !  !we calculate the position of the cell crossed
                 !  xcell(1)= xg(igrid,1) + xc(ii,1)
                 !  xcell(2)= xg(igrid,2) + xc(ii,2)
                 !  xcell(3)= xg(igrid,3) + xc(ii,3)
                 !  call get_mn(xpos, xcell, m, n)
                 !  if(m .NE. Mdirection(ind,ii) .OR. n .NE. Ndirection(ind,ii)) write(*,*) ind,ii
                 !  if(myid .EQ.1 .AND. ii .EQ. 8  ) write(*,*) ind, i, ii, "     m,n:", m, n
                 !-------------------------------------------------
                 

                 ! knowing the index of the target cell and the sibling cell we can find 
                 ! the closest direction from the target cell center to the sibling cell center.
                
!                 if(isnan(uold(indc2,neulS+1))) write(*,*) "WARNING LOC: nH2 is NaN: ", xpos/dx    !i, ii, uold(indc2,1), uold(indc2,neulS+1), xpos
                 m = Mdirection(ind,ii)
                 n = Ndirection(ind,ii)
                 
                 ! Here we make a loop around the direction m,n in order to treat all the 
                 ! concerned directions by the sibling cell
                 
                 do mloop = max(1, m-deltam), min(NdirExt_m, m+deltam)  !avoid forbidden intervals
                    do nloop = -deltan1(mloop) -1, deltan2(mloop) -1
                       ! the value of nloop is cyclic
                       nl = 1+ mod(n+ nloop+ NdirExt_n, NdirExt_n)
                       
                       ! Here we sum up the contribution to the column density. The distance crossed
                       ! through the cell can be found using the corrective factor Mdx_cross_loc, that
                       ! depends on the relative positions in the oct and on the direction
                       
                       column_dens_loc(ind_ll,mloop,nl) = column_dens_loc(ind_ll,mloop,nl) + dx*Mdx_cross_loc(ind,ii,mloop,nl)*uold(indc2,1)        
#if NEXTINCT>1

                       H2column_dens_loc(ind_ll,mloop,nl) = H2column_dens_loc(ind_ll,mloop,nl) + dx*Mdx_cross_loc(ind,ii,mloop,nl)*uold(indc2,neulS+1)  

!                       if(isnan(uold(indc2,neulS+1))) write(*,*) "WARNING LOC: nH2 is NaN: ", uold(indc2,1), uold(indc2,neulS+1), Mdx_cross_loc(ind,ii,mloop,nl) 
!                       if(Mdx_cross_loc(ind,ii,mloop,nl) .NE. 0.0) H2column_dens_loc(ind_ll,mloop,nl) = H2column_dens_loc(ind_ll,mloop,nl) + dx*Mdx_cross_loc(ind,ii,mloop,nl)*uold(indc2,neulS+1)  
!                       if(Mdx_cross_loc(ind,ii,mloop,nl) .NE. 0.0) H2column_dens_loc(ind_ll,mloop,nl) = H2column_dens_loc(ind_ll,mloop,nl) + dx*Mdx_cross_loc(ind,ii,mloop,nl)*uold(indc2,1)        
!                       if(isnan(H2column_dens_loc(ind_ll,mloop,nl))) write(*,*) "WARNING: LOC", uold(indc2,neulS+1), Mdx_cross_loc(ind,ii,mloop,nl), mloop, nloop, nl, dx, ind, ii
!                       if(isnan(uold(indc2,neulS+1)) .AND. Mdx_cross_loc(ind,ii,mloop,nl) .NE. 0.0) write(*,*) "WARNING LOC: nH2 is NaN and Mdx ne 0 !!!", uold(indc2,1), uold(indc2,neulS+1), Mdx_cross_loc(ind,ii,mloop,nl) 
#endif
                    end do
                 end do
                 
              end if      !ii ne ind
           end do         !ii
  
           !--- HERE THE VALUE IN COLUMN_DENS_LOC IS THE TOTAL VALUE ---

           !-- WRITE FILES for a test ---        
           if(writing .and. mod(nstep_coarse,foutput)==0) then
              
              ! we calculate here the value of the extinction just to write the files. This value is not used here and its calculated after.

              v_extinction=0.
              do index_m=1,NdirExt_m
                 do index_n=1,NdirExt_n
                    v_extinction= v_extinction+ exp(-column_dens_loc(ind_ll,index_m,index_n)*coef)
                 end do
              end do
              v_extinction= v_extinction/(NdirExt_m*NdirExt_n)
              
              ! we calculate the closest directions to the +/- cartesian directions 
              mmmm=NINT((NdirExt_m-1.)/2.)+1     ! (5-1)/2 +1 = 3 ok
              nnnn=NINT(NdirExt_n/2.)+1          !     8/2 +1 = 5 ok  

              if(abs(xpos(1)-x0) .LT. 0.5*dx) write(uleidx,296) ilevel, xpos(2), xpos(3), column_dens_loc(ind_ll,mmmm,1), column_dens_loc(ind_ll,mmmm,NdirExt_n/2+1), v_extinction   !Write column density
              if(abs(xpos(2)-y0) .LT. 0.5*dx) write(uleidy,296) ilevel, xpos(1), xpos(3), column_dens_loc(ind_ll,mmmm,NdirExt_n/4+1), column_dens_loc(ind_ll,mmmm,3*NdirExt_n/4+1), v_extinction
              if(abs(xpos(3)-z0) .LT. 0.5*dx) write(uleidz,296) ilevel, xpos(1), xpos(2), column_dens_loc(ind_ll,1,1), column_dens_loc(ind_ll,NdirExt_m,nnnn), v_extinction 
              
           end if

           do index_m=1,NdirExt_m
              do index_n=1,NdirExt_n
                 vcolumn_dens(index_m,index_n)=column_dens_loc(ind_ll,index_m,index_n)          
              end do
           end do

           !---  we calculate the extinction using the column density   ----
           extinct=0.0_dp
           
           !! Loop in theta and phi 
           do index_m=1,NdirExt_m
              do index_n=1,NdirExt_n
                 ! now take the exponential and sum over directions 
                 extinct = extinct + exp(-vcolumn_dens(index_m,index_n)*coef) 
              end do
           end do
           coeff_chi  = extinct/(NdirExt_m*NdirExt_n)

           !store extinction in variable firstindex_extinct+1
           uold(ind_leaf(i),firstindex_extinct+1) = coeff_chi

#if NEXTINCT>1
           cst2 = scale_l* boxlen
           coeff = -2.d-21*cst2  !-1.3d-21*cst2   
           ind_ll=ind_leaf_loc(i)
           kph = 0.0
           m_kph = 0.0
           do index_m=1,NdirExt_m
              do index_n=1,NdirExt_n
                 xshield = H2column_dens_loc(ind_ll,index_m,index_n)*cst2/(5d14)   !Glover & MacLow 2007, Eq (32)                                  
                 if(isnan(xshield)) write(*,*) "WARNING: xshield is NaN", index_m, index_n, H2column_dens_loc(ind_ll,index_m,index_n)
                 !--------------------------------------------------------                                                                         
                 !!!-----   Now it is not possible to use bthermal   -----                                                                         
                 !!---- Doppler parameter using bturb and bthermal--------                                                                         
                 !!b5 = (2d0*kbcgs*Temp/mH2)                                                                                                       
                 !!b5 = b5 + bturb**2                                                                                                              
                 !!b5 = 1d-5*sqrt(b5)                                                                                                              
                 !---- using only bturb ---------------------------------                                                                          
                 b5 = 1d-5*bturb
                 !-------------------------------------------------------                                                                          
                 valexp = max(cst1*(one + xshield)**0.5, small_arg)
                 e_valexp = exp(valexp)
                 if(isnan(exp(valexp))) write(*,*) "WARNING: exp(valexp) is nan: valexp=", valexp
!                 if(exp(valexp) .LT. small_exp) write(*,*) "WARNING: exp(valexp).LT. small_exp",exp(valexp), e_valexp                             
                 fshield = 0.035*e_valexp / (one + xshield)**0.5
                 fshield = fshield + 0.965/(one + xshield/b5)**2.0
                 valexp = max(coeff*column_dens_loc(ind_ll,index_m,index_n), small_arg)
                 kph = exp(valexp) !coeff*column_dens_loc(ind_ll,index_m,index_n))                                                                 
                 kph = kph*fshield                    !kph = fshield(N_H2)*exp(-tau_{d,1000}) *kph0                                                
                 m_kph = m_kph + kph                  !mean kph                                                                                    
              end do
           end do
           uold(ind_leaf(i),firstindex_extinct+2) = m_kph/(NdirExt_m*NdirExt_n)
!           kdest = kph0*uold(ind_leaf(i),neulS+2)
#endif



     end do
#endif


  end do
  ! End loop over cells


end subroutine extinctionfine1
!###########################################################                                          
!###########################################################                                          
!###########################################################                                          
!########################################################### 

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
!!  use chemo_parameters
!!  use thermochemistry  , verbose_thermo=>verbose, mype_thermo=>mype
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
  integer                      :: neulS=8+nrad+nextinct, neulP=5
  real(dp)                     :: testf = 1.d-8

  logical:: once


  once=.true.

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
        
!#if NSCHEM !=0                                                                                                                                                        

#if NEXTINCT>1
                                                                                                                                                       
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
        if(isnan(uold(ind_leaf(i),firstindex_extinct+2))) write(*,*) "WARNING (simple_chemostep1 BEFORE solve) uold( ,neulS+2) is nan, i, ind_leaf(i)",i, ind_leaf(i) 
        k2 = kph0 * uold(ind_leaf(i),firstindex_extinct+2) ! kph_0*<fshield*exp(-tau_(d,1000) )> (eq 31, Glover&MacLow2006) 
        k1 = k1_0*sqrt(T2(i)/100.0_dp)/(1.0_dp + (T2(i)/464.0_dp)**1.5_dp)        ! k1 = k1_0*sqrt(T/100K)*S(T)

!        if(isnan(k2))  write(*,*) "WARNING (simple_chemostep1) k2 is NaN: uold=", uold(ind_leaf(i),neulS+2)
!        k2 = kph0
                                             
        ! 2nd: I solve the evolution of H2 for the next step
!        if(myid .EQ. 1) write(*,*) "CALLING SOLVE2_H2FORM"
        call solve2_H2form(nH2, ntot, k1, k2, dt_ilev)   !!! xH and xH2 must be in/out variables.

!        call solve2_H2form(nH2, ntot, k1_0, 0., dt_ilev)   !!! xH and xH2 must be in/out variables.


!        if(once) then
!           write(*,*) nh2,ntot,k1_0,dt_ilev*scale_t,t*scale_t

!           once=.false.
!        endif


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

