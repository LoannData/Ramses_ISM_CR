!================================================================
!================================================================
!================================================================
! This routine generates initial conditions for RAMSES.
! Positions are in user units:
!   x(i,1:3) are in [0,boxlen]**ndim.
! U is the conservative variable vector. Conventions are here:
!   U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E, U(i,>5): d.scalar
! Q is the primitive variable vector. Conventions are here:
!   Q(i,1): d, Q(i,2:4):u,v,w, Q(i,5): P, Q(i,>5): scalar
! If nvar > 5, remaining variables (6:nvar) are treated as passive
! scalars in the hydro solver.
!   U(:,:) and Q(:,:) are in user units.
!
!================================================================
!================================================================
!================================================================
subroutine condinit(x,u,dx,nn,first_lmax)
  use amr_parameters
  use amr_commons
  use hydro_parameters
  use poisson_parameters
  implicit none

  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:nvar+3)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  logical ::first_lmax

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  integer::ivar,i,ii,jj,kk,i_rad
  real(dp),dimension(1:nvector,1:nvar+3),save::q   ! Primitive variables
  real(dp)::x0,y0,z0,xc0,yc0,zc0
  real(dp)::xc,yc,zc,xl,yl,zl,xr,yr,zr,R0,Rin,Rout,Rout2,Hout
  real(dp)::V0,P0,P1,pi,f,rr,abar,alpha,beta
  real(dp)::m,power,pm1,pp1,thetamax,GM,GM_code,All,Alr,Arl,Arr,B0,By,Bz
  real(dp)::theta,g,D0,D1,dens,T0,r_cyl,v_cyl,r_sph,exp_arg,Rin0,Rout0,hR_clip
  real(dp)::pmin,R_integral,R_left,R_right,tmini,temp_zero,temp_minus1
  real(dp)::dR,dz,r_loop,z_loop,rho,rho_l,rho_r,temp,temp_l,temp_r,pres,pres_mid, &
             f_cent_out,f_grav_inR,f_grav_inZ,R200,R1,H1,c,Reff,p_tot,pm_tot,dens_tot
  real(dp)::Rmin_grav,Rmax_grav,dx_max,pres_l,pres_r,dR_scale,dZ_scale,dPdR,dPdz
  real(dp)::r_cyl_l,r_cyl_r,r_loop_l,r_loop_r,f_grav_inZ_l,f_grav_inZ_r,dPdz_R
  real(dp)::dPdz_R_gal,sech,sech_term, dPdR2_z0,dPdR2_z,dPdz2_R,pres_int,pres_intm1
  real(dp)::pres_p1,pres_l_p1,pres_r_p1,P_fact_cen0,P_fact_cen,soft,Tmaxi
  real(dp)::v_turb(3),mom_tot(2),ztrans1,ztrans2,rtrans1,dztrans,dz_here,ftrans1,max_z, &
            maxz_fact,maxz_fact2,maxz_fact0,alpha2,aspect_ratio
  real(dp)::exp_arg_l,exp_arg_r,r_bound,eps,dPdR_z0,rom_int,f_table,fgrav_R
  real(dp)::k5,k4,k3,delta_trans,theta_trans,eta_trans,dx_ss
  real(dp)::rcyl_ratio,rsph_ratio,z1,rsph1,pres1,pres0,rho1,rho0,rho10,vcyl2,delta_pres
  real(dp),dimension(20)::func_params
  integer::nR_loop,iR,nZ_loop,iZ,iiz,iter,iter_tot,max_iter,i_table
  real(dp)::C0,C1,C2,C3,Q1,Q2,Q0,pres_zero,pres_minus1,xscl,r_sph_l,r_sph_r
  real(dp) :: sum_dust
#if NDUST>0
  integer :: idust
  real(dp):: epsilon_0
  real(dp),dimension(1:ndust):: dustMRN
#endif

  real(dp)::Mass_of_R,MPrime_of_R,Rho_of_R,RhoPrime_of_R,Pres_of_R,PPrime_of_R,Z_of_R,ZPrime_of_R
  external::Mass_of_R,MPrime_of_R,Rho_of_R,RhoPrime_of_R,Pres_of_R,PPrime_of_R,Z_of_R,ZPrime_of_R
  external :: dPdR_z0, dPdz_R, fgrav_R, dPdz_R_gal, sech, dPdR2_z0,dPdR2_z,dPdz2_R

  pi=ACOS(-1.0d0)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  epsilon_0 = dust_ratio(1)

  dx_max = boxlen/2.d0**nlevelmax  ! max here means for maximum refinement level
  dR_scale = 1.d0
  dZ_scale = 40.d0
  eps = 1.d-9

! Gressel values ... simpler disk.  --> type ####1
! gamma = 7./5.
! R0 = 5 AU
! T0 = 108 K
! abar = 2.024 --> Chosen to Make T0 = Pres / ( rho R_idealGas abar ) = 108 K @ R0
! Sigma = 150 g / cm^2 @ R0
!  This gives D0 = Sigma / (2 * R0 * .0427878) ! I have double checked this!
!     .0427878 = integral of (exp(400(1/sqrt(1+x^2) - 1))) from 0 to .05.
!       --> D0 = 2.3434e-11 (only for unflared model, q = -1)
! theta_max = 0.4 radians above / below midplane. (only for unflared model, q = -1)
! R on = 0.5 - 5.5 AU
! res is 1024 x 384
!  v_theta = sqrt(G Msun / R ) * sqrt( (p+q)/400 + (1+q) - qR/r)
!    p = -3/2   q = -1 , gives: (note in this code I use alpha and beta instead of p and q)
!     v_theta = sqrt(G Msun / R) * r/R * sqrt( max(R/r - 1/240, 0) )
!  aspect ratio = H/R = 0.05.
!  beta_p0 = 2p/B^2 = 1e5 in midplane. I add this.
!   Bz = sqrt(beta_p0/2/p)
!

! Bai values for below: --> type 2
!   alpha = 2.0
!   D0 = 1
!   T0 = 1
!   GM = 1
!   m = 0.5
!   beta0 = 1e4 = q5/(B^2/8 pi) --> B0^2 + Bz^2 = ?

! Setting up an ideal MHD wedge version of Bai & Stone 2017
! Note, best way to get B field is use vector potential on
! corners and take first order finite difference. This well
! have left / right states be the same on neighbouring cells
!

! No B, Keplarian disk --> type 3

! Saving parameters to local variables to tidy up their use below
  Rin  =disk_innerR
  Rin0 =disk_innerR0
  Rout =disk_outerR
  Rout0=disk_outerR0
  Rout2=disk_outerR2 ! This is a boundary condition relevant R, see HTL2017, where beyond R>0.9 d and v are set.

  R0=disk_R0
  R1=disk_R1 ! scale radius
  H1=disk_H1 ! scale height

  Hout=Rout*H1/R1
  R200=disk_R200 ! virial radius, applicable for virial halo sims
  c=disk_conc ! halo concentration

  ! alpha = rho power law index ; beta = T power law index
  alpha = disk_rho_pow
  beta  = disk_cs_pow
  abar  = disk_abar
  power = (1.d0 - alpha)/2.d0 ! this is the flaring index.
  pm1 = power - 1.d0
  pp1 = power + 1.d0
  m = disk_m ! for Bai

  ! disk_xc etc in unis of dx. How many / what fraction of dx offset to I want to put the center. Usually 0 but I've looked at dx/2.
  x0 = boxlen/2.d0 + disk_xc*dx_max
  y0 = boxlen/2.d0 + disk_yc*dx_max
  z0 = boxlen/2.d0 + disk_zc*dx_max

  B0 = disk_B0
  By = disk_By
  Bz = disk_Bz

  ! Values
  V0 = disk_V0
  D0 = disk_D0
  D1 = disk_D1
  T0 = disk_T0
  P0 = 8.3144598d7 * T0 * (D0/scale_d) / abar / (scale_l**2 / scale_t**2)
  P1 = 8.3144598d7 * T0 * (D1/scale_d) / abar / (scale_l**2 / scale_t**2)
  pmin = P1/1000.d3 ! arbitrarily setting lower pressure, factor 1000

  Tmini = pmin/(D1/scale_d) ! Here and below I call T = P/rho; it's really cs^2

  thetaMax=disk_thetaMax
  hR_clip = 4.6054 ! means at boundary you are 1./100. that of clipping point

  GM = disk_Mass ! Msun
  GM = GM*6.674d-8*1.98855d33 ! cgs
  GM_code = GM/scale_l**3*scale_t**2 ! GM in code units

 
    ! Here I assume the sun has a softening out to Rin, where within Rin the mass is distributed so v has continuous 2nd derivatives.
    !  HSE solved assuming this, with Vrot = Vrot on boundary of Sun and going to zero at centre.

    ! Similar, I interpolate P and rho to be constant third derivatives. This puts limit on how much the dens and pres can increase.

    ! Polynomials for r <= R are:
    !   v(r) = v(R) * [  (99/8)(r/R)^3  -  (77/4)(r/R)^4  + (63/8)(r/R)^5  ] == v(R) * P5(r/R)
    !     --> same as if M(r<R) = Mtot * (r/R) * P5(r/R)^2
    !   rho & pres = rho(R) & pres(R) * [ Delta - P7(r/R;Delta,alpha) ]
    !          --> rho : Delta = 2 ; alpha = -1.5
    !               --> P7(r/R;2,-1.5) = (31/16)(r/R)^4 + (69/16)(r/R)^5 - (143/16)(r/R)^6 + (59/16)(r/R)^7
    !          --> Pres : Delta = 3 ; alpha = -2.5
    !               --> P7(r/R;3,-2.5) = (275/48)(r/R)^4 + (87/16)(r/R)^5 - (265/16)(r/R)^6 + (355/48)(r/R)^7

    ! 26001 --> Still setting P = rho/R GM/alpha2, but subsampling within cell to get more average value.
    !  alpha2 is now a function of R that allows it to grow / shrink, set with disk_delta_alpha.

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

    if (disk_T1 .lt. 0.d0) then
      Tmaxi = abs(disk_T1)*(Rin/R0)**beta*T0
    else
      Tmaxi = disk_T1
    endif

    if (Rmax_grav .eq. 0.d0) Rmax_grav = 1000.d0*boxlen

    P_fact_cen0 = ( GM_code )*disk_aspect**2

    if (D1 .gt. 0.d0) then
      disk_exp_limit = log(D0/D1) ! leaving disk_exp_limit as positive!!
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit*disk_aspect**2)**2 - 1.d0)
    else
      maxz_fact = sqrt( 1.d0/(1.d0 - disk_exp_limit*disk_aspect**2)**2 - 1.d0)
    endif

    ! Note, self-grav off
    dx_ss = dx/2.d0**disk_nsubsample
    dR = dx_ss

    do i=1,nn
        xc0=x(i,1) - x0
        yc0=x(i,2) - y0
        zc0=x(i,3) - z0

        dens_tot = 0.d0
        mom_tot = 0.d0
        p_tot = 0.d0
        pm_tot = 0.d0

        do kk=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
          zc = zc0 + kk*dx_ss/2.d0
          do jj=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
            yc = yc0 + jj*dx_ss/2.d0
            do ii=-2**disk_nsubsample+1,2**disk_nsubsample-1,2
              xc = xc0 + ii*dx_ss/2.d0

              if (disk_testR .gt. 0.d0) then
                xc = disk_testR
                yc = 0.d0
                r_cyl_l = disk_testR
              elseif (disk_testZ .gt. 0.d0) then
                r_cyl_l=sqrt(xc**2+yc**2)
                zc = disk_testz
              else
                r_cyl_l=max(sqrt(xc**2+yc**2),dR)
              endif

              do i_Rad=-1,3,2 ! want to do i_rad=0 last so values are saved at end of loop.
                r_cyl = max(r_cyl_l + i_Rad*dR,dR)
                if (i_Rad .eq. 3) r_cyl = r_cyl_l

                r_sph=sqrt(r_cyl**2+zc**2)

                rcyl_ratio = r_cyl/Rin

                if (rcyl_ratio .ge. 1.d0) then
                  ! Recalculate full density - first midplane
                  rho = max(D0 * (R_cyl/R0)**(alpha),D1)/scale_d
                  pres_mid = rho*P_fact_cen0/R_cyl/Bz
                  temp = (pres_mid*Bz*disk_abar/8.314d7/rho)*(scale_l/scale_t)**2
                  if (temp .gt. Tmaxi) then
                    pres_mid = pres_mid * Tmaxi/temp
                  endif

                  max_z = maxz_fact * R_cyl
                  if (disk_trans_fact .lt. 1.0d0) then
                    ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
                    ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
                    dztrans = ztrans2 - ztrans1
                  else
                    ztrans1 = max_z
                    ztrans2 = max_z
                    dztrans = 0.d0
                  endif
                  if (abs(zc) .le. ztrans1) then
                    if (r_sph .gt. 0.d0) then
                      exp_arg = 400.d0*(1.d0 - R_cyl/r_sph)*(R_cyl/R0)**(-beta-1.d0)
                      rho = rho*exp(-exp_arg)
                    endif
                  else if (abs(zc) .ge. ztrans2) then
                    rho = rho*exp(-disk_exp_limit)
                  else
                    rtrans1 = sqrt(R_cyl**2 + ztrans1**2)
                    ftrans1 = -400.d0*(1.d0 - R_cyl/rtrans1)
                    delta_trans = ftrans1 + disk_exp_limit

                    theta_trans = -400.d0 * R_cyl * ztrans1 / rtrans1 ** 3 * dztrans
                      eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

                    k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
                    k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
                    k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

                    dz_here = ztrans2 - abs(zc)

                    exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
                    rho = rho*exp(-exp_arg)
                  endif

                  rho = max(rho,D1/scale_D)
                  pres = rho/R_cyl * P_fact_cen0

                else
                  rsph_ratio = r_sph/Rin
                  P_fact_cen = P_fact_cen0 * Mass_of_r(rcyl_ratio) * rcyl_ratio

        !          z1 = zc/Z_of_r(rcyl_ratio)
        !          rsph1 = sqrt(Rin**2 + z1**2)
                  ! Note pinching flaring anymore ...
                  z1 = max(abs(zc),dR/4.d0) ! Softening a bit 
                  rsph1 = max(r_sph,sqrt(r_cyl**2 + z1**2))

                  ! Recalculate full density - first midplane
                  rho = max(D0 * (Rin/R0)**(alpha)/scale_d, D1/scale_D) ! rho(R=Rin,0)
                  pres1 = rho/Rin * P_fact_cen0 ! P(Rin,0)

                  rho10 = rho ! rho(Rin,0) saved for later
                  rho0  = max(rho*Rho_of_r(rcyl_ratio), D1/scale_D) ! rho(R<Rin,0)
                  pres0 = pres1*Pres_of_r(rcyl_ratio) ! pres(R<Rin,0)
                  delta_pres = (disk_raiseP_factor*rho - rho)*P_fact_cen/R_cyl ! Note P_fact_cen has a R_cyl/Rin factor
                  pres_mid = (disk_raiseP_factor*pres0 - delta_pres)/Bz
                  temp = (pres_mid*Bz*disk_abar/8.314d7/rho)*(scale_l/scale_t)**2
                  if (temp .gt. Tmaxi) then
                    pres_mid = pres_mid * Tmaxi/temp
                  endif
        !          max_z = maxz_fact * Rin

                  alpha2 = (1.d0/disk_aspect**2)*(1.d0 - (1.d0 - disk_delta_alpha)*(1 - Rcyl_ratio)**3)**2
                  maxz_fact2 = sqrt( 1.d0/(1.d0 - disk_exp_limit/alpha2)**2 - 1.d0)

                  max_z = maxz_fact2 * R_cyl

                  if (disk_trans_fact .lt. 1.0d0) then
                    ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
                    ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
                    dztrans = ztrans2 - ztrans1
                  else
                    ztrans1 = max_z
                    ztrans2 = max_z
                    dztrans = 0.d0
                  endif

                  if (abs(z1) .le. ztrans1) then
                    if (rsph1 .gt. 0.d0) then
                      exp_arg = alpha2*(1.d0 - R_cyl/rsph1)*(R_cyl/R0)**(-beta-1.d0)
                      rho = rho*exp(-exp_arg)
                    endif
                  else if (abs(z1) .ge. ztrans2) then
                    rho = rho*exp(-disk_exp_limit)
                  else
                    rtrans1 = sqrt(R_cyl**2 + ztrans1**2)
                    ftrans1 = -alpha2*(1.d0 - R_cyl/rtrans1)
                    delta_trans = ftrans1 + disk_exp_limit

                    theta_trans = -alpha2 * R_cyl * ztrans1 / rtrans1 ** 3 * dztrans
                      eta_trans = theta_trans / ztrans1 * (1.d0 - 3.d0 * ztrans1**2 / rtrans1**2) * dztrans

                    k5 = ( 6.d0*delta_trans + 3.d0*theta_trans + 0.5d0*eta_trans)/dztrans**5
                    k4 =-(15.d0*delta_trans + 7.d0*theta_trans +       eta_trans)/dztrans**4
                    k3 = (10.d0*delta_trans + 4.d0*theta_trans + 0.5d0*eta_trans)/dztrans**3

                    dz_here = ztrans2 - abs(z1)

                    exp_arg = -(k5*dz_here**5 + k4*dz_here**4 + k3*dz_here**3) + disk_exp_limit
                    rho = rho*exp(-exp_arg)
                  endif

                  rho1 = max(rho,D1/scale_D) ! rho(Rin,z_tilde)
                  rho  = max(rho*Rho_of_r(rcyl_ratio),D1/scale_d)

        !  Can't do the originally smooth P version, since it sets P above a background, so high T above disk.
        !   However, can raise it a bit to try and soften central shear. Basically, pres0 and rho0 terms are scaled by disk_raiseP_factor
                  delta_pres = (disk_raiseP_factor*rho0 - rho)*P_fact_cen/disk_aspect**2/alpha2/max(R_cyl,dR) ! Note P_fact_cen has a R_cyl/Rin factor
                  pres = disk_raiseP_factor*pres0 - delta_pres
                  rho = rho*Mass_of_r(rcyl_ratio)/Mass_of_r(rsph_ratio)

        !          pres = rho*P_fact_cen/R_cyl

                endif

                temp = (pres*disk_abar/8.314d7/rho)*(scale_l/scale_t)**2
                if (temp .gt. Tmaxi) then
                  pres = pres * Tmaxi/temp
                endif
                if (i_Rad .eq. -1) then
                  pres_l = pres
                endif
                if (i_Rad .eq. 1) then
                  pres_r = pres
                endif

      ! Done loop
              enddo

              ! Now to get Pressure gradient
              dPdR = (pres_R - pres_L)/(r_cyl + dR - max(r_cyl - dR,dR) )

              if (disk_testR .lt. 0.d0) then
                if (r_sph .gt. 0.d0) then
                  v_cyl = sqrt(max(GM_code*Mass_of_R(r_sph/Rin)*r_cyl**2/r_sph**2/max(Rin,r_sph) + disk_shear_forcing*r_cyl*dPdR/rho,0.d0))
                else
                  v_cyl = 0.d0
                endif
              else
                v_cyl = 0.d0
              endif



              ! Add to totals
              dens_tot = dens_tot + rho
              p_tot = p_tot + pres
              pm_tot = pm_tot + pres_mid
              if (r_cyl .gt. 0.d0) then
                mom_tot(1) = mom_tot(1) - v_cyl*yc/r_cyl*rho
                mom_tot(2) = mom_tot(2) + v_cyl*xc/r_cyl*rho
              endif

              ! Finish subsample loops
            enddo ! ii
          enddo ! jj
        enddo ! kk

        ! average totals
        mom_tot = mom_tot / dens_tot
        rho = dens_tot/8**disk_nsubsample
        pres = p_tot/8**disk_nsubsample
        pres_mid = pm_tot/8**disk_nsubsample

        q(i,1)=  rho
        call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
        ! note, first get theta on 0 : 2pi
        !       then get z on -1 to 1
        !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
        !  To me this looks like it would favour the poles.
        !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
        q(i,2)= mom_tot(1) + v_turb(1)
        q(i,3)= mom_tot(2) + v_turb(2)
        q(i,4)= v_turb(3)
        sum_dust=0.0d0
#if NDUST>0
        if(mrn) call init_dust_ratio(epsilon_0, dustMRN)
        do idust =1,ndust
           q(i, firstindex_ndust+idust)= dust_ratio(idust)/(1.0d0+dust_ratio(idust))
           if(mrn) q(i, firstindex_ndust+idust) = dustMRN(idust)
           sum_dust = sum_dust + q(i, firstindex_ndust+idust)
        end do   
#endif
        q(i,5)=  pres*(1.0d0-sum_dust)

        ! Left B fields - Assuming Bz is beta for z, pres is total pressure. P_B = B^2/2
        q(i,6)     = 0.d0
        q(i,7)     = 0.d0
        q(i,8)     = sqrt(2.d0 * pres_Mid )

        ! Right B fields. Div * B is zero with linear operator (Bxr - Bxl)/dx ...
        q(i,nvar+1)= 0.d0
        q(i,nvar+2)= 0.d0
        q(i,nvar+3)= q(i,8)

    enddo


  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
  ! kinetic energy
  u(1:nn,5)=0.0d0
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,2)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,3)**2
  u(1:nn,5)=u(1:nn,5)+0.5*q(1:nn,1)*q(1:nn,4)**2
  ! pressure -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+q(1:nn,5)/(gamma-1.0d0)
  ! magnetic energy -> total fluid energy
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,6)+q(1:nn,nvar+1))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,7)+q(1:nn,nvar+2))**2
  u(1:nn,5)=u(1:nn,5)+0.125d0*(q(1:nn,8)+q(1:nn,nvar+3))**2
  u(1:nn,6:8)=q(1:nn,6:8)
  u(1:nn,nvar+1:nvar+3)=q(1:nn,nvar+1:nvar+3)
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do irad=1,nener
     u(1:nn,8+irad)=q(1:nn,8+irad)/(gamma_rad(irad)-1.0d0)
     u(1:nn,5)=u(1:nn,5)+u(1:nn,8+irad)
  enddo
#endif
#if NVAR>8+NENER
  ! passive scalars
  do ivar=9+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif
#if NDUST>0
     ! dust
     do idust=1,ndust
        u(1:nn,firstindex_ndust+idust)=q(1:nn,1)*q(1:nn,firstindex_ndust+idust)
     end do
#endif 
end subroutine condinit
!================================================================
!================================================================
!================================================================

function Mass_of_r(x)
  implicit none
  double precision :: x, Mass_of_r

  if (x .lt. 1.d0) then
    Mass_of_r = ((99.d0*x**3 - 154.d0*x**4  + 63.d0*x**5)**2)/64.d0
  else
    Mass_of_r = 1.d0
  endif

  return
end function Mass_of_r

function MPrime_of_r(x)
  implicit none
  double precision :: x, MPrime_of_r

  if (x .lt. 1.d0) then
    MPrime_of_r = ( 99.d0*x**3 - 154.d0*x**4  +  63.d0*x**5)* &
                  (297.d0*x**2 - 616.d0*x**3  + 315.d0*x**4)/32.d0
  else
    MPrime_of_r = -1.d0
  endif

  return
end function MPrime_of_r

function Rho_of_r(x)
  implicit none
  double precision :: x, Rho_of_r

  ! Assumes constant 1st, 2nd, 3rd derivatives, delta = 2, alpha = -1.5
  if (x .lt. 1.d0) then
    Rho_of_r = 2.d0 - (31.d0*x**4 + 69.d0*x**5 - 143.d0*x**6 + 59.d0*x**7)/16.d0
  else
    Rho_of_r = 1.d0
  endif

  return
end function Rho_of_r

function RhoPrime_of_r(x)
  implicit none
  double precision :: x, RhoPrime_of_r

  ! Assumes constant 1st, 2nd, 3rd derivatives, delta = 2, alpha = -1.5
  if (x .lt. 1.d0) then
    RhoPrime_of_r = -(124.d0*x**3 + 345.d0*x**4 - 858.d0*x**5 + 413.d0*x**6)/16.d0
  else
    RhoPrime_of_r = 1.5d0
  endif

  return
end function RhoPrime_of_r

function Pres_of_r(x)
  implicit none
  double precision :: x, Pres_of_r

  ! Assumes constant 1st, 2nd, 3rd derivatives, delta = 3, alpha = -2.5
  if (x .lt. 1.d0) then
    Pres_of_r = 3.d0 - (275.d0*x**4 + 261.d0*x**5 - 795.d0*x**6 + 355.d0*x**7)/48.d0
  else
    Pres_of_r = 1.d0
  endif

  return
end function Pres_of_r

function PPrime_of_r(x)
  implicit none
  double precision :: x, PPrime_of_r

  ! Assumes constant 1st, 2nd, 3rd derivatives, delta = 3, alpha = -2.5
  if (x .lt. 1.d0) then
    PPrime_of_r = 3.d0 - (275.d0*x**4 + 261.d0*x**5 - 795.d0*x**6 + 355.d0*x**7)/48.d0
  else
    PPrime_of_r = 2.5d0
  endif

  return
end function PPrime_of_r

function Z_of_r(x)
  implicit none
  double precision :: x, Z_of_r

  ! Assumes a Z squish factor of 0.5, so z0/z1 = 0.5.
  if (x .lt. 1.d0) then
    Z_of_r = 0.5 + x**3 - 0.5d0*x**4
  else
    Z_of_r = 1.d0
  endif

  return
end function Z_of_r

function ZPrime_of_r(x)
  implicit none
  double precision :: x, ZPrime_of_r

  ! Assumes a Z squish factor of 0.5, so z0/z1 = 0.5.
  if (x .lt. 1.d0) then
    ZPrime_of_r = 3.d0*x**2 - 2.d0*x**3
  else
    ZPrime_of_r = 1.d0
  endif

  return
end function ZPrime_of_r

    !   v(r) = v(R) * [  (99/8)(r/R)^3  -  (77/4)(r/R)^4  + (63/8)(r/R)^5  ] == v(R) * P5(r/R)
    !     --> same as if M(r<R) = Mtot * (r/R) * P5(r/R)^2
    !   rho & pres = rho(R) & pres(R) * [ Delta - P7(r/R;Delta,alpha) ]
    !          --> rho : Delta = 2 ; alpha = -1.5
    !               --> P7(r/R;2,-1.5) = (31/16)(r/R)^4 + (69/16)(r/R)^5 - (143/16)(r/R)^6 + (59/16)(r/R)^7
    !          --> Pres : Delta = 3 ; alpha = -2.5
    !               --> P7(r/R;3,-2.5) = (275/48)(r/R)^4 + (87/16)(r/R)^5 - (265/16)(r/R)^6 + (355/48)(r/R)^7

!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy,zz,vx,vy,vz,rr,tt,omega,aa,twopi

  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)
     yy=x(i,2)
     zz = 0.d0
#if NDIM>2
     zz=x(i,3)
#endif
     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

     v(i,1)=vx
     v(i,2)=vy
     v(i,3)=vz

  end do


end subroutine velana



function f(theta)

  implicit none

  double precision :: f, theta

  return

end function f



function g(theta)

  implicit none

  double precision :: g, theta
  double precision :: del_mid, del_cor, theta_trans, del, pi

  pi=ACOS(-1.0d0)

  del_mid = 0.1d0
  del_cor = 0.3d0
  theta_trans = 0.3d0

  if (theta .lt. pi/2.d0) then
    del = (pi/2. - theta_trans) - theta
  else
    del = theta - (pi/2. + theta_trans)
  endif

  g = del/del_mid
  g = (exp(2*g) - 1) / (exp(2*g) + 1)
  g = 0.5d0*(g + 1.d0)
  g = g * ( del_cor - del_mid + (0.5 - del_cor)*max(del,0.d0)/(pi/2. - theta_trans) )
  g = (del_mid + g)**2

  return
end function g

function dPdz_R(z,params)
  implicit none

  ! Note, need dPdz to be continuous, even at Rin
  real(8) :: dPdz_R, z, alpha2, R0, D0, D1, T0, R, Rin, p,q, R_IG, GMsun,rho,grav_z,abar
  real(8) :: r_bound, exp_limit
  real(8), dimension(20) :: params

  D0     = params(1)
  GMsun  = params(2)
  exp_limit = params(3)
  R0     = params(4)
  alpha2 = params(5)
  p      = params(6)
  q      = params(7)
  Abar   = params(8)
  R      = params(9)
  D1     = params(10)
  Rin    = params(11)


  if (R .lt. Rin) then
    R_bound = sqrt(Rin**2+z**2)
    rho = max(D0*(Rin/R0)**p*exp(-min(alpha2*(1.d0 - Rin/R_bound)*(Rin/R0)**(-q-1.0d0),exp_limit)),D1)
  else
    R_bound = sqrt(R**2+z**2)
    rho = max(D0*(R/R0)**p*exp(-min(alpha2*(1.d0 - R/R_bound)*(R/R0)**(-q-1.0d0),exp_limit)),D1)
  endif
  grav_z = -(GMsun*z)/R_bound**3
  dPdz_R = rho*grav_z

  return

end function dPdz_R

function dPdz2_R(z,params)
  implicit none

  ! Note, need dPdz to be continuous, even at Rin
  real(8) :: dPdz2_R, z, alpha2, R0, D0, D1, T0, R, Rin, p,q, R_IG, GMsun,rho,grav_z,abar
  real(8) :: r_bound, exp_limit,Reff,Rsph,testR,testZ

  real(8) :: grav_angle,grav_width,zmax,zmin,dz,gm_rescale,pi
  real(8), dimension(20) :: params

  pi = acos(-1.d0)

  D0     = params(1)
  GMsun  = params(2)
  exp_limit = params(3)
  R0     = params(4)
  alpha2 = params(5)
  p      = params(6)
  q      = params(7)
  Abar   = params(8)
  R      = params(9)
  D1     = params(10)
  Rin    = params(11)
  testR = params(12)
  testz = params(13)
  grav_angle = params(14)
  grav_width = params(15)

  if (testz .lt. 0.d0) then
    Reff = max(R, sqrt(Rin**2 - z**2) )
    Rsph = R**2 + z**2
    if (Rsph .eq. 0.d0) then
      Rsph = 1.d0
    else
      Rsph = sqrt(Rsph)
    endif

    rho = max(D0*(Reff/R0)**p*exp(-min(alpha2*(1.d0 - Reff/max(Rsph,Rin))*(Reff/R0)**(-q-1.0d0),exp_limit)),D1)

    if (grav_angle .gt. 0.d0) then
      zmin = grav_angle*(1.d0 - 0.5d0*grav_width)*Reff
      zmax = grav_angle*(1.d0 + 0.5d0*grav_width)*Reff
      dz = zmax - zmin
      gm_rescale = GMsun * (cos( min(1.d0,max(0.d0,(abs(z)-zmin)/dz))*pi ) + 1.d0)/2.d0
    else
      gm_rescale = GMsun
    endif
    Reff = max(Rsph, Rin)
    grav_z = -gm_rescale*z/Reff**3

    dPdz2_R = rho*grav_z
  else
    dPdz2_R = 0.d0
  endif

  return

end function dPdz2_R

function dPdR_z0(R,params)
  implicit none

  real(8) :: dPdR_z0, z, R0, D0, D1,R,Rin,p,q, R_IG,GMsun,rho,grav_R,alpha2,abar
  real(8) :: fcent,PI=3.1415926535897932384626433832795028841971694d0

  real(8), dimension(20) :: params
  D0     = params(1)
  GMsun  = params(2)
  R0     = params(4)
  alpha2 = params(5)
  p      = params(6)
  q      = params(7)
  Abar   = params(8)
  D1     = params(10)
  Rin    = params(11)

  if (R .lt. Rin) then
    dPdR_z0 = 0.d0
  else
    rho = max(D0*(R/R0)**p,D1)
    grav_R = -(GMsun)/R**2
    fcent = -grav_R*(1.d0 + (p+q)/alpha2)
    dPdR_z0 = rho*(fcent + grav_R)
  endif

  return

end function dPdR_z0

function dPdR2_z0(R,params)
  implicit none

  real(8) :: dPdR2_z0, z, R0, D0, D1,R,Rin,p,q, R_IG,GMsun,rho,grav_R,alpha2,abar,testR,testz
  real(8) :: fcent,PI=3.1415926535897932384626433832795028841971694d0,Reff,yZero,yOne,hR_atan

  real(8), dimension(20) :: params
  D0     = params(1)
  GMsun  = params(2)
  R0     = params(4)
  alpha2 = params(5)
  p      = params(6)
  q      = params(7)
  Abar   = params(8)
  D1     = params(10)
  Rin    = params(11)
  testR = params(12)
  testz = params(13)

  Reff = max(R,Rin)

  rho = max(D0*(Reff/R0)**p,D1)
  grav_R = -(GMsun)/Reff**2
  hR_atan = 10.d0
  yZero = atan(-hR_atan*Rin/2.d0) ! so non Keplarian goes to zero at Rin/2.
  yOne  = atan( hR_atan*Rin     ) ! so non Keplarian at full effect at 2Rin
  fcent = -grav_R*(1.d0 + (p+q)/alpha2*min(max((atan((R-Rin)*hR_atan) - yZero)/(yOne - yZero),0.d0), 1.d0) )
  dPdR2_z0 = rho*(fcent + grav_R)

  ! Note, I'm integrating as if dPdR != 0 so I can compare with the actual disk.

  return

end function dPdR2_z0

function dPdR2_z(R,params)
  implicit none

  real(8) :: dPdR2_z, z, R0, D0, D1,R,Rin,p,q, R_IG,GMsun,rho,grav_R,alpha2,abar,testR,testz,rsph
  real(8) :: fcent,Reff,yZero,yOne,hR_atan
  real(8) :: grav_angle,grav_width,zmax,zmin,dz,gm_rescale,pi

  real(8), dimension(20) :: params

  pi = acos(-1.d0)
  D0     = params(1)
  GMsun  = params(2)
  R0     = params(4)
  alpha2 = params(5)
  p      = params(6)
  q      = params(7)
  Abar   = params(8)
  D1     = params(10)
  Rin    = params(11)
  z      = params(9)
  testR = params(12)
  testz = params(13)
  grav_angle = params(14)
  grav_width = params(15)

  rsph = sqrt(R**2 + z**2)

  Reff = max(R,Rin)
  rho = max(D0*(Reff/R0)**p,D1)

  if (rsph .lt. Rin) then
    Reff = sqrt(Rin**2 - z**2)
  else
    Reff = R
  endif
  rsph = max(rsph,Rin)

  if (grav_angle .gt. 0.d0) then
    zmin = grav_angle*(1.d0 - 0.5d0*grav_width)*Reff
    zmax = grav_angle*(1.d0 + 0.5d0*grav_width)*Reff
    dz = zmax - zmin
    gm_rescale = GMsun * (cos( min(1.d0,max(0.d0,(abs(z)-zmin)/dz))*pi ) + 1.d0)/2.d0
  else
    gm_rescale = GMsun
  endif
  grav_R = -gm_rescale*Reff/Rsph**3

  hR_atan = 10.d0
  yZero = atan(-hR_atan*Rin/2.d0) ! so non Keplarian goes to zero at Rin/2.
  yOne  = atan( hR_atan*Rin     ) ! so non Keplarian at full effect at 2Rin
  fcent = -grav_R*(1.d0 + ( (p+q)/alpha2 + q*(1.d0 - Reff/sqrt(Reff**2 + z**2)))*min(max((atan((R-Rin)*hR_atan) - yZero)/(yOne - yZero),0.d0), 1.d0) )
  dPdR2_z = rho*(fcent + grav_R)

  ! Note, I'm integrating as if dPdR != 0 so I can compare with the actual disk.

  return

end function dPdR2_z



function dPdz_R_gal(z,params)
  implicit none

  ! Note, need dPdz to be continuous, even at Rin
  real(8) :: dPdz_R_gal, z, D0, D1, R_cyl, rho,hR,hz,r_sph
  real(8) :: fgrav_z,sech,Rout,Hout,alpha
  real(8), dimension(20) :: params

  external :: fgrav_z,sech

  D0     = params(1)
  hR     = params(2)
  hz     = params(3)
  r_cyl  = params(4)
  D1     = params(5)
  Rout   = params(9)
  Hout   = Rout*hz/hR
  alpha  = params(10)

  r_sph = sqrt(r_cyl**2+z**2)
  if ( sqrt( (r_cyl/Rout)**2 + (z/Hout)**2) .gt. 1.d0) then
    rho = D1*(r_sph/Rout)**alpha
  else
    rho = max(D0*sech(r_cyl/hR)*sech(z/hz),D1*(r_sph/Rout)**alpha)
  endif
  dPdz_R_gal = rho*fgrav_z(z,r_cyl,r_sph,params) ! remember, fgrav in scale_l/scale_t**2

  return

end function dPdz_R_gal



function sech(x)
  implicit none

  real(8) :: sech,x

  sech = 2.d0/( exp(x) + exp(-x) )

  return

end function sech



function fgrav_R(zc,r_cyl,r_sph,params)
  implicit none
  real(8) :: fgrav_R, zc, r_cyl, r_sph, fr, fz
  real(8),dimension(20) :: params
  real(8) :: GM,R200,c,fofx,dummy,scale_l,scale_t
  real(8) :: Msun=1.989d33, kpc=3.085678d21

  real(8), save :: dr_b, dr_d, dz_d

  integer, parameter :: nbulge=1000, nzdisk=3334, nrdisk=84
  real(8), save, dimension(0:nbulge) :: save_bulge
  real(8), save, dimension(0:nzdisk,0:nrdisk) :: save_disk
  integer :: ir, iz, ir1,iz1,ir2,iz2
  logical, save :: first=.true.

  real(8) :: acc_bulge, acc_disk, acc_NFW

  external :: fofx
  ! File units are in kpc , (km/s)**2/kpc, but length_units = kpc, mass_unit = 1e10 Msun ???, and then
  !  time_units = such that G = 1 in code units is 3.8433374223e10. Not t_unit ~ 1/sqrt(d_unit) ~ 1/sqrt(m_unit)
  if (first) then ! read data
    first = .false.
    scale_l = params(11)
    scale_t = params(12)

    ! Read Bulge info: first read r into dr_b, second into fr. Bulge is spherical
    open(42,file='Bulge_acc_1d.dat',status='old')
    read(42,*) dr_b, save_bulge(0)
    read(42,*) fr, save_bulge(1)
    do ir=2,nbulge
      read(42,*) dummy, save_bulge(ir)
    enddo
    ! Units are (km/s)**2 / kpc
    save_bulge = save_bulge * 1d10 / kpc
    ! Units are cm/s**2
    save_bulge = save_bulge/scale_l*scale_t**2
    ! Units are in code_length now

    ! set dr / dz
    dr_b = fr - dr_b
    close(42)

    ! Read Disk info: first read r into dr_b, second into fr
    open(42,file='Gas_acc.dat',status='old')
    read(42,*) dr_d, dz_d, save_disk(0,0)
    read(42,*) dummy, fz, save_disk(1,0)
    do iz=2,nzdisk
      read(42,*) dummy, dummy, save_disk(iz,0)
    enddo
    read(42,*) fr, dummy, save_disk(0,1)
    do iz=1,nzdisk
      read(42,*) dummy, dummy, save_disk(iz,1)
    enddo
    do ir=2,nrdisk
      do iz=0,nzdisk
        read(42,*) dummy, dummy, save_disk(iz,ir)
      enddo
    enddo
    ! Units are (km/s)**2 / kpc
    save_disk = save_disk * 1d10 / kpc
    ! Units are cm/s**2
    save_disk = save_disk/scale_l*scale_t**2
    ! Units are in code_length now

    ! set dr / dz
    dr_d = fr - dr_d
    dz_d = fz - dz_d
    close(42)

    ! Halo NFW: analytic
  endif

  ! Assumes inputs are in kpc!
  ! bulge:
  fr = r_sph/dr_b
  ir1 = int(fr)
  ir2 = min(ir1+1,nbulge)
  acc_bulge = (fr*save_bulge(ir2) + (1.d0-fr)*save_bulge(ir1))*r_cyl/r_sph

  ! disk
  fr = r_cyl/dr_d
  ir1 = int(fr)
  ir2 = min(ir1+1,nrdisk)
  fz = abs(zc)/dz_d
  iz1 = int(fz)
  iz2 = min(iz1+1,nzdisk)
  acc_disk = ( fz*(fr*save_disk(iz2,ir2) + (1.d0-fr)*save_disk(iz2,ir1)) + &
        (1.d0-fz)*(fr*save_disk(iz1,ir2) + (1.d0-fr)*save_disk(iz1,ir1)) )*r_cyl/r_sph

  ! acc_NFW
  GM     = params(6)
  R200   = params(7)
  c      = params(8)
  ! GM has units scale_l**3 / scale_t**2
  !  fofx is unitless
  ! (r_cyl / r_sph**3) has units 1/scale_l**2
  ! We have units scale_l / scale_t**2 --> correct units for acceleration in code length

  ! My opinion here is the onus is on the input of the files to match these units.
  acc_NFW = -GM*fofx(r_sph/(R200/c))/fofx(c)*r_cyl/max(r_sph,1d-10)**3

  ! return fgrav_R
  fgrav_R = acc_bulge + acc_disk + acc_NFW

  return
end function fgrav_R



function fgrav_z(zc,r_cyl,r_sph,params)
  implicit none
  real(8) :: fgrav_z, zc, r_cyl, r_sph, fr, fz,scale_l,scale_t
  real(8),dimension(20) :: params
  real(8) :: GM,R200,c,fofx,dummy
  real(8) :: Msun=1.989d33, kpc=3.085678d21

  real(8), save :: dr_b, dr_d, dz_d

  integer, parameter :: nbulge=1000, nzdisk=3334, nrdisk=84
  real(8), save, dimension(0:nbulge) :: save_bulge
  real(8), save, dimension(0:nzdisk,0:nrdisk) :: save_disk
  integer :: ir, iz, ir1,iz1,ir2,iz2
  logical, save :: first=.true.

  real(8) :: acc_bulge, acc_disk, acc_NFW

  external :: fofx
  ! File units are in kpc , (km/s)**2/kpc
  if (first) then ! read data
    first = .false.
    scale_l = params(11)
    scale_t = params(12)
    ! Read Bulge info: first read r into dr_b, second into fr. Bulge is spherical
    open(42,file='Bulge_acc_1d.dat',status='old')
    read(42,*) dr_b, save_bulge(0)
    read(42,*) fr, save_bulge(1)
    do ir=2,nbulge
      read(42,*) dummy, save_bulge(ir)
    enddo
    ! Units are (km/s)**2 / kpc
    save_bulge = save_bulge * 1d10 / kpc
    ! Units are cm/s**2
    save_bulge = save_bulge/scale_l*scale_t**2
    ! Units are in code_length now

    ! set dr / dz
    dr_b = fr - dr_b
    close(42)

    ! Read Disk info: first read r into dr_b, second into fr
    open(42,file='Gas_acc.dat',status='old')
    read(42,*) dr_d, dz_d, dummy, save_disk(0,0)
    read(42,*) dummy, fz, dummy, save_disk(1,0)
    do iz=2,nzdisk
      read(42,*) dummy, dummy, dummy, save_disk(iz,0)
    enddo
    read(42,*) fr, dummy, dummy, save_disk(0,1)
    do iz=1,nzdisk
      read(42,*) dummy, dummy, dummy, save_disk(iz,1)
    enddo
    do ir=2,nrdisk
      do iz=0,nzdisk
        read(42,*) dummy, dummy, dummy, save_disk(iz,ir)
      enddo
    enddo
    ! Units are (km/s)**2 / kpc
    save_disk = save_disk * 1d10 / kpc
    ! Units are cm/s**2
    save_disk = save_disk/scale_l*scale_t**2
    ! Units are in code_length now

    ! set dr / dz
    dr_d = fr - dr_d
    dz_d = fz - dz_d
    close(42)

    ! Halo NFW: analytic
  endif

  ! Assumes inputs are in kpc!
  ! bulge:
  fr = r_sph/dr_b
  ir1 = int(fr)
  ir2 = min(ir1+1,nbulge)
  acc_bulge = (fr*save_bulge(ir2) + (1.d0-fr)*save_bulge(ir1))*zc/max(r_sph,min(dr_b,dr_d)/10.d0)

  ! disk
  fr = r_cyl/dr_d
  ir1 = int(fr)
  ir2 = min(ir1+1,nrdisk)
  fz = abs(zc)/dz_d
  iz1 = int(fz)
  iz2 = min(iz1+1,nzdisk)
  acc_disk = ( fz*(fr*save_disk(iz2,ir2) + (1.d0-fr)*save_disk(iz2,ir1)) + &
        (1.d0-fz)*(fr*save_disk(iz1,ir2) + (1.d0-fr)*save_disk(iz1,ir1)) )*zc/max(r_sph,min(dr_b,dr_d)/10.d0)

  ! acc_NFW
  GM     = params(6)
  R200   = params(7)
  c      = params(8)
  acc_NFW = -GM*Msun*fofx(r_sph/(R200/c))/fofx(c)*zc/max(r_sph,min(dr_b,dr_d)/10.d0)**3

  ! return fgrav_R
  fgrav_z = acc_bulge + acc_disk + acc_NFW

  return
end function fgrav_z



function fofx(x)
  implicit none
  real(8) :: fofx, x

  fofx = log(1+x) - x/(1 + x)
  return
end function fofx




subroutine Romberg(a,b,eps,x,ri,func,func_params)

  implicit none

  integer :: k, j
  real(8) ::  a, b, eps, x, ri, ri0, func, h, fa, fb, t1, t2, s, s1, s2, c1, c2, ff, err
  real(8), dimension(20) :: func_params

  integer, parameter :: maxN=20
  real(8), dimension(maxN,maxN) :: rom

  external func

  h = b-a
  fa = func(a,func_params)
  fb = func(b,func_params)
  t1 = h*(fa+fb)/2.0d0
  k = 1
  j=1
  rom = 0.d0
  rom(1,1) = t1

  do k=2,maxN

     s = 0.0d0
     x = a+0.5d0*h

     do while((x-b)*(x-h-b) .GT. 0.0d0)
       ff = func(x,func_params)
       s = s+ff
       x = x+h
     enddo

     rom(k,1) = (rom(k-1,1)+h*s)/2.0d0
     h = h/2.d0
     do j=2,k
       rom(k,j) = rom(k,j-1) + (rom(k,j-1)-rom(k-1,j-1))/(4.d0**(j-1)-1.d0)
     enddo

     err = ABS(rom(k,k)-rom(k-1,k-1))/max(ABS(rom(k,k)),eps/1.d3)
     if ( err .LT. eps .or. k .EQ. maxN) then
       ri = rom(k,k)
       exit
     endif
  enddo

  return

end subroutine Romberg

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

subroutine mag_compute(ilevel)
  use amr_commons
  !use pm_commons
  use hydro_commons
  !use random
  implicit none
  integer::i,j,ilevel,icell
  logical::nogal
  real(dp)::Axdl,Axdr,Axul,Axur
  real(dp)::Aydl,Aydr,Ayul,Ayur
  real(dp)::Azdl,Azdr,Azul,Azur
  real(dp)::Bxl,Bxr,Byl,Byr,Bzl,Bzr
  real(dp),dimension(1:3)::pos,cell_center

  real(dp)::dx,dxhalf,dxmin,dxminhalf,scale,dx_loc,vol_loc
  integer::nx_loc,ind,ix,iy,iz,iskip,nfine
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  dxhalf = 0.5D0*dx
  dxmin = 0.5D0**nlevelmax
  dxminhalf = 0.5D0*dxmin
  nfine = 2**(nlevelmax-ilevel)

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! compute field
  do ind=1,twotondim
    iskip=ncoarse+(ind-1)*ngridmax
    do j=1,active(ilevel)%ngrid
      icell=active(ilevel)%igrid(j)
      cell_center = xg(icell,:)+xc(ind,:)-skip_loc(:)

      ! edge coordinates
      ! Ax
      pos = cell_center
      pos(1) = pos(1) - dxhalf + dxminhalf
      pos(2) = pos(2) - dxhalf
      pos(3) = pos(3) - dxhalf
      Axdl=0.0;Axdr=0.0;Axul=0.0;Axur=0.0
      do i=1,nfine
        CALL mag_potential(pos,1,Axdl)
        pos(1) = pos(1) + dxmin
      enddo
      pos(2) = pos(2) + dx
      do i=1,nfine
        pos(1) = pos(1) - dxmin
        CALL mag_potential(pos,1,Axdr)
      enddo
      pos(3) = pos(3) + dx
      do i=1,nfine
        CALL mag_potential(pos,1,Axur)
        pos(1) = pos(1) + dxmin
      enddo
      pos(2) = pos(2) - dx
      do i=1,nfine
        pos(1) = pos(1) - dxmin
        CALL mag_potential(pos,1,Axul)
      enddo
      ! Ay
      pos = cell_center
      pos(1) = pos(1) - dxhalf
      pos(2) = pos(2) - dxhalf + dxminhalf
      pos(3) = pos(3) - dxhalf
      Aydl=0.0;Aydr=0.0;Ayul=0.0;Ayur=0.0
      do i=1,nfine
        CALL mag_potential(pos,2,Aydl)
        pos(2) = pos(2) + dxmin
      enddo
      pos(1) = pos(1) + dx
      do i=1,nfine
        pos(2) = pos(2) - dxmin
        CALL mag_potential(pos,2,Aydr)
      enddo
      pos(3) = pos(3) + dx
      do i=1,nfine
        CALL mag_potential(pos,2,Ayur)
        pos(2) = pos(2) + dxmin
      enddo
      pos(1) = pos(1) - dx
      do i=1,nfine
        pos(2) = pos(2) - dxmin
        CALL mag_potential(pos,2,Ayul)
      enddo
      ! Az
      pos = cell_center
      pos(1) = pos(1) - dxhalf
      pos(2) = pos(2) - dxhalf
      pos(3) = pos(3) - dxhalf + dxminhalf
      Azdl=0.0;Azdr=0.0;Azul=0.0;Azur=0.0
      do i=1,nfine
        CALL mag_potential(pos,3,Azdl)
        pos(3) = pos(3) + dxmin
      enddo
      pos(1) = pos(1) + dx
      do i=1,nfine
        pos(3) = pos(3) - dxmin
        CALL mag_potential(pos,3,Azdr)
      enddo
      pos(2) = pos(2) + dx
      do i=1,nfine
        CALL mag_potential(pos,3,Azur)
        pos(3) = pos(3) + dxmin
      enddo
      pos(1) = pos(1) - dx
      do i=1,nfine
        pos(3) = pos(3) - dxmin
        CALL mag_potential(pos,3,Azul)
      enddo

      ! Bx left
      Bxl = ((Azul - Azdl) - (Ayul - Aydl))/dx / nfine
      uold(icell+iskip,6)      = uold(icell+iskip,6) + Bxl
      ! By left
      Byl = ((Axul - Axdl) - (Azdr - Azdl))/dx / nfine
      uold(icell+iskip,7)      = uold(icell+iskip,7) + Byl
      ! Bz left
      Bzl = ((Aydr - Aydl) - (Axdr - Axdl))/dx / nfine
      uold(icell+iskip,8)      = uold(icell+iskip,8) + Bzl

      ! Bx right
      Bxr = ((Azur - Azdr) - (Ayur - Aydr))/dx / nfine
      uold(icell+iskip,nvar+1) = uold(icell+iskip,nvar+1) + Bxr
      ! By right
      Byr = ((Axur - Axdr) - (Azur - Azul))/dx / nfine
      uold(icell+iskip,nvar+2) = uold(icell+iskip,nvar+2) + Byr
      ! Bz right
      Bzr = ((Ayur - Ayul) - (Axur - Axul))/dx / nfine
      uold(icell+iskip,nvar+3) = uold(icell+iskip,nvar+3) + Bzr
    end do
  end do

end subroutine mag_compute

subroutine mag_potential(pos,dir,Amag)
  use amr_parameters, ONLY: boxlen,dp
  use hydro_parameters
  implicit none
  real(dp)::r,h
  real(dp)::Ah,Amag,Bpow,AR
  integer::i,dir
  real(dp),dimension(1:3)::pos,gcenter,gaxis,xrel,grad,cross
  real(dp)::sBR,sBh,sR,sH,soft

  ! Note : toroidal field comes from a Bz term alone, particular one with
  !   Az ~ f(z) * exp(-R/h) --> Bx/y = -By/x = Az/h/R

  ! Note : vertical field comes from a B_R term alone, where Ax/y = - Ay/x
  !   A ~ A0 / R^m (-y/R,x/R,0) --> Bz = A0 / R^(m+1) *(1-m)

  gcenter(1) = 0.5d0 ! disk_cenX
  gcenter(2) = 0.5d0 ! disk_cenY
  gcenter(3) = 0.5d0 ! disk_cenZ

  gaxis(1) = 0.d0 ! disk_axisX
  gaxis(2) = 0.d0 ! disk_axisY
  gaxis(3) = 1.d0 ! disk_axisZ

  sR = disk_mag_scaleR / boxlen
  sH = disk_mag_scaleH / boxlen
  sBh = disk_mag_scaleBh
  sBr = disk_mag_scaleBR

  Bpow = disk_mag_pow

  soft = disk_soft

  ! coordinates in galaxy frame
  xrel = pos - gcenter
  h = DOT_PRODUCT(xrel,gaxis)
  grad = xrel - h*gaxis
  r = max(NORM2(grad),soft)

  call CROSS_PRODUCT(gaxis,xrel,cross)
  ! toroidal field
  Ah = sBh * sR * exp(-r/sR) * exp(-ABS(h)/sH)

  ! vertical field
  AR = sBr /r**(Bpow+1)

  ! vector in cartesian frame
  Amag = Amag + Ah*gaxis(dir) + AR*cross(dir)

end subroutine mag_potential

! FUNCTION DOT_PRODUCT(a, b)
!   real(8) :: dot_product
!   real(8), DIMENSION(3), INTENT(IN) :: a, b

!   dot_product = a(1) * b(1) + a(2) * b(2) + a(3) * b(3)
!   return
! END FUNCTION DOT_PRODUCT


subroutine CROSS_PRODUCT(a, b, c)
  real(8), DIMENSION(3) :: c
  real(8), DIMENSION(3) :: a, b

  c(1) = a(2) * b(3) - a(3) * b(2)
  c(2) = a(3) * b(1) - a(1) * b(3)
  c(3) = a(1) * b(2) - a(2) * b(1)
  return
END subroutine CROSS_PRODUCT
