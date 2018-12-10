!================================================================
!================================================================
!================================================================
! This routine reinforces initial conditions for RAMSES for some(x,y,z).
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
!
!  IMPORTANT: This DOES NOT reset the magnetic field or passive variables.
!   The former is because I don't know how to do that while preserving Div B = 0
!
!   The latter is because at this point we are not using any.
!
!
!================================================================
!================================================================
!================================================================
subroutine relax_condinit(xc,yc,zc,dx,rho,vel,pres)

  use amr_parameters
  use amr_commons
  use hydro_parameters
  use poisson_parameters
  implicit none

  real(dp)::xc,yc,zc,dx
  real(dp)::rho,vel(ndim),pres

  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  real(dp)::scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2

  integer::ivar,i,ii,jj,kk,i_rad
  real(dp)::x0,y0,z0,xc0,yc0,zc0
  real(dp)::xl,yl,zl,xr,yr,zr,R0,Rin,Rout,Rout2,Hout
  real(dp)::V0,P0,P1,pi,f,rr,abar,alpha,beta
  real(dp)::m,power,pm1,pp1,thetamax,GM,GM_code
  real(dp)::theta,g,D0,D1,dens,T0,r_cyl,v_cyl,r_sph,exp_arg,Rin0,Rout0,hR_clip
  real(dp)::pmin,R_integral,R_left,R_right,tmini,temp_zero,temp_minus1
  real(dp)::dR,dz,r_loop,z_loop,rho_l,rho_r,temp,temp_l,temp_r,&
             f_cent_out,f_grav_inR,f_grav_inZ,R200,R1,H1,c,Reff,p_tot,dens_tot
  real(dp)::Rmin_grav,Rmax_grav,dx_max,pres_l,pres_r,dR_scale,dZ_scale,dPdR,dPdz
  real(dp)::r_cyl_l,r_cyl_r,r_loop_l,r_loop_r,f_grav_inZ_l,f_grav_inZ_r,dPdz_R
  real(dp)::dPdz_R_gal,sech,sech_term, dPdR2_z0,dPdR2_z,dPdz2_R,pres_int,pres_intm1
  real(dp)::pres_p1,pres_l_p1,pres_r_p1,P_fact_cen0,P_fact_cen,soft
  real(dp)::v_turb(3),mom_tot(2),ztrans1,ztrans2,rtrans1,dztrans,dz_here,ftrans1,max_z, &
            maxz_fact,maxz_fact2,maxz_fact0,alpha2,aspect_ratio
  real(dp)::exp_arg_l,exp_arg_r,r_bound,eps,dPdR_z0,rom_int,f_table,fgrav_R
  real(dp)::k5,k4,k3,delta_trans,theta_trans,eta_trans,dx_ss
  real(dp)::rcyl_ratio,rsph_ratio,z1,rsph1,pres1,pres0,rho1,rho0,rho10,vcyl2,delta_pres
  real(dp),dimension(20)::func_params
  integer::nR_loop,iR,nZ_loop,iZ,iiz,iter,iter_tot,max_iter,i_table
  real(dp)::C0,C1,C2,C3,Q1,Q2,Q0,pres_zero,pres_minus1,xscl,r_sph_l,r_sph_r


  real(dp)::Mass_of_R,MPrime_of_R,Rho_of_R,RhoPrime_of_R,Pres_of_R,PPrime_of_R,Z_of_R,ZPrime_of_R
  external::Mass_of_R,MPrime_of_R,Rho_of_R,RhoPrime_of_R,Pres_of_R,PPrime_of_R,Z_of_R,ZPrime_of_R
  external :: dPdR_z0, dPdz_R, fgrav_R, dPdz_R_gal, sech, dPdR2_z0,dPdR2_z,dPdz2_R

  pi=ACOS(-1.0d0)
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

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

  elseif ( disk_setup == 26001 ) then ! Gressel
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

    Rmin_grav = gravity_params(6)
    Rmax_grav = gravity_params(7)
    gravity_params(8) = 0.d0 ! Buffer for R for dPdz calc

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

    dens_tot = 0.d0
    mom_tot = 0.d0
    p_tot = 0.d0

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

            max_z = maxz_fact * R_cyl
            ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
            ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
            dztrans = ztrans2 - ztrans1
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
            z1 = max(abs(zc),dR)
            rsph1 = max(r_sph,sqrt(r_cyl**2 + z1**2))

            ! Recalculate full density - first midplane
            rho = max(D0 * (Rin/R0)**(alpha)/scale_d, D1/scale_D) ! rho(R=Rin,0)
            pres1 = rho/Rin * P_fact_cen0 ! P(Rin,0)

            rho10 = rho ! rho(Rin,0) saved for later
            rho0  = max(rho*Rho_of_r(rcyl_ratio), D1/scale_D) ! rho(R<Rin,0)
            pres0 = pres1*Pres_of_r(rcyl_ratio) ! pres(R<Rin,0)
  !          max_z = maxz_fact * Rin

            alpha2 = 400.d0*(1.d0 - 0.75d0*(1 - Rcyl_ratio)**3)**2
            maxz_fact2 = sqrt( 1.d0/(1.d0 - disk_exp_limit/alpha2)**2 - 1.d0)

            max_z = maxz_fact2 * R_cyl

            ztrans1 = max(dx_max/2.d0,min(max_z - 2.5d0*dx_max,disk_trans_fact*max_z))
            ztrans2 = max(ztrans1+5.d0*dx_max,(2.d0 - disk_trans_fact)*max_z)
            dztrans = ztrans2 - ztrans1

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
            delta_pres = (disk_raiseP_factor*rho0 - rho)*P_fact_cen*400.d0/alpha2/max(R_cyl,dR) ! Note P_fact_cen has a R_cyl/Rin factor
            pres = disk_raiseP_factor*pres0 - delta_pres
            rho = rho*Mass_of_r(rcyl_ratio)/Mass_of_r(rsph_ratio)

  !          pres = rho*P_fact_cen/R_cyl

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

    call get_vturb(disk_vrms,sqrt(gamma*pres/rho),v_turb)
    ! note, first get theta on 0 : 2pi
    !       then get z on -1 to 1
    !       (x,y,z) is [ sqrt(1-z**2)*cos(theta), sqrt(1-z**2)*sin(theta), z ]
    !  To me this looks like it would favour the poles.
    !  Actually might not. at high z, a set width in z corresponds to a large dphi, low z is smaller dphi.
    vel(1)= mom_tot(1) + v_turb(1)
    vel(2)= mom_tot(2) + v_turb(2)
    vel(3)= v_turb(3)

end subroutine relax_condinit
!================================================================
!================================================================
!================================================================
