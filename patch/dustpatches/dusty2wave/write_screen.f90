subroutine write_screen
  use amr_commons
  use hydro_commons
  use pm_commons
  use poisson_commons
  implicit none
#if NDIM==1  
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::info
#endif
  integer::igrid,jgrid,ind,icpu
  integer::i,j,icell,ncell,ilevel,ncache
  integer::icellmin,nx_loc
  real(dp)::dx,scale,smallp,ddd,ppp

  integer     ,dimension(:),allocatable::ind_grid,ind_cell,ind_sort,ll,ll_all
  real(kind=8),dimension(:),allocatable::rr,et,ei,dd,uu,mm,gg,dtot
  real(kind=8),dimension(:),allocatable::ek,em,vv,ww,AA,BB,CC
  real(kind=8),dimension(:),allocatable::rr_all,et_all,ei_all,ek_all,em_all
  real(kind=8),dimension(:),allocatable::dd_all,uu_all,vv_all,ww_all
  real(kind=8),dimension(:),allocatable::AA_all,BB_all,CC_all
  real(kind=8),dimension(:),allocatable::mm_all,gg_all,dtot_all
#if NENER>NGRP
  real(kind=8),dimension(:,:),allocatable::prad_all,prad
#endif
#if NENER>0
  integer::irad
#endif  
#if NGRP>0
  real(kind=8),dimension(:),allocatable::er,er_all
  real(kind=8),dimension(:,:),allocatable::EE,EE_all
#if USE_M_1==1
  real(kind=8),dimension(:),allocatable::fr,fr_all
  real(kind=8),dimension(:,:),allocatable::FF,FF_all
#endif
#endif
#if NDUST>0
  real(kind=8),dimension(:,:),allocatable::eps_dust_all,eps_dust, vdust, vdust_all
#endif

  integer,dimension(1:ncpu)::iskip,ncell_loc,ncell_all
#endif  
  if(ndim>1)return

#if NDIM==1
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif
  
  ncell=0
  do ilevel=1,nlevelmax
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Count leaf cells
        do ind=1,twotondim
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))== 0)then
                 ncell=ncell+1
              end if
           end do
        end do
        deallocate(ind_grid, ind_cell)
     end if
  end do

  ncell_loc=0
  ncell_all=0
  ncell_loc(myid)=ncell
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(ncell_loc,ncell_all,ncpu,MPI_INTEGER,MPI_SUM,&
       & MPI_COMM_WORLD,info)
#endif
#ifdef WITHOUTMPI
  ncell_all=ncell_loc
#endif

  ncell=0
  iskip=0
  do icpu=1,ncpu
     iskip(icpu)=ncell
     ncell=ncell+ncell_all(icpu)
  end do

  if(myid==1)write(*,114)ncell

  if(ncell>0)then

  allocate(rr(1:ncell),mm(1:ncell),dd(1:ncell),dtot(1:ncell))
  allocate(et(1:ncell),ei(1:ncell),em(1:ncell),ek(1:ncell))
  allocate(uu(1:ncell),ll(1:ncell),gg(1:ncell))
  allocate(vv(1:ncell),ww(1:ncell),AA(1:ncell))
  allocate(BB(1:ncell),CC(1:ncell))
  allocate(rr_all(1:ncell),mm_all(1:ncell),dd_all(1:ncell),dtot_all(1:ncell))
  allocate(et_all(1:ncell),ei_all(1:ncell),em_all(1:ncell),ek_all(1:ncell))
  allocate(uu_all(1:ncell),ll_all(1:ncell),gg_all(1:ncell))
  allocate(vv_all(1:ncell),ww_all(1:ncell),AA_all(1:ncell))
  allocate(BB_all(1:ncell),CC_all(1:ncell))
  rr=0.0D0; mm=0.0D0; dd=0.0D0; dtot=0.0D0
  et=0.0D0; ei=0.0D0; em=0.0D0; ek=0.0D0
  uu=0.0D0; vv=0.0D0; ww=0.0D0; gg=0.0D0; ll=0
  AA=0.0D0; BB=0.0D0; CC=0.0D0
  rr_all=0.0D0; mm_all=0.0D0; dd_all=0.0D0; dtot_all=0.0D0
  et_all=0.0D0; ei_all=0.0D0; ek_all=0.0D0; em_all=0.0D0
  uu_all=0.0D0; vv_all=0.0D0; ww_all=0.0D0; gg_all=0.0D0; ll_all=0
  AA_all=0.0D0; BB_all=0.0D0; CC_all=0.0D0
#if NENER>NGRP
  allocate(prad(1:ncell,1:nent),prad_all(1:ncell,1:nent))
  prad=0.0D0; prad_all=0.0D0
#endif
#if NGRP>0
  allocate(er(1:ncell),er_all(1:ncell),EE(1:ncell,1:ngrp),EE_all(1:ncell,1:ngrp))
  er=0.0D0; er_all=0.0D0; EE=0.0D0; EE_all=0.0D0
#if USE_M_1==1
  allocate(fr(1:ncell,1:ngrp),fr_all(1:ncell,1:ngrp),FF(1:ncell,1:ngrp),FF_all(1:ncell,1:ngrp))
  fr=0.0D0; fr_all=0.0D0; FF=0.0D0; FF_all=0.0D0
#endif
#endif
#if NDUST>0
  allocate(eps_dust_all(1:ncell,1:ndust),eps_dust(1:ncell,1:ndust),vdust(1:ncell,1:ndust),vdust_all(1:ncell,1:ndust))
  eps_dust_all=0.0d0;eps_dust=0.0d0;vdust=0.0d0;vdust_all=0.0d0
#endif

  icell=iskip(myid)
  do ilevel=1,nlevelmax
     icellmin=icell
     ncache=numbl(myid,ilevel)
     if(ncache > 0)then
        dx=0.5D0**ilevel
        allocate(ind_grid(1:ncache),ind_cell(1:ncache))
        ! Gather all grids
        igrid=headl(myid,ilevel)
        do jgrid=1,ncache
           ind_grid(jgrid)=igrid
           igrid=next(igrid)
        end do
        ! Gather variables
        icell=icellmin
        do ind=1,twotondim
           do i=1,ncache
              ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
           end do
           do i=1,ncache
              if(son(ind_cell(i))==0)then
                 icell=icell+1
                 rr(icell)=xg(ind_grid(i),1)+(dble(ind)-1.5D0)*dx
                 ll(icell)=ilevel
              end if
           end do
        end do
        if(hydro)then
           icell=icellmin
           do ind=1,twotondim
              do i=1,ncache
                 ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
              end do
              do i=1,ncache
                 if(son(ind_cell(i))==0)then
                    icell=icell+1
                    dd(icell)=uold(ind_cell(i),1)
                    mm(icell)=dd(icell)
                    uu(icell)=uold(ind_cell(i),2)/dd(icell)
                    vv(icell)=uold(ind_cell(i),3)/dd(icell)
                    ww(icell)=uold(ind_cell(i),4)/dd(icell)
                    et(icell)=uold(ind_cell(i),5)
                    AA(icell)=0.5d0*(uold(ind_cell(i),6)+uold(ind_cell(i),nvar+1))
                    BB(icell)=0.5d0*(uold(ind_cell(i),7)+uold(ind_cell(i),nvar+2))
                    CC(icell)=0.5d0*(uold(ind_cell(i),8)+uold(ind_cell(i),nvar+3))
                    em(icell)=0.5d0*(AA(icell)**2+BB(icell)**2+CC(icell)**2)
                    ek(icell)=0.5d0*(uu(icell)**2+vv(icell)**2+ww(icell)**2)
                    ei(icell)=et(icell)-em(icell)-dd(icell)*ek(icell)
#if NENER>NGRP
                    do irad=1,nent
                       ei(icell)=ei(icell)-uold(ind_cell(i),8+irad)
                       prad(icell,irad)=(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),8+irad)
                    end do
#endif
#if NGRP>0
                    er(icell)=0.0d0
                    do irad=1,ngrp
                       EE(icell,irad)=uold(ind_cell(i),firstindex_er+irad)
                       ei(icell)=ei(icell)-EE(icell,irad)
                       er(icell)=er(icell)+EE(icell,irad)
                    end do
#if USE_M_1==1
                    fr(icell)=0.0d0
                    do j=1,ngrp
                       FF(icell,j)=uold(ind_cell(i),firstindex_er+ngrp+j)
                       fr(icell)=fr(icell)+uold(ind_cell(i),firstindex_er+ngrp+j)
                    enddo
#endif
#endif
#if NDUST>0
                    do j=1,ndust
                       eps_dust(icell,j)=uold(ind_cell(i),firstindex_ndust+j)
                       vdust(icell,j)=v_dust(ind_cell(i),j,1)
                    end do
#endif
                 end if
              end do
           end do
        end if
        if(poisson)then
           icell=icellmin
           do ind=1,twotondim
              do i=1,ncache
                 ind_cell(i)=ncoarse+(ind-1)*ngridmax+ind_grid(i)
              end do
              do i=1,ncache
                 if(son(ind_cell(i))==0)then
                    icell=icell+1
                    dtot(icell)=rho(ind_cell(i))
                    gg(icell)=f(ind_cell(i),1)
                 end if
              end do
           end do
        end if
        deallocate(ind_grid, ind_cell)
     end if
  end do

#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rr,rr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(mm,mm_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dd,dd_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(dtot,dtot_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(et,et_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ei,ei_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(em,em_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ek,ek_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(uu,uu_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vv,vv_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ww,ww_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(AA,AA_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(BB,BB_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(CC,CC_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(gg,gg_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(ll,ll_all,ncell,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(rr,rr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  rr=rr_all; mm=mm_all; dd=dd_all; dtot=dtot_all
  et=et_all; ek=ek_all; em=em_all; ei=ei_all
  uu=uu_all; vv=vv_all; ww=ww_all
  AA=AA_all; BB=BB_all; CC=CC_all
  gg=gg_all; ll=ll_all 
#if NENER>NGRP
  call MPI_ALLREDUCE(prad,prad_all,ncell*nent,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  prad=prad_all
#endif

#if NGRP>0
  call MPI_ALLREDUCE(er,er_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(EE,EE_all,ncell*ngrp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  er=er_all
  EE=EE_all
#if USE_M_1==1
  call MPI_ALLREDUCE(fr,fr_all,ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(FF,FF_all,ncell*ngrp,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  fr=fr_all
  FF=FF_all
#endif
#endif
#if NDUST>0
  call MPI_ALLREDUCE(eps_dust,eps_dust_all,ndust*ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(vdust,vdust_all,ndust*ncell,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)

  eps_dust=eps_dust_all
  vdust=vdust_all
  
#endif

#endif

  if(myid==1)then
     open(172,file='res/ampl.dat')
!!$#if NENER>0
!!$     write(*,116)'lev      x           d          u          v           w          Pnt      P          A          B          C'
!!$#else
!!$     write(*,116)'lev      x           d          u          v           w          P          A          B          C'
!!$#endif
#if NENER>NGRP
     write(*,116)'===================================================================================================================='
     write(*,116)'lev       x            d          u          v          w          Pnt        P          A          B          C      '
#else     
#if USE_FLD==0 && USE_M_1==0
     write(*,116)'===================================================================================================================='
     write(*,116)'lev       x            d          u          v          w          P          A          B          C      '
#endif
#endif     
#if USE_FLD==1
     write(*,116)'========================================================================================================================================'
     write(*,116)'lev       x            d          u          v          w          P          A          B          C       Er_tot    Er(igrp)'
#endif
#if USE_M_1==1
     write(*,116)'============================================================================================================================================================='
     write(*,116)'lev       x            d          u          v          w          P          A          B          C       Er_tot    Er(igrp)  Fr_tot   Fr(igrp)'
#endif
#if NDUST>0
     write(*,116)'===================================================================================================================================='
     write(*,116)'lev       x            d          u          v          w          P          A          B          C      d_dust(idust) vdust'
#endif
     ! Sort radius
     allocate(ind_sort(1:ncell))
     call quick_sort(rr,ind_sort,ncell)
     ! Write results to screen
     smallp=smallc**2/gamma
     nx_loc=icoarse_max-icoarse_min+1
     scale=boxlen/dble(nx_loc)
     ! Prevent underflow for velocity
     do i=1,ncell
        if(ABS(uu(i))<smallc)uu(i)=0.0D0
        if(ABS(vv(i))<smallc)vv(i)=0.0D0
        if(ABS(ww(i))<smallc)ww(i)=0.0D0
        if(ABS(AA(i))<smallc*sqrt(smallr))AA(i)=0.0D0
        if(ABS(BB(i))<smallc*sqrt(smallr))BB(i)=0.0D0
        if(ABS(CC(i))<smallc*sqrt(smallr))CC(i)=0.0D0
     end do
     do i=1,ncell
        ddd=MAX(dd(ind_sort(i)),smallr)
        ppp=MAX((gamma-1.0)*ei(ind_sort(i)),ddd*smallp)
        write(*,*) &
             & ll(ind_sort(i)),  &
             & (rr(i)-dble(icoarse_min))*scale, &
             & ddd , &
             & uu(ind_sort(i)), &
             & vv(ind_sort(i)), &
             & ww(ind_sort(i)), &
#if NENER>NGRP
             & (prad(ind_sort(i),j),j=1,nent), &
#endif
             & ppp, &
             & AA(ind_sort(i)),  &
             & BB(ind_sort(i)),  &
#if USE_FLD==0 && USE_M_1==0
#if NDUST>0
             & CC(ind_sort(i)), &
             & (eps_dust(ind_sort(i),j),j=1,ndust) , &
             & (vdust(ind_sort(i),j),j=1,ndust)
#else
             & CC(ind_sort(i))
#endif
#else
             & CC(ind_sort(i)),  &
#endif
#if USE_FLD==1
             & er(ind_sort(i)),  &
             & (EE(ind_sort(i),j),j=1,ngrp)
#endif
#if USE_M_1==1
             & er(ind_sort(i)),  &
             & (EE(ind_sort(i),j),j=1,ngrp), &
             & fr(ind_sort(i)),  &
             & (FF(ind_sort(i),j),j=1,ngrp)
#endif


     end do
     deallocate(ind_sort)
     write(*,116)'================================================================================================================'

     write(172,*) maxval(abs(uu)),maxval(abs(dd)),maxval(abs(vdust(:,1))),maxval(abs(vdust(:,2))),maxval(abs(eps_dust_all(:,1)/dd)),maxval(abs(eps_dust_all(:,2)/dd)),t
  end if

  ! Deallocate local arrays
  deallocate(mm,rr,dd,dtot,et,ei,ek,em,uu,vv,ww,AA,BB,CC,ll,gg)
  deallocate(mm_all,rr_all,dd_all,dtot_all)
  deallocate(et_all,ei_all,em_all,ek_all)
  deallocate(uu_all,vv_all,ww_all,ll_all,gg_all)
  deallocate(AA_all,BB_all,CC_all)
#if NENER>NGRP
  deallocate(prad,prad_all)
#endif
#if NGRP>0
  deallocate(er,er_all)
#if USE_M_1==1
  deallocate(fr,fr_all)
#endif
#endif
#if NGRP>0
  deallocate(eps_dust,eps_dust_all,vdust,vdust_all)
#endif
  end if
 
#ifndef WITHOUTMPI
  call MPI_BARRIER(MPI_COMM_WORLD,info)
#endif

#endif

111 format(2(1pe12.5,1x))
112 format(i3,1x,1pe10.3,1x,8(1pe10.3,1x))
113 format(i3,1x,1pe12.5,1x,100(1pe10.3,1x))
114 format(' Output ',i5,' cells')
115 format(' Output ',i5,' parts')
116 format(100(A))

end subroutine write_screen
