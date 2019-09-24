module sink_particle_tracer
  use amr_parameters, only : dp, ndim, nvector, verbose, aexp, nlevelmax, boxlen, icoarse_max, icoarse_min, nlevelsheld
  use amr_parameters, only : X_floor
  use amr_commons, only: ncpu, MC_tracer, metal, myid, star, sink, twondim, use_initial_mass, &
  &     write_stellar_densities, communicator
  ! Particle list and stuff
  use pm_commons

  ! Util function
  use utils, only : cgs, cell2grid, grid_level, normalize_position, pi, draw_normal, find_grid_containing_max
  use random, only : ranf, gaussdev




contains

 
  subroutine fill_tracer_comm(ind_part, ind_com, ind_list, np, ilevel, icpu, &
       reception)
    ! Fill the communicator array and remove particles from the linked list
    integer,  intent(in) :: np, ilevel, icpu
    integer, intent(in), dimension(1:nvector) :: ind_part, ind_com, ind_list
    integer :: current_property
    integer :: i, idim
    logical, dimension(1:nvector), save :: ok=.true.

    type(communicator), intent(inout) :: reception(:, :)

    ! Gather particle level and identity
    do i=1,np
       reception(icpu,ilevel)%fp(ind_com(i),2)=levelp(ind_part(i))
       reception(icpu,ilevel)%fp(ind_com(i),3)=idp   (ind_part(i))
    end do

    ! Gather particle position and velocity
    do idim=1,ndim
       do i=1,np
          reception(icpu,ilevel)%up(ind_com(i),idim     )=xp(ind_part(i),idim)
          reception(icpu,ilevel)%up(ind_com(i),idim+ndim)=vp(ind_part(i),idim)
       end do
    end do

    current_property = twondim+1
    ! Gather particle mass
    do i=1,np
       reception(icpu,ilevel)%up(ind_com(i),current_property)=mp(ind_part(i))
    end do
    current_property = current_property+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
    ! Gather particle potential
    do i=1,np
       reception(icpu,ilevel)%up(ind_com(i),current_property)=ptcl_phi(ind_part(i))
    end do
    current_property = current_property+1
#endif

    ! Gather particle birth epoch
    if(star.or.sink)then
       do i=1,np
          reception(icpu,ilevel)%up(ind_com(i),current_property)=tp(ind_part(i))
       end do
       current_property = current_property+1
       if(metal)then
          do i=1,np
             reception(icpu,ilevel)%up(ind_com(i),current_property)=zp(ind_part(i))
          end do
          current_property = current_property+1
       end if
       if(write_stellar_densities) then
          do i=1,np
             reception(icpu,ilevel)%up(ind_com(i),current_property)  =st_n_tp(ind_part(i))
             reception(icpu,ilevel)%up(ind_com(i),current_property+1)=st_n_SN(ind_part(i))
             reception(icpu,ilevel)%up(ind_com(i),current_property+2)=st_e_SN(ind_part(i))
          end do
          current_property = current_property+3
       endif
       if(use_initial_mass)then
          do i=1,np
             reception(icpu,ilevel)%up(ind_com(i),current_property)=mp0(ind_part(i))
          end do
          current_property = current_property+1
       end if
#ifdef NTRACEGROUPS
       do i=1,np
          reception(icpu,ilevel)%up(ind_com(i),current_property)=ptracegroup(ind_part(i))
       end do
       current_property = current_property + 1
#endif
    end if
    ! Family
    do i=1,np
       reception(icpu,ilevel)%fp(ind_com(i), 4) = part2int(typep(ind_part(i)))
    end do

    if (MC_tracer) then
       do i=1,np
          reception(icpu, ilevel)%fp(ind_com(i), 5) = &
               partp(ind_part(i))
       end do
       ! Only useful if you add extra properties
       current_property = current_property + 2
    end if


    ! Remove particles
    call remove_list(ind_part, ind_list, ok, np)
    call add_free(ind_part, np)

  end subroutine fill_tracer_comm

  subroutine empty_tracer_comm(ind_com, np, ilevel, icpu, emission)
    ! Empty communicators for tracer particles and attach particles to the grid
#ifdef NTRACEGROUPS
    use pm_commons, only: ptracegroup
#endif
    integer, intent(in) ::np, icpu, ilevel
    integer, intent(in), dimension(1:nvector)::ind_com
    type(communicator), intent(in) :: emission(:, :)

    integer :: i, idim
    integer, dimension(1:nvector), save :: ind_list, ind_part, ind_cell, ind_level
    logical, dimension(1:nvector), save :: ok=.true.
    integer :: current_property
    real(dp), dimension(1:nvector, 1:ndim), save :: pos

    ! Compute parent grid index
    do idim = 1, 3
       do i = 1, np
          pos(i, idim) = emission(icpu, ilevel)%up(ind_com(i), idim)
       end do
    end do

    ! Particles will be attached at most as deep as the current level
    call find_grid_containing_max(pos, ind_list, ind_cell, ind_level, np, ilevel)

    ! Add particle to parent linked list
    call remove_free(ind_part, np)
    call add_list(ind_part, ind_list, ok, np)

    ! Scatter particle level and identity
    do i=1,np
       levelp(ind_part(i)) = emission(icpu, ilevel)%fp(ind_com(i), 2)
       idp   (ind_part(i)) = emission(icpu, ilevel)%fp(ind_com(i), 3)
    end do

    ! Scatter particle position and velocity
    do idim = 1, ndim
       do i = 1, np
          xp(ind_part(i), idim) = emission(icpu, ilevel)%up(ind_com(i), idim     )
          vp(ind_part(i), idim) = emission(icpu, ilevel)%up(ind_com(i), idim+ndim)
       end do
    end do

    current_property = twondim + 1

    ! Scatter particle mass
    do i = 1, np
       mp(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
    end do
    current_property = current_property + 1

#ifdef OUTPUT_PARTICLE_POTENTIAL
    ! Scatter particle phi
    do i = 1, np
       ptcl_phi(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
    end do
    current_property = current_property+1
#endif

    ! Scatter particle birth epoch
    if(star.or.sink)then
       do i = 1, np
          tp(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
       end do
       current_property = current_property+1
       if(metal)then
          do i = 1, np
             zp(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
          end do
          current_property = current_property+1
       end if
       if(write_stellar_densities) then
          do i = 1, np
             st_n_tp(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)   !SD
             st_n_SN(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property+1) !SD
             st_e_SN(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property+2) !SD
          end do
          current_property = current_property+3
       endif
       if(use_initial_mass)then
          do i = 1, np
             mp0(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
          end do
          current_property = current_property+1
       end if
#ifdef NTRACEGROUPS
       do i = 1, np
          ptracegroup(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
       end do
       current_property = current_property+1
#endif
    end if
    ! Add family
    do i = 1, np
       typep(ind_part(i)) = int2part(emission(icpu, ilevel)%fp(ind_com(i), 4))
       typep(ind_part(i))%tag = typep(ind_part(i))%tag + 10
    end do
    ! MC Tracer
    if (MC_tracer) then
       do i = 1, np
          partp(ind_part(i)) = emission(icpu, ilevel)%fp(ind_com(i), 5)
          move_flag(ind_part(i)) = 1
       end do
    end if
  end subroutine empty_tracer_comm

end module sink_particle_tracer

subroutine tracer2othersink(ind_tracer, isink_new_part, xsink_loc, np)
  ! Transfer tracers from one sink to another one
  use pm_commons
  implicit none

  integer, dimension(1:nvector), intent(in) :: ind_tracer, isink_new_part
  integer, intent(in) :: np
  real(dp), dimension(1:nvector, 1:ndim), intent(in) :: xsink_loc

  integer :: j,  isink_new
  do j = 1, np
     ! Retrieve the index of the new sink
     isink_new = isink_new_part(j)

     ! If the sink has changed, reattach to the new one
     if (isink_new /= 0) then
        partp(ind_tracer(j)) = isink_new

        ! Do not move tracers there, they will be moved in move_fine.f90.
        ! xp(ind_tracer(j), :) = xsink(isink_new, :)
     end if
  end do

end subroutine tracer2othersink

subroutine tracer2sink(ind_tracer, proba, xsink_loc, isink, nattach, dx_loc)
  ! Attach tracers to a sink (accretion onto the BH)
  use amr_commons
  use random, only : ranf
  use pm_commons, only : tracer_seed, typep, partp, move_flag, FAM_TRACER_CLOUD
  implicit none

  integer, intent(in) :: nattach
  integer, dimension(1:nvector), intent(in) :: ind_tracer, isink
  real(dp), dimension(1:nvector), intent(in) :: proba
  real(dp), dimension(1:nvector, 1:3), intent(in) :: xsink_loc
  real(dp), intent(in) :: dx_loc

  logical, dimension(1:nvector), save :: attach = .false.
  integer :: i
  real(dp) :: r

  do i = 1, nattach
     call ranf(tracer_seed, r)
     attach(i) = r < proba(i)
  end do

  ! Change particle pointer and kind
  do i = 1, nattach
     if (attach(i)) then
        partp(ind_tracer(i)) = isink(i)
        typep(ind_tracer(i))%family = FAM_TRACER_CLOUD
     end if
  end do

  ! Do not move tracers there, they will be moved in move_fine.f90 specifically.

  ! do idim = 1, ndim
  !    do i = 1, nattach
  !       if (attach(i)) then
  !          ! Move at most 2*dx in the direction
  !          xtmp = xp(ind_tracer(i), idim)
  !          if (abs(xsink_loc(i, idim) - xtmp) > 2*dx_loc) then
  !             xp(ind_tracer(i), idim) = xtmp + (xsink_loc(i, idim)-xtmp) / dx_loc * 2
  !          else
  !             xp(ind_tracer(i), idim) = xsink_loc(i, idim)
  !          end if
  !       end if
  !    end do
  ! end do

  ! Save state of the particle
  do i = 1, nattach
     if (attach(i)) then
        ! Set to 1 to prevent further moves
        move_flag(ind_tracer(i)) = 1
     end if
  end do

end subroutine tracer2sink
