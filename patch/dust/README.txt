This patch was develloped by Ugo Lebreuilly to treat dust dynamics in the diffusion approximation.


Instructions :

1)It is valid in the diffusion and terminal velocity approximation. i.e when tstop << tdyn.
2)It uses an explicit scheme but it is supposed to be stable as long as the previous approximation is physically acceptable.
3)Dust ratio is treated as a passive scalar contained in uold and unew for ivar in [ firstindex_ndust, firstindex_ndust + NDUST].
4) You can run it if you take a Ndust > 0 in the Makefile, and you must specify dust_diffusion =.true. in the namelist
5) The dust region variable is dust_region, it is a 2D matrix of dim (n_regions, ndust)
6) The dust boundary  variable is dust_bound, it is also a 2D matrix of dim (n_regions, ndust)
7) Several parameters like the grain size (grain_size) or the grain density  (grain_dens) can be specified in the namelist in HYDRO_PARAMS they must be chosen wisely (cf 1 and 2)
8) It is now possible to do multigrain simulation, the two previous parameters are vectors of dimension ndust
9) It is possible to use a simple stopping time, however non-physical for test with K_drag=.true., the drag coefficient K_dust replace grain_size and grain_dens
