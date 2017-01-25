#############################################################################
# If you have problems with this makefile, contact Romain.Teyssier@gmail.com
#############################################################################
# Compilation time parameters
NVECTOR = 32
QUADHILBERT = 16
NPRE = 8

NDIM = 3
SOLVER = mhd
NIMHD= 0

NENT = 0

USE_FLD=0
USE_M_1=0
NGRP = 0

NEXTINCT = 0
NPSCAL = 0

NENER = $(NENT)+$(NGRP)

#for RT
NGROUPS = 3
NIONS = 3

PATCH = ../patch/rt/frig_main
EXEC = ramses

MPI = 1
##COMP = "GNU"
COMP = "INTEL"



ifeq ($(NIMHD),1)
   NPSCAL := $(NPSCAL)+$(NIONS)+4 # for internal energy + current
else
   NPSCAL := $(NPSCAL)+$(NIONS)+1 # for internal energy
endif
ifeq ($(USE_FLD),1)
    NVAR=8+$(NENER)+$(NEXTINCT)+$(NPSCAL)
 else
    NVAR=8+$(NENER)+$(NDIM)*$(NGRP)+$(NEXTINCT)+$(NPSCAL)
endif


##note be vary careful with the variables allocation
## for now it is assumed that
## 1  is density
## 2:4 velocity
## 5 pressure
## 6:8 magnetic field
## 9:9+nener : trapped photons (look at definitions in mhd/read_hydro_param)
## 9+nener9+nener+1:9+nener+nextinct : for extenction (look at definition of firstindex_pscal in mhd/read_hydro_param)
## 9+nener+nextinct+1 : 9+nener+nextinct+1+nions : the ions (rt_init, cooling_fine) 
## note defined through imetal=firstindex_pscal+1 in read_hydro_param then idelay=imetal ixion=idelay ichem=ixion 
## tehn in rt_init nvar_count=ichem-1 and iIons = nvar_count+1 (this is terrible, I know....)
## nvar-4:nvar-1 : current (used in godfine1 if IMHD==1) unew(....,nvar-4+idim)
## nvar-1:nvar : internal energy used in set_uold if (energy_fix)

#if you want to add your own scalar this has to be done between : 
#9+nener+nextinct+1+nions+1 and nvar-4

#then do not forgot that within the code the magnetic variable are dupplicated and stored
#nvar+1:nvar+3

#finally note that ichem, idelay and imetal have not been checked carefully

#############################################################################
GITBRANCH = $(shell git rev-parse --abbrev-ref HEAD)
GITHASH = $(shell git log --pretty=format:'%H' -n 1)
GITREPO = $(shell git config --get remote.origin.url)
BUILDDATE = $(shell date +"%D-%T")
DEFINES = -DNVECTOR=$(NVECTOR) -DNDIM=$(NDIM) -DNPRE=$(NPRE) -DNVAR=$(NVAR) -DNENER=$(NENER) \
          -DSOLVER$(SOLVER) -DNGRP=$(NGRP) -DNPSCAL=$(NPSCAL) -DNIMHD=$(NIMHD) \
	  -DNEXTINCT=$(NEXTINCT) -DRT  -DNIONS=$(NIONS) -DNGROUPS=$(NGROUPS)  #-DQUADHILBERT=$(QUADHILBERT) 
ifneq (,$(filter 1,$(USE_FLD) $(USE_M_1)))
DEFINES += -DUSE_FLD=$(USE_FLD) -DUSE_M_1=$(USE_M_1)
endif
#############################################################################
# Fortran compiler options and directives

# GFORTRAN
ifeq ($(COMP),"GNU")
   ifeq ($(MPI),1)
      F90 = mpif90
   else
      F90 = gfortran
      FFLAGS = -DWITHOUTMPI
   endif
   FFLAGS += -x f95-cpp-input -O3 -frecord-marker=4 -ffree-line-length-none $(DEFINES) #-freal-4-real-8 
#   FFLAGS += -fbounds-check -fbacktrace -Wuninitialized -Wall
#   FFLAGS += -ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
#   FFLAGS += -finit-real=nan
endif

# INTEL
ifeq ($(COMP),"INTEL")
   ifeq ($(MPI),1)
      F90 = mpif90
#      FFLAGS = -cpp -fast $(DEFINES) -DNOSYSTEM
      FFLAGS = -cpp  $(DEFINES) -DNOSYSTEM
   else
      F90 = ifort
      FFLAGS = -cpp $(DEFINES) -DWITHOUTMPI
   endif
#   FFLAGS += -warn all -O0 -g -traceback -fpe0 -ftrapuv -check bounds
   FFLAGS +=  -O1
endif

# INTEL
ifeq ($(COMP),"FTN")
   F90 = ftn
   FFLAGS = -xAVX -g -traceback -fpp -fast $(DEFINES) -DNOSYSTEM #-DRT
#   FFLAGS = -O3 -g -traceback -fpe0 -ftrapuv -cpp $(DEFINES) -DNOSYSTEM #-DRT
endif

# PGF90
ifeq ($(COMP),"PGF90")
   ifeq ($(MPI),1)
      F90 = mpif90
      FFLAGS = -O3 -Mpreprocess $(DEFINES)
   else
      F90 = pgf90
      FFLAGS = -Mpreprocess $(DEFINES) -DWITHOUTMPI
   endif
endif

# TAU
ifeq ($(COMP),"TAU")
   F90 = tau_f90.sh
   FFLAGS = -optKeepFiles -optPreProcess -optCPPOpts=$(DEFINES) -DWITHOUTMPI
endif

# XLF
ifeq ($(COMP),"XLF")
   F90 = xlf
   FFLAGS = -WF,-DNDIM=$(NDIM),-DNPRE=$(NPRE),-DNVAR=$(NVAR),-DSOLVER$(SOLVER),-DWITHOUTMPI -qfree=f90 -qsuffix=f=f90 -qsuffix=cpp=f90
endif

# F90
ifeq ($(COMP),"F90")
   F90 = f90
   FFLAGS = -cpp $(DEFINES) -DWITHOUTMPI
endif

#############################################################################
MOD = mod
#############################################################################
# MPI librairies
LIBMPI = 
#LIBMPI = -lfmpi -lmpi -lelan

# --- CUDA libraries, for Titane ---
LIBCUDA = -L/opt/cuda/lib  -lm -lcuda -lcudart

LIBGRACKLE = -lgrackle -lhdf5 -lz
LIBS = $(LIBMPI)

# Lapack
ifeq ($(USE_M_1),1)
LIBS += -llapack
else
ifneq ($(NGRP),$(filter $(NGRP),0 1))
LIBS += -llapack
endif
endif

#############################################################################
# Sources directories are searched in this exact order
BINDIR = .
#PATHSTRING = $(subst ../,$(BINDIR)/../,$(PATCH):../$(SOLVER):../fld/:../nimhd:../aton:../hydro:../pm:../poisson:../amr)
#VPATH = $(PATHSTRING)$(B
ifeq ($(strip $(PATCH)),)
VPATH = $(BINDIR)/../$(SOLVER):$(BINDIR)/../fld/:$(BINDIR)/../nimhd:$(BINDIR)/../aton:$(BINDIR)/../rt:$(BINDIR)/../hydro:$(BINDIR)/../pm:$(BINDIR)/../poisson:$(BINDIR)/../amr
else
VPATH = $(BINDIR)/$(PATCH):$(BINDIR)/../$(SOLVER):$(BINDIR)/../fld/:$(BINDIR)/../nimhd:$(BINDIR)/../aton:$(BINDIR)/../rt:$(BINDIR)/../hydro:$(BINDIR)/../pm:$(BINDIR)/../poisson:$(BINDIR)/../amr
endif
#############################################################################
# All objects
MODOBJ = amr_parameters.o amr_commons.o random.o pm_parameters.o pm_commons.o poisson_parameters.o \
         poisson_commons.o hydro_parameters.o hydro_commons.o cooling_module.o bisection.o sparse_mat.o \
         clfind_commons.o gadgetreadfile.o write_makefile.o write_patch.o write_gitinfo.o \
         radiation_parameters.o rt_parameters.o rt_hydro_commons.o coolrates_module.o \
         rt_spectra.o  rt_metal_cooling_module.o rt_cooling_module.o rt_flux_module.o feedback_module.o cloud_module.o
ifeq ($(NIMHD),1)
MODOBJ +=  variables_X.o
endif

AMROBJ = read_params.o init_amr.o init_time.o init_refine.o adaptive_loop.o amr_step.o update_time.o \
         output_amr.o flag_utils.o physical_boundaries.o virtual_boundaries.o refine_utils.o nbors_utils.o \
         hilbert.o load_balance.o title.o sort.o cooling_fine.o units.o light_cone.o movie.o
# Particle-Mesh objects
PMOBJ = init_part.o output_part.o rho_fine.o synchro_fine.o move_fine.o newdt_fine.o particle_tree.o \
        add_list.o remove_list.o star_formation.o sink_particle.o feedback.o clump_finder.o clump_merger.o \
        flag_formation_sites.o init_sink.o output_sink.o  stellar_particle.o init_stellar.o output_stellar.o
# Poisson solver objects
POISSONOBJ = init_poisson.o phi_fine_cg.o interpol_phi.o force_fine.o multigrid_coarse.o multigrid_fine_commons.o \
             multigrid_fine_fine.o multigrid_fine_coarse.o gravana.o boundary_potential.o rho_ana.o output_poisson.o
# Hydro objects
HYDROOBJ = init_hydro.o init_flow_fine.o write_screen.o output_hydro.o courant_fine.o godunov_fine.o \
           uplmde.o umuscl.o interpol_hydro.o godunov_utils.o condinit.o hydro_flag.o hydro_boundary.o \
           boundana.o read_hydro_params.o synchro_hydro_fine.o \
           multigroup.o m1_utils.o \
           radiative_transfer_bicg.o radiation_boundary.o
# RT objects
RTOBJ = rt_init_hydro.o rt_init_xion.o rt_init.o rt_init_flow_fine.o rt_output_hydro.o rt_godunov_fine.o \
        rt_interpol_hydro.o rt_godunov_utils.o rt_condinit.o rt_hydro_flag.o rt_hydro_boundary.o rt_boundana.o \
        rt_read_hydro_params.o rt_units.o
PATCHOBJ = cooling_module_frig.o sink_RT_feedback.o ## \
        velocity_fine_frig.o write_cool.o output_clumps.o
ifeq ($(NIMHD),1)
HYDROOBJ += uplmde_sts_ohm.o uplmde_sts_ohm_dtu.o
endif
# Turbulence objects
#MODOBJ += turb_parameters.o turb_commons.o
#TURBOBJ = turb_force_utils.o read_turb_params.o init_turb.o mpi_share_turb_fields.o \
#          turb_next_field.o turb_check_time.o write_turb_fields.o read_turb_fields.o \
#          add_turb_forcing.o write_pymses.o libfftw3.a
# All objects
AMRLIB = $(AMROBJ) $(HYDROOBJ) $(PMOBJ) $(POISSONOBJ) $(RTOBJ) $(PATCHOBJ) # $(TURBOBJ)  
# ATON objects
ATON_MODOBJ = timing.o radiation_commons.o rad_step.o
ATON_OBJ = observe.o init_radiation.o rad_init.o rad_boundary.o rad_stars.o rad_backup.o $(BINDIR)/../aton/atonlib/libaton.a
#############################################################################
all:
	for file in $(MODOBJ); do unset MAKELEVEL ; $(MAKE) $$file; done
	$(MAKE) ramses
ramses:	$(MODOBJ) $(AMRLIB) ramses.o
	$(F90) $(MODOBJ) $(AMRLIB) ramses.o -o $(EXEC)$(NDIM)d $(LIBS)
	rm write_makefile.f90
	rm write_patch.f90
ramses_aton: $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o
	$(F90) $(MODOBJ) $(ATON_MODOBJ) $(AMRLIB) $(ATON_OBJ) ramses.o -o $(EXEC)$(NDIM)d $(LIBS) $(LIBCUDA)
	rm write_makefile.f90
	rm write_patch.f90
#############################################################################
write_gitinfo.o: FORCE
	$(F90) $(FFLAGS) -DPATCH=\'$(PATCH)\' -DGITBRANCH=\'$(GITBRANCH)\' -DGITHASH=\'"$(GITHASH)"\' \
 -DGITREPO=\'$(GITREPO)\' -DBUILDDATE=\'"$(BUILDDATE)"\' -c $(BINDIR)/../amr/write_gitinfo.f90 -o $@		
write_makefile.o: FORCE
	$(BINDIR)/../utils/scripts/cr_write_makefile.sh $(MAKEFILE_LIST)
	$(F90) $(FFLAGS) -c write_makefile.f90 -o $@
write_patch.o: FORCE
	$(BINDIR)/../utils/scripts/cr_write_patch.sh $(PATCH)
	$(F90) $(FFLAGS) -c write_patch.f90 -o $@
%.o:%.f90
	$(F90) $(FFLAGS) -c $^ -o $@
FORCE:
#############################################################################
clean :
	rm *.o *.$(MOD) *__genmod.f90
#############################################################################
