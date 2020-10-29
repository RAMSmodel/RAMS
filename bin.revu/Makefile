#Makefile
###############################################################################
# Include definitions
###############################################################################
include ../include.mk
###############################################################################
# Compiler commands
###############################################################################
INCLUDES  = $(UTILS_INCS) $(PAR_INCS) $(HDF5_INCS)
F_CD1 = $(F_COMP) -c $(F_OPTS1) $(INCLUDES)
F_CD2 = $(F_COMP) -c $(F_OPTS2) $(INCLUDES)
C_CD  = $(C_COMP) -c $(C_OPTS)  $(INCLUDES) $(PAR_DEFS) $(HDF5_DEFS)

C_CM =$(C_CD)  -D$(CMACH) $< && $(ARCH) $@ $(<F:.c=.o)   && rm -f $(<F:.c=.o)
F_CM1=$(F_CD1) $<            && $(ARCH) $@ $(<F:.f90=.o) && rm -f $(<F:.f90=.o)
F_CM2=$(F_CD2) $<            && $(ARCH) $@ $(<F:.f90=.o) && rm -f $(<F:.f90=.o)
F_CM3=$(F_CD2) -D$(CMACH) $< && $(ARCH) $@ $(<F:.F90=.o) && rm -f $(<F:.F90=.o)

################################################################################
## File extension rules
################################################################################
$(ARC)(%.o): %.f90; $(F_CM2)

$(ARC)(%.o): %.F90; $(F_CM3)

$(ARC)(%.o): %.c;   $(C_CM)

################################################################################
# Define objects
################################################################################
OBJ = $(ARC)($(MODEL)/lib/an_header.o)  \
      $(ARC)($(MODEL)/lib/rconstants.o) \
      $(ARC)($(MODEL)/memory/grid_dims.o) \
      $(ARC)($(MODEL)/micro_bin/micro_prm.o) \
      $(ARC)($(MODEL)/micro/micphys.o) \
      $(ARC)($(MODEL)/kpp/kpp_params.o) \
      $(ARC)($(MODEL)/io/io_params.o) \
      $(ARC)($(MODEL)/core/ref_sounding.o) \
      $(ARC)($(MODEL)/memory/var_tables.o) \
      $(ARC)($(MODEL)/memory/mem_leaf.o) \
      $(ARC)($(MODEL)/memory/mem_scratch.o) \
      $(ARC)($(MODEL)/memory/mem_grid.o) \
      $(ARC)($(MODEL)/memory/mem_varinit.o) \
      $(ARC)($(MODEL)/memory/vtab_fill.o) \
      $(ARC)($(MODEL)/memory/grid_struct.o) \
      $(ARC)($(MODEL)/memory/mem_kpp.o) \
      $(ARC)($(MODEL)/mpi/node_mod.o) \
      $(ARC)($(MODEL)/surface/leaf_coms.o) \
      $(ARC)($(MODEL)/isan/isan_coms.o) \
      $(ARC)($(MODEL)/lib/hdf5_f2c.o) \
      $(ARC)($(MODEL)/lib/hdf5_utils.o) \
      $(ARC)($(MODEL)/lib/charutils.o) \
      $(ARC)($(MODEL)/lib/dateutils.o) \
      $(ARC)($(MODEL)/lib/interp_lib.o) \
      $(ARC)($(MODEL)/lib/map_proj.o) \
      $(ARC)($(MODEL)/lib/numutils.o) \
      $(ARC)($(MODEL)/lib/polarst.o) \
      $(ARC)($(MODEL)/lib/rnamel.o) \
      $(ARC)($(MODEL)/lib/therm_lib.o) \
      $(ARC)($(MODEL)/lib/utils_f.o) \
      $(ARC)($(MODEL)/lib/rget.o) \
      $(ARC)($(MODEL)/lib/rsys.o) \
      $(ARC)($(MODEL)/lib/parlib.o) \
      $(ARC)($(MODEL)/lib/utils_c.o) \
      $(ARC)($(MODEL)/micro/mic_gamma.o) \
      $(ARC)($(MODEL)/revu/revu_utils.o) \
      $(ARC)($(MODEL)/revu/rhdf5_utils.o) \
      $(ARC)($(MODEL)/revu/rhdf5_f2c.o) \
      $(ARC)($(MODEL)/io/rcio.o) \
      $(ARC)($(MODEL)/revu/rcommons.o) \
      $(ARC)($(MODEL)/revu/rbat.o) \
      $(ARC)($(MODEL)/revu/getvar.o) \
      $(ARC)($(MODEL)/revu/textout.o) \
      $(ARC)($(MODEL)/revu/hvlib.o) \
      $(ARC)($(MODEL)/revu/iplt.o) \
      $(ARC)($(MODEL)/revu/rainit.o) \
      $(ARC)($(MODEL)/revu/rcomp.o) \
      $(ARC)($(MODEL)/revu/revugrads.o) \
      $(ARC)($(MODEL)/revu/vplt.o) \
      $(ARC)($(MODEL)/revu/RAMS_hdf5.o)

###############################################################################
# Define archive and executable names
###############################################################################
BASE=revu
EXE=$(BASE)-$(RAMS_VERSION)
ARC=$(BASE)-$(RAMS_VERSION).a

# Define main source.
MAIN_OBJ = ./revu_main.o
MAIN = $(MODEL)/revu/revu_main.f90

# Define targets.
all: $(EXE)

$(EXE): $(ARC) $(MAIN) FORCE
	@echo ""
	$(F_COMP) -o $(EXE) $(MAIN_OBJ) $(LOADER_OPTS) $(ARC) \
	$(HDF5_LIBS) $(LIBS) $(PAR_LIBS)
	rm -f *.o
	@echo ""
	@echo Finished building === $(EXE)
	@echo ""

$(MAIN): FORCE
	@echo ""
	$(F_CD2) $@

$(ARC): $(OBJ)

FORCE:

check: FORCE
	@echo ""
	check

clean:
	@echo ""
	rm -f $(ARC) $(EXE) $(BASE) *.o *.mod *.f
	@echo ""

