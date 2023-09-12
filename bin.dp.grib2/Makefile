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
OBJ = $(ARC)($(MODEL)/lib/griber_grb2.o)  \
      $(ARC)($(MODEL)/lib/rget.o)

###############################################################################
# Define archive and executable names
###############################################################################
BASE=dgrib
EXE=$(BASE)-$(RAMS_VERSION)
ARC=$(BASE)-$(RAMS_VERSION).a

# Define main source.
MAIN_OBJ = ./dgrib2_main.o
MAIN = $(MODEL)/dprep/dgrib2_main.f90

# Define targets.
all: wgrib2 $(EXE)

wgrib2:
	@echo "wgrib2 executable does not exist!"
	@echo "Precompile wgrib2 and place executable in this directory!"
	exit 1

$(EXE): $(ARC) $(MAIN) FORCE
	@echo ""
	$(F_COMP) -o $(EXE) $(MAIN_OBJ) $(LOADER_OPTS) $(ARC) \
        $(PAR_LIBS)
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

