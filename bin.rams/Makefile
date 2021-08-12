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
$(ARC)($(MODEL)/cuparm/kf_cuparm%o): $(MODEL)/cuparm/kf_cuparm%f90; $(F_CM1)

$(ARC)(%.o): %.f90; $(F_CM2)

$(ARC)(%.o): %.F90; $(F_CM3)

$(ARC)(%.o): %.c;   $(C_CM)

################################################################################
# Define objects
################################################################################
OBJ = $(ARC)($(MODEL)/lib/an_header.o) \
      $(ARC)($(MODEL)/lib/rconstants.o) \
      $(ARC)($(MODEL)/memory/grid_dims.o) \
      $(ARC)($(MODEL)/micro_bin/micro_prm.o) \
      $(ARC)($(MODEL)/micro/micphys.o) \
      $(ARC)($(MODEL)/kpp/kpp_params.o) \
      $(ARC)($(MODEL)/io/io_params.o) \
      $(ARC)($(MODEL)/core/ref_sounding.o) \
      $(ARC)($(MODEL)/memory/var_tables.o) \
      $(ARC)($(MODEL)/memory/mem_basic.o) \
      $(ARC)($(MODEL)/memory/mem_cuparm.o) \
      $(ARC)($(MODEL)/memory/mem_leaf.o) \
      $(ARC)($(MODEL)/memory/mem_sib.o) \
      $(ARC)($(MODEL)/memory/mem_micro.o) \
      $(ARC)($(MODEL)/memory/mem_radiate.o) \
      $(ARC)($(MODEL)/memory/mem_tracer.o) \
      $(ARC)($(MODEL)/memory/mem_scratch.o) \
      $(ARC)($(MODEL)/memory/mem_turb.o) \
      $(ARC)($(MODEL)/mpi/node_mod.o) \
      $(ARC)($(MODEL)/memory/mem_tend.o) \
      $(ARC)($(MODEL)/memory/mem_grid.o) \
      $(ARC)($(MODEL)/memory/mem_varinit.o) \
      $(ARC)($(MODEL)/memory/mem_nestb.o) \
      $(ARC)($(MODEL)/memory/mem_mksfc.o) \
      $(ARC)($(MODEL)/memory/mem_oda.o) \
      $(ARC)($(MODEL)/memory/mem_kpp.o) \
      $(ARC)($(MODEL)/memory/mem_all.o) \
      $(ARC)($(MODEL)/micro_bin/module_hujisbm.o) \
      $(ARC)($(MODEL)/memory/vtab_fill.o) \
      $(ARC)($(MODEL)/memory/alloc.o)   \
      $(ARC)($(MODEL)/memory/grid_struct.o) \
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
      $(ARC)($(MODEL)/gdf/gdf_input.o) \
      $(ARC)($(MODEL)/gdf/gdf_read_sfc.o) \
      $(ARC)($(MODEL)/gdf/gdf_read_upa.o) \
      $(ARC)($(MODEL)/lib/parlib.o) \
      $(ARC)($(MODEL)/lib/utils_c.o) \
      $(ARC)($(MODEL)/micro/aero_include.o) \
      $(ARC)($(MODEL)/core/rtimh.o) \
      $(ARC)($(MODEL)/core/rams_model.o) \
      $(ARC)($(MODEL)/core/rthrm.o) \
      $(ARC)($(MODEL)/core/coriolis.o) \
      $(ARC)($(MODEL)/core/raco.o) \
      $(ARC)($(MODEL)/core/radvc.o) \
      $(ARC)($(MODEL)/core/rexev.o) \
      $(ARC)($(MODEL)/core/rtimi.o) \
      $(ARC)($(MODEL)/core/modsched.o) \
      $(ARC)($(MODEL)/core/newbase.o) \
      $(ARC)($(MODEL)/init/rhhi.o)  \
      $(ARC)($(MODEL)/init/rinit.o)  \
      $(ARC)($(MODEL)/init/rams_grid.o) \
      $(ARC)($(MODEL)/init/rdint.o)  \
      $(ARC)($(MODEL)/init/gridset.o)  \
      $(ARC)($(MODEL)/bc/cyclic_mod.o) \
      $(ARC)($(MODEL)/bc/rbnd.o) \
      $(ARC)($(MODEL)/bc/rbnd_nonscalar.o) \
      $(ARC)($(MODEL)/bc/rbnd_trsetsns.o) \
      $(ARC)($(MODEL)/kpp/kpp_init.o) \
      $(ARC)($(MODEL)/kpp/kpp_typexfer.o) \
      $(ARC)($(MODEL)/kpp/kpp_timestep.o) \
      $(ARC)($(MODEL)/kpp/kpp_mix.o) \
      $(ARC)($(MODEL)/kpp/kpp_ocint.o) \
      $(ARC)($(MODEL)/kpp/kpp_coeffs.o) \
      $(ARC)($(MODEL)/kpp/kpp_override.o) \
      $(ARC)($(MODEL)/sib/sib2_init.o) \
      $(ARC)($(MODEL)/sib/sib2_co2.o) \
      $(ARC)($(MODEL)/surface/sfc_driver.o) \
      $(ARC)($(MODEL)/surface/leaf3_init.o) \
      $(ARC)($(MODEL)/surface/leaf3_prep.o) \
      $(ARC)($(MODEL)/surface/leaf3_canopy.o) \
      $(ARC)($(MODEL)/surface/leaf3_slayer.o) \
      $(ARC)($(MODEL)/surface/leaf3_soilveg.o) \
      $(ARC)($(MODEL)/surface/ruser.o) \
      $(ARC)($(MODEL)/mpi/mpass_dtl.o) \
      $(ARC)($(MODEL)/mpi/mpass_init.o) \
      $(ARC)($(MODEL)/mpi/mpass_lbc.o)  \
      $(ARC)($(MODEL)/mpi/mpass_nest.o) \
      $(ARC)($(MODEL)/mpi/mpass_cyclic.o) \
      $(ARC)($(MODEL)/mpi/para_init.o) \
      $(ARC)($(MODEL)/mpi/par_decomp.o) \
      $(ARC)($(MODEL)/mpi/rnode.o) \
      $(ARC)($(MODEL)/nesting/nest_drivers.o) \
      $(ARC)($(MODEL)/nesting/nest_intrp.o) \
      $(ARC)($(MODEL)/nesting/nest_feed.o) \
      $(ARC)($(MODEL)/radiate/rrad3.o) \
      $(ARC)($(MODEL)/radiate/rad_driv.o) \
      $(ARC)($(MODEL)/radiate/rad_aero.o) \
      $(ARC)($(MODEL)/radiate/rad_mclat.o) \
      $(ARC)($(MODEL)/radiate/rad_ccmp.o) \
      $(ARC)($(MODEL)/radiate/rad_stable.o) \
      $(ARC)($(MODEL)/radiate/rrad2.o) \
      $(ARC)($(MODEL)/init/data_init.o) \
      $(ARC)($(MODEL)/micro/mic_aero.o) \
      $(ARC)($(MODEL)/micro/mic_chknan.o) \
      $(ARC)($(MODEL)/micro/mic_init.o) \
      $(ARC)($(MODEL)/micro/mic_coll.o) \
      $(ARC)($(MODEL)/micro/mic_driv.o) \
      $(ARC)($(MODEL)/micro/mic_misc.o) \
      $(ARC)($(MODEL)/micro/mic_adj.o) \
      $(ARC)($(MODEL)/micro/mic_vap.o) \
      $(ARC)($(MODEL)/micro/mic_nuc.o) \
      $(ARC)($(MODEL)/micro/mic_nuctab.o) \
      $(ARC)($(MODEL)/micro/mic_nucpre.o) \
      $(ARC)($(MODEL)/micro/mic_tabs.o) \
      $(ARC)($(MODEL)/micro/mic_gamma.o) \
      $(ARC)($(MODEL)/micro/aero_sources.o) \
      $(ARC)($(MODEL)/micro/aero_deposit.o) \
      $(ARC)($(MODEL)/micro_bin/microphysics.o) \
      $(ARC)($(MODEL)/io/anal_extra.o) \
      $(ARC)($(MODEL)/io/history_start.o) \
      $(ARC)($(MODEL)/io/anal_write.o) \
      $(ARC)($(MODEL)/io/rcio.o) \
      $(ARC)($(MODEL)/io/inithis.o) \
      $(ARC)($(MODEL)/io/opspec.o) \
      $(ARC)($(MODEL)/io/rname.o) \
      $(ARC)($(MODEL)/io/ranlavg.o) \
      $(ARC)($(MODEL)/io/rprnt.o) \
      $(ARC)($(MODEL)/mksfc/sst_read.o) \
      $(ARC)($(MODEL)/mksfc/ndvi_read.o) \
      $(ARC)($(MODEL)/mksfc/mksfc_driver.o) \
      $(ARC)($(MODEL)/mksfc/mksfc_sfc.o) \
      $(ARC)($(MODEL)/mksfc/mksfc_top.o) \
      $(ARC)($(MODEL)/mksfc/mksfc_sst.o) \
      $(ARC)($(MODEL)/mksfc/mksfc_ndvi.o) \
      $(ARC)($(MODEL)/mksfc/geodat.o) \
      $(ARC)($(MODEL)/mksfc/landuse_input.o) \
      $(ARC)($(MODEL)/mksfc/nest_geosst.o) \
      $(ARC)($(MODEL)/fdda/oda_read.o) \
      $(ARC)($(MODEL)/fdda/oda_krig.o) \
      $(ARC)($(MODEL)/fdda/oda_nudge.o) \
      $(ARC)($(MODEL)/fdda/oda_proc_obs.o) \
      $(ARC)($(MODEL)/fdda/oda_sta_count.o) \
      $(ARC)($(MODEL)/fdda/oda_sta_input.o) \
      $(ARC)($(MODEL)/fdda/varf_read.o) \
      $(ARC)($(MODEL)/fdda/varf_update.o) \
      $(ARC)($(MODEL)/fdda/nud_analysis.o) \
      $(ARC)($(MODEL)/turb/turb_k.o) \
      $(ARC)($(MODEL)/turb/turb_ke.o) \
      $(ARC)($(MODEL)/turb/turb_diff.o) \
      $(ARC)($(MODEL)/turb/rgrad.o) \
      $(ARC)($(MODEL)/isan/aobj.o) \
      $(ARC)($(MODEL)/isan/asgen.o) \
      $(ARC)($(MODEL)/isan/asti2.o) \
      $(ARC)($(MODEL)/isan/asti.o) \
      $(ARC)($(MODEL)/isan/astp.o) \
      $(ARC)($(MODEL)/isan/avarf.o) \
      $(ARC)($(MODEL)/isan/file_inv.o) \
      $(ARC)($(MODEL)/isan/isan_name.o) \
      $(ARC)($(MODEL)/isan/isan_io.o) \
      $(ARC)($(MODEL)/isan/refstate.o) \
      $(ARC)($(MODEL)/isan/v_interps.o) \
      $(ARC)($(MODEL)/isan/write_varf.o) \
      $(ARC)($(MODEL)/cuparm/conv_coms.o) \
      $(ARC)($(MODEL)/cuparm/kf_cuparm.o) \
      $(ARC)($(MODEL)/cuparm/kf_rconv.o) \
      $(ARC)($(MODEL)/cuparm/kf_driver.o) \
      $(ARC)($(MODEL)/cuparm/rconv.o) \
      $(ARC)($(MODEL)/micro/write_scm.o)

###############################################################################
# Define archive and executable names
###############################################################################
BASE=rams
EXE=$(BASE)-$(RAMS_VERSION)
ARC=$(BASE)-$(RAMS_VERSION).a

# Define main source.
MAIN_OBJ = ./rammain.o
MAIN = $(MODEL)/core/rammain.f90

# Define targets.

all: $(EXE)

$(EXE): $(ARC) $(MAIN) FORCE
	@echo ""
	$(F_COMP) -o $(EXE) $(MAIN_OBJ) $(LOADER_OPTS) $(ARC) \
	$(HDF5_LIBS) $(LIBS) $(PAR_LIBS)
	rm -f *.o $(BASE)
	ln -s $(EXE) $(BASE)
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
