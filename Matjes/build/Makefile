include ../make.sys

.SUFFIXES: .o .f90 .F90 .f

OBJS = ext_param_rw.o energy_matrix.o simu-parameters.o basic_types.o \
 io_files_utils.o solver.o setup-Hamiltonian.o rw_motif.o check_inputs.o \
 tightbinding.o Beff_sd.o choosespin.o entropic.o topohall.o topocharge_local.o topocharge_all.o \
 fit.o derived_types.o constant.o vector.o sym_utils.o Correlation.o main.o rw_lattice.o \
 MC_fft.o sampling.o Energy.o relaxationtype.o topocharge.o \
 create_lattice.o   fft.o init_rand_seed.o gauge_util.o Energy_density.o coordinatesv4.o setup_neigh.o \
 qorienplot.o   CalculateAverages.o WriteSpinAndCorrFile.o topocharge_sd.o \
 SkX_utils.o  CreateSpinFile.o spstm.o inp_rw.o rw_dyna.o av_sd.o randist.o spindynamics.o setup_DM.o \
 mtprng.o  MCstep.o InitSpin.o Table_dist.o total_Energy.o DeriveValue.o Relaxation.o topoplot.o \
 setup-simu.o arrange_neihg.o mapping.o indexation.o order_lattice.o user_def_struct.o qorien.o \
 stdtypes.o Efield_sd.o error.o setup_zdir.o structure.o rw_efield.o \
 rw_zdir.o MonteCarlo.o check_restart.o store_relaxation.o \
 initialize_temperature.o cal_local_energy.o welcome.o parallel_tempering.o paratemp.o field-eff_sd.o mpi_prop.o \
 initial-configuration.o init-spiral.o init-DW.o init-skyrmions.o init-Sk-lattice.o init-heavyside.o get-position.o \
 IO-modes.o io-utils.o get-random.o simu-parameters.o init-variables.o convert.o kmesh.o dyna-utils.o pointgrp.o

Matjes: $(OBJS)
	$(LD) -o $@ $^ $(LDFLAGS) $(LIBRARIES)

%.o: %.f90
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) -o $@
	
%.o: %.F90
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) $(INCLUDES) -o $@

%.o: %.f
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) $(INCLUDES) -o $@
	
.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

.F.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o
	
clean:
	$(RM) $(OBJS) *.mod core *.f90 *.F90 *.f

ext_param_rw.o: derived_types.o io-utils.o io_files_utils.o
energy_matrix.o: basic_types.o
basic_types.o:
pointgrp.o: basic_types.o
dyna-utils.o:
kmesh.o: io_files_utils.o
convert.o:
simu-parameters.o: basic_types.o
get-random.o:
IO-modes.o: derived_types.o
get-position.o: derived_types.o
init-DW.o:
init-skyrmions.o: io-utils.o vector.o constant.o create_lattice.o
init-Sk-lattice.o: vector.o init-skyrmions.o
init-heavyside.o:
init-spiral.o: vector.o get-position.o io-utils.o
initial-configuration.o: init-spiral.o init-DW.o init-skyrmions.o init-Sk-lattice.o init-heavyside.o io-utils.o
simu-parameters.o: basic_types.o
io_files_utils.o: derived_types.o
io-utils.o: constant.o convert.o derived_types.o simu-parameters.o vector.o
setup-Hamiltonian.o: derived_types.o inp_rw.o
rw_motif.o: derived_types.o structure.o vector.o
check_inputs.o:
structure.o: vector.o
gauge_util.o:
entropic.o: cal_local_energy.o relaxationtype.o sampling.o
mpi_prop.o:
constant.o:
field-eff_sd.o: inp_rw.o sym_utils.o rw_efield.o
Beff_sd.o: field-eff_sd.o
welcome.o:
CreateSpinFile.o: get-position.o constant.o
topocharge_local.o: vector.o create_lattice.o inp_rw.o rw_lattice.o
topohall.o: vector.o constant.o create_lattice.o rw_lattice.o
fit.o:
user_def_struct.o: vector.o inp_rw.o create_lattice.o
derived_types.o: basic_types.o
rw_zdir.o:
vector.o:
rw_efield.o: create_lattice.o vector.o
topoplot.o: rw_lattice.o inp_rw.o constant.o create_lattice.o
av_sd.o: rw_lattice.o create_lattice.o Energy.o vector.o
WriteSpinAndCorrFile.o: CreateSpinFile.o rw_lattice.o create_lattice.o constant.o
init_rand_seed.o: main.o
MakeLattice.o: main.o
relaxationtype.o:
fft.o: kmesh.o rw_lattice.o create_lattice.o vector.o constant.o
randist.o: rw_dyna.o mtprng.o
Energy_density.o: main.o Energy.o MakeLattice.o rw_lattice.o create_lattice.o inp_rw.o Table_dist.o
rw_dyna.o: main.o
DM_vector.o: main.o
topocharge_sd.o: vector.o rw_lattice.o create_lattice.o constant.o inp_rw.o create_lattice.o
coordinatesv4.o: rw_lattice.o create_lattice.o gauge_util.o
sampling.o: vector.o
InitSpin.o: get-random.o initial-configuration.o IO-modes.o get-position.o create_lattice.o derived_types.o rw_lattice.o order_lattice.o mtprng.o
createSpinFile.o: rw_lattice.o constant.o inp_rw.o
computeDM.o: main.o
#field_sd.o: rw_lattice.o create_lattice.o
rw_lattice.o: vector.o derived_types.o create_lattice.o constant.o
MC_fft.o: vector.o constant.o rw_lattice.o create_lattice.o
Table_dist.o: vector.o derived_types.o
setup_neigh.o: rw_lattice.o create_lattice.o vector.o
spstm.o: rw_lattice.o create_lattice.o vector.o
setup_DM.o: vector.o derived_types.o Table_dist.o
create_lattice.o: derived_types.o
error.o:
mapping.o: create_lattice.o inp_rw.o vector.o derived_types.o
indexation.o: derived_types.o vector.o
#dispBZ.o: inp_rw.o constant.o rw_lattice.o vector.o create_lattice.o constant.o sym_utils.o
order_lattice.o: derived_types.o
pathImc.o:
qorienplot.o: rw_lattice.o inp_rw.o constant.o create_lattice.o
SkX_utils.o: rw_lattice.o vector.o
stdtypes.o:
mtprng.o: stdtypes.o
Efield_sd.o: rw_lattice.o create_lattice.o inp_rw.o Energy.o
#minimize.o: Beff_sd.o inp_rw.o WriteSpinAndCorrFile.o solver.o av_sd.o
DeriveValue.o: inp_rw.o total_Energy.o
#gneb.o: Energy.o minimize.o gneb_utils.o
DeriveValue.o: inp_rw.o total_Energy.o
MonteCarlo.o: initialize_temperature.o inp_rw.o mtprng.o check_restart.o Correlation.o vector.o store_relaxation.o parallel_tempering.o check_restart.o
Relaxation.o: store_relaxation.o mtprng.o Correlation.o MCstep.o
check_restart.o:
store_relaxation.o: constant.o Correlation.o
topocharge.o: vector.o topoplot.o inp_rw.o
parallel_tempering.o: paratemp.o MCstep.o check_restart.o store_relaxation.o initialize_temperature.o CreateSpinFile.o
paratemp.o: fit.o mtprng.o constant.o inp_rw.o
#gneb_utils.o: WriteSpinAndCorrFile.o Beff_sd.o av_sd.o
initialize_temperature.o:
cal_local_energy.o: Energy.o

arrange_neihg.o: create_lattice.o inp_rw.o sym_utils.o constant.o
MCstep.o: choosespin.o cal_local_energy.o main.o create_lattice.o inp_rw.o rw_lattice.o  Energy.o sampling.o relaxationtype.o mtprng.o
solver.o: constant.o create_lattice.o vector.o mtprng.o randist.o rw_dyna.o
Relaxation1000.o: main.o  inp_rw.o create_lattice.o mtprng.o topocharge_all.o
setup_zdir.o: rw_lattice.o create_lattice.o
inp_rw.o: rw_zdir.o rw_lattice.o constant.o mpi_prop.o
spindynamics.o: dyna-utils.o av_sd.o Efield_sd.o createSpinFile.o topocharge_sd.o rw_lattice.o create_lattice.o inp_rw.o rw_dyna.o vector.o av_sd.o randist.o solver.o Beff_sd.o  mtprng.o cal_local_energy.o
total_Energy.o: vector.o inp_rw.o
qorien.o: vector.o inp_rw.o rw_lattice.o
Energy.o: main.o rw_efield.o vector.o rw_lattice.o create_lattice.o inp_rw.o
topocharge.o: rw_lattice.o create_lattice.o vector.o create_lattice.o inp_rw.o
CalculateAverages.o: vector.o rw_lattice.o create_lattice.o topocharge_all.o
main.o: io-utils.o init-variables.o welcome.o WriteSpinAndCorrFile.o topohall.o fft.o qorien.o topocharge.o vector.o Correlation.o constant.o MC_fft.o create_lattice.o inp_rw.o check_inputs.o rw_lattice.o CalculateAverages.o mtprng.o check_restart.o
Correlation.o: vector.o constant.o derived_types.o create_lattice.o
sym_utils.o: vector.o
topocharge_all.o: rw_lattice.o topocharge_local.o

setup-simu.o: setup-Hamiltonian.o rw_motif.o user_def_struct.o topocharge.o arrange_neihg.o indexation.o mapping.o derived_types.o create_lattice.o inp_rw.o rw_lattice.o setup_DM.o main.o setup_neigh.o sym_utils.o Table_dist.o mtprng.o
Energy.o:
cal_local_energy.o: Energy.o
choosespin.o: mtprng.o
