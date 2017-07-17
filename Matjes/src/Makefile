include ../make.sys

.SUFFIXES: .F .f .f90 .F90

DIRSRC= ./SRC
DIRBIN= ../bin

SRCS= MC_fft constant main rw_lattice rw_zdir vector inp_rw Relaxation field-eff setup_zdir 	sampling rw_efield Energy SignatureFile relaxationtype topocharge rw_skyrmion setup_neigh 	mpi_gatherorsend fft gauge_util Energy_density coordinatesv4 structure demag-matrix 	CalculateAverages WriteSpinAndCorrFile CreateSpinFile total_Energy DeriveValue CreationDestructionSkyrmion 	init_rand_seed spstm MCstep InitSpin rw_dyna av_dynamics spindynamics topoplot qorien demag-tensor        dispBZ randist setup-simu setup_DM topocharge_sd field_sd Correlation sym_utils parallel_tempering

CPP_OBJS=

MPI  =

OBJS = pimc.o choosespin.o Torana.o entropic.o topohall.o topocharge_local.o topocharge_all.o calculate_Beff.o fit.o derived_types.o constant.o vector.o sym_utils.o Correlation.o main.o rw_lattice.o 	 field-eff.o MC_fft.o sampling.o Energy.o SignatureFile.o relaxationtype.o topocharge.o rw_skyrmion.o create_lattice.o 	 fft.o init_rand_seed.o gauge_util.o Energy_density.o coordinatesv4.o setup_neigh.o qorienplot.o 	 CalculateAverages.o WriteSpinAndCorrFile.o CreationDestructionSkyrmion.o topocharge_sd.o SkX_utils.o 	 CreateSpinFile.o spstm.o inp_rw.o rw_dyna.o av_dynamics.o randist.o spindynamics.o setup_DM.o mtprng.o 	 MCstep.o InitSpin.o Table_dist.o total_Energy.o DeriveValue.o Relaxation.o topoplot.o dispBZ.o setup-simu.o field_sd.o 	 arrange_neihg.o mapping.o indexation.o order_lattice.o user_def_struct.o qorien.o solver.o stdtypes.o 	 error_correction_SD.o Efield_sd.o error.o setup_zdir.o structure.o rw_efield.o rw_zdir.o gneb.o minimize.o 	 MonteCarlo.o check_restart.o store_relaxation.o gneb_utils.o initialize_temperature.o cal_local_energy.o welcome.o parallel_tempering.o paratemp.o

DIP =

Matjes : $(OBJS) $(MPI) $(DIP)
	$(LD) -o $@ $(LDFLAGS) $(FIRST) $(OBJS) $(DIP) $(MPI) $(LIBRARIES)

.f.o:
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) -o $@

.F.o:
	$(FC) -c $< $(FFLAGS) -o $@

.f77.o:
	$(FC) -c $< $(FFLAGS) -o $@

.f90.o:
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) -o $@

.F90.o:
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) $(INCLUDES) -o $@

.f.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

.F.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

.f77.a:
	$(FC) -c $(FFLAGS) $<
	$(AR) $(ARFLAGS) $@ $*.o
	$(RM) $*.o

clean:
	$(RM) $(FIRST) $(OBJS) $(MPI) $(DIP) *.mod core *.o

touchsrcs:
	touch $(SRCS)

tests:
	cd tests && $(make)

pimc.o: mtprng.o
welcome.o:
CreateSpinFile.o: constant.o rw_lattice.o inp_rw.o
topocharge_local.o: vector.o create_lattice.o inp_rw.o rw_lattice.o
topohall.o: vector.o constant.o create_lattice.o rw_lattice.o
calculate_Beff.o: field-eff.o inp_rw.o create_lattice.o randist.o create_lattice.o mtprng.o
fit.o:
user_def_struct.o: vector.o inp_rw.o create_lattice.o
derived_types.o:
rw_zdir.o:
vector.o:
rw_efield.o: create_lattice.o vector.o
topoplot.o: rw_lattice.o inp_rw.o constant.o create_lattice.o
av_dynamics.o: rw_lattice.o create_lattice.o Energy.o vector.o
WriteSpinAndCorrFile.o: CreateSpinFile.o rw_lattice.o create_lattice.o constant.o
init_rand_seed.o: main.o
CreationDestructionSkyrmion.o: rw_skyrmion.o rw_lattice.o create_lattice.o
MakeLattice.o: main.o
relaxationtype.o: field-eff.o
fft.o: rw_lattice.o create_lattice.o vector.o constant.o
randist.o: rw_dyna.o mtprng.o
Energy_density.o: main.o Energy.o MakeLattice.o rw_lattice.o create_lattice.o inp_rw.o Table_dist.o
rw_dyna.o: main.o
DM_vector.o: main.o
topocharge_sd.o: vector.o rw_lattice.o create_lattice.o constant.o inp_rw.o create_lattice.o
coordinatesv4.o: rw_lattice.o create_lattice.o gauge_util.o
SignatureFile.o: create_lattice.o inp_rw.o rw_lattice.o rw_efield.o
sampling.o: vector.o
InitSpin.o: create_lattice.o derived_types.o rw_lattice.o order_lattice.o mtprng.o
createSpinFile.o: rw_lattice.o constant.o inp_rw.o
computeDM.o: main.o
field_sd.o: field-eff.o rw_lattice.o create_lattice.o
rw_skyrmion.o: main.o SkX_utils.o
rw_lattice.o: vector.o derived_types.o create_lattice.o
MC_fft.o: vector.o constant.o rw_lattice.o create_lattice.o
Table_dist.o: vector.o derived_types.o
setup_neigh.o: rw_lattice.o create_lattice.o vector.o
spstm.o: rw_lattice.o create_lattice.o vector.o
setup_DM.o: vector.o derived_types.o Table_dist.o
create_lattice.o: derived_types.o
error.o:
mapping.o: create_lattice.o inp_rw.o vector.o derived_types.o
indexation.o: derived_types.o vector.o
dispBZ.o: inp_rw.o constant.o rw_lattice.o vector.o create_lattice.o constant.o sym_utils.o
order_lattice.o: derived_types.o
pathImc.o:
qorienplot.o: rw_lattice.o inp_rw.o constant.o create_lattice.o
SkX_utils.o: rw_lattice.o vector.o
stdtypes.o:
mtprng.o: stdtypes.o
Efield_sd.o: rw_lattice.o create_lattice.o inp_rw.o Energy.o
minimize.o: calculate_Beff.o inp_rw.o WriteSpinAndCorrFile.o solver.o av_dynamics.o
DeriveValue.o: inp_rw.o total_energy.o
gneb.o: field-eff.o Energy.o minimize.o gneb_utils.o
minimize.o: calculate_Beff.o inp_rw.o WriteSpinAndCorrFile.o solver.o av_dynamics.o
DeriveValue.o: inp_rw.o total_energy.o
MonteCarlo.o: initialize_temperature.o inp_rw.o mtprng.o check_restart.o Correlation.o vector.o store_relaxation.o parallel_tempering.o check_restart.o
Relaxation.o: store_relaxation.o mtprng.o Correlation.o MCstep.o
check_restart.o:
store_relaxation.o: constant.o Correlation.o
topocharge.o: vector.o topoplot.o inp_rw.o
parallel_tempering.o: paratemp.o MCstep.o check_restart.o store_relaxation.o initialize_temperature.o CreateSpinFile.o
paratemp.o: fit.o mtprng.o constant.o inp_rw.o
gneb_utils.o: WriteSpinAndCorrFile.o calculate_Beff.o av_dynamics.o
initialize_temperature.o:
cal_local_energy.o: Energy.o

arrange_neihg.o: create_lattice.o inp_rw.o sym_utils.o constant.o
MCstep.o: choosespin.o cal_local_energy.o main.o create_lattice.o inp_rw.o rw_lattice.o  Energy.o sampling.o relaxationtype.o mtprng.o
solver.o: constant.o create_lattice.o vector.o mtprng.o randist.o rw_dyna.o
Relaxation1000.o: main.o  inp_rw.o create_lattice.o mtprng.o topocharge_all.o
setup_zdir.o: rw_lattice.o create_lattice.o
inp_rw.o: rw_zdir.o rw_lattice.o constant.o mpi_prop.o
spindynamics.o: Torana.o Efield_sd.o createSpinFile.o error_correction_SD.o topocharge_sd.o rw_lattice.o create_lattice.o inp_rw.o field-eff.o rw_dyna.o vector.o av_dynamics.o randist.o solver.o calculate_Beff.o  mtprng.o cal_local_energy.o
total_energy.o: vector.o inp_rw.o
qorien.o: vector.o inp_rw.o rw_lattice.o
Energy.o: main.o rw_efield.o vector.o rw_lattice.o create_lattice.o inp_rw.o
topocharge.o: rw_lattice.o create_lattice.o vector.o create_lattice.o inp_rw.o
CalculateAverages.o: vector.o rw_lattice.o create_lattice.o topocharge_all.o
main.o: pimc.o welcome.o WriteSpinAndCorrFile.o topohall.o fft.o qorien.o topocharge.o vector.o Correlation.o constant.o MC_fft.o create_lattice.o inp_rw.o rw_lattice.o CalculateAverages.o mtprng.o check_restart.o
Correlation.o: vector.o constant.o derived_types.o create_lattice.o
sym_utils.o: vector.o
error_correction_SD.o: create_lattice.o solver.o vector.o constant.o inp_rw.o
topocharge_all.o: rw_lattice.o topocharge_local.o

setup-simu.o: user_def_struct.o topocharge.o arrange_neihg.o indexation.o mapping.o derived_types.o create_lattice.o inp_rw.o rw_lattice.o setup_DM.o main.o setup_neigh.o sym_utils.o Table_dist.o mtprng.o split_work.o

field-eff.o: inp_rw.o sym_utils.o rw_efield.o
error_correction_SD.o: create_lattice.o solver.o vector.o constant.o inp_rw.o
entropic.o: CalculateAverages.o sampling.o relaxationtype.o inp_rw.o mtprng.o choosespin.o vector.o constant.o cal_local_energy.o
Torana.o: calculate_Beff.o solver.o constant.o
choosespin.o: mtprng.o
