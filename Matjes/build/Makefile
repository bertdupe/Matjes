include ../make.sys

.SUFFIXES: .o .f90 .F90 .f *cpp

OBJS = tbessk.o tbessj.o tbessy.o summerfeld_exp.o random_lib.o TB_types.o rw_TB.o tightbinding.o exchange_tight_binding.o chemical_potential_TB.o \
 total_Hamiltonian_TB.o envelope.o eigen_interface.o eigen_sparse_reduce.o eigen_electron.o matrix_power.o matrix.o Hamiltonian_variables.o modes_variables.o tracker.o solver_order.o \
 null.o zeeman.o anisotropy_heisenberg.o solver_commun.o propagator.o couplage_ME.o exchange_Heisenberg.o total_Heisenberg_Ham.o \
 EMwave.o BTeff_sd.o update_time.o TPulse.o heavyside.o rampe.o topo_commons.o internal_fields.o ext_param_rw.o simu-parameters.o basic_types.o \
 operator_pointer_utils.o io_files_utils.o solver.o rw_motif.o info-dynamics.o symmetry_operators.o precision.o \
 tb_params.o energy_set_real_sc.o energy_set_real_sparse.o dos.o dos_sc.o energy_solve_dense.o\
 Beff_sd.o choosespin.o entropic.o topocharge_local.o topocharge_all.o fermi.o energy_set_real.o energy_k.o energy_r.o tightbinding_r.o tightbinding_k.o wavefunction.o DOS.o bandstructure.o \
 highsymmetry.o fit.o derived_types.o constant.o vector.o sym_utils.o Correlation.o main.o rw_lattice.o minimize.o J_sd_exchange.o \
 sampling.o relaxationtype.o external_fields.o init_constant.o dipolar-energy.o dipolar-field.o \
 create_lattice.o forces.o init_rand_seed.o gauge_util.o derivative.o coordinatesv4.o eigen_val_vec.o invert.o determinant.o \
 qorienplot.o rw_MC.o  CalculateAverages.o WriteSpinAndCorrFile.o topocharge_sd.o shape_excitations.o \
 SkX_utils.o  CreateSpinFile.o spstm.o inp_rw.o rw_dyna.o randist.o spindynamics.o setup_DM.o \
 mtprng.o  MCstep.o InitSpin.o Table_dist.o total_Energy.o Relaxation.o topoplot.o plot_FFT.o sort.o \
 print_Beff.o setup-simu.o arrange_neihg.o mapping.o indexation.o order_lattice.o qorien.o path.o \
 stdtypes.o Energy_distrib.o error.o rw_efield.o energy_commons.o excitations.o rotation.o projection.o spline.o \
 rw_zdir.o MonteCarlo.o check_restart.o measure_temp.o store_relaxation.o torques.o tangent.o gneb.o io-gneb.o gneb_utils.o gneb_parameters.o initialize_path.o \
 initialize_temperature.o cal_local_energy.o user_information.o parallel_tempering.o paratemp.o  mpi_prop.o \
 initial-configuration.o init-spiral.o init_random_config.o init-DW.o init-skyrmions.o init-Sk-lattice.o init-heavyside.o get-position.o \
 IO-modes.o io-utils.o get-random.o simu-parameters.o init-variables.o convert.o kmesh.o dyna-utils.o pointgrp.o fftw.o

Matjes: $(OBJS)
	$(LD) -o $@ $^ $(LDFLAGS) $(LIBRARIES)

%.o: %.f90
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) $(INCLUDES) -o $@
	
%.o: %.F90
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) $(INCLUDES) -o $@

%.o: %.f
	$(FC) -cpp $(CPPS) -c $< $(FFLAGS) $(INCLUDES) -o $@
	
%.o: %.cpp
	$(CXX) -cpp $(CPPS) -c $< $(CXXFLAGS) $(C_INCLUDES) -o $@

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

J_sd_exchange.o: Hamiltonian_variables.o symmetry_operators.o rw_TB.o create_lattice.o
eigen_interface.o:
eigen_sparse_reduce.o:
precision.o:
eigen_val_vec.o:
invert.o:
determinant.o:
TB_types.o:
tb_params.o: rw_TB.o
rw_TB.o: io_files_utils.o io-utils.o convert.o TB_types.o
tightbinding_r.o: tb_params.o io_files_utils.o io-utils.o energy_r.o basic_types.o dos.o fermi.o
tbessk.o:
tbessj.o:
tbessy.o:
summerfeld_exp.o: tbessk.o tbessj.o tbessy.o Table_dist.o symmetry_operators.o convert.o create_lattice.o Hamiltonian_variables.o io_files_utils.o io-utils.o
random_lib.o:
total_Heisenberg_Ham.o: Hamiltonian_variables.o zeeman.o anisotropy_heisenberg.o exchange_Heisenberg.o  
exchange_tight_binding.o: io_files_utils.o io-utils.o convert.o Hamiltonian_variables.o symmetry_operators.o rw_TB.o create_lattice.o
chemical_potential_TB.o: io_files_utils.o io-utils.o convert.o Hamiltonian_variables.o symmetry_operators.o create_lattice.o
total_Hamiltonian_TB.o: J_sd_exchange.o exchange_tight_binding.o chemical_potential_TB.o
fermi.o:  energy_commons.o sort.o
energy_set_real.o: energy_commons.o basic_types.o eigen_interface.o tb_params.o
energy_set_real_sc.o: energy_commons.o basic_types.o eigen_interface.o energy_set_real.o tb_params.o
energy_solve_dense.o:
energy_r.o: basic_types.o fermi.o tb_params.o energy_set_real_sparse.o energy_set_real.o energy_set_real_sc.o dos_sc.o energy_solve_dense.o
tightbinding_k.o: basic_types.o get-position.o derived_types.o energy_k.o highsymmetry.o energy_r.o fermi.o dos.o
tightbinding.o: basic_types.o tb_params.o derived_types.o operator_pointer_utils.o get-position.o bandstructure.o DOS.o wavefunction.o fermi.o energy_set_real.o energy_k.o energy_r.o tightbinding_k.o tightbinding_r.o highsymmetry.o
bandstructure.o: basic_types.o derived_types.o fftw.o energy_commons.o
DOS.o: io_files_utils.o io-utils.o
wavefunction.o: basic_types.o
energy_k.o: invert.o eigen_val_vec.o energy_commons.o fftw.o sort.o energy_r.o
tightbinding.o: basic_types.o derived_types.o operator_pointer_utils.o get-position.o bandstructure.o DOS.o wavefunction.o energy_k.o highsymmetry.o

total_Heisenberg_Ham.o: Hamiltonian_variables.o zeeman.o anisotropy_heisenberg.o exchange_Heisenberg.o  
print_Beff.o: convert.o io_files_utils.o
envelope.o:
sort.o:
matrix_power.o:
matrix.o: basic_types.o
Hamiltonian_variables.o: basic_types.o
modes_variables.o: basic_types.o
tracker.o: invert.o envelope.o sort.o topoplot.o topo_commons.o convert.o io_files_utils.o cal_local_energy.o derived_types.o derivative.o
solver_order.o: convert.o
plot_FFT.o: io-utils.o io_files_utils.o convert.o fftw.o
derivative.o: basic_types.o operator_pointer_utils.o io-utils.o io_files_utils.o
rw_MC.o: constant.o io-utils.o io_files_utils.o
null.o: derived_types.o
dipolar-field.o: fftw.o derived_types.o get-position.o io-utils.o io_files_utils.o basic_types.o vector.o
dipolar-energy.o: dipolar-field.o
path.o: rotation.o vector.o constant.o mtprng.o basic_types.o
spline.o:
rotation.o: vector.o
projection.o:
tangent.o: projection.o
shape_excitations.o: vector.o basic_types.o 
zeeman.o: create_lattice.o symmetry_operators.o io-utils.o io_files_utils.o derived_types.o constant.o
anisotropy_heisenberg.o: create_lattice.o symmetry_operators.o io-utils.o io_files_utils.o derived_types.o
exchange_Heisenberg.o: create_lattice.o symmetry_operators.o io-utils.o io_files_utils.o derived_types.o
couplage_ME.o: create_lattice.o io-utils.o io_files_utils.o symmetry_operators.o derived_types.o
propagator.o: BTeff_sd.o io-utils.o io_files_utils.o constant.o vector.o torques.o
solver_commun.o: solver_order.o update_time.o propagator.o Beff_sd.o solver.o derived_types.o info-dynamics.o
init_constant.o:
symmetry_operators.o: invert.o arrange_neihg.o basic_types.o matrix.o pointgrp.o derived_types.o
info-dynamics.o:
torques.o: derived_types.o vector.o io-utils.o io_files_utils.o
EMwave.o: constant.o vector.o io-utils.o io_files_utils.o
BTeff_sd.o: randist.o random_lib.o
update_time.o: constant.o vector.o io-utils.o io_files_utils.o
TPulse.o: constant.o io-utils.o
init_random_config.o: derived_types.o
heavyside.o: constant.o basic_types.o
rampe.o: constant.o basic_types.o 
excitations.o: shape_excitations.o EMwave.o TPulse.o simu-parameters.o heavyside.o rampe.o derived_types.o io-utils.o io_files_utils.o
external_fields.o: create_lattice.o operator_pointer_utils.o derived_types.o get-position.o
measure_temp.o: constant.o vector.o
operator_pointer_utils.o: Hamiltonian_variables.o create_lattice.o null.o derived_types.o get-position.o modes_variables.o
topo_commons.o: vector.o get-position.o derived_types.o operator_pointer_utils.o
internal_fields.o: energy_commons.o get-position.o derived_types.o operator_pointer_utils.o
energy_commons.o: total_Hamiltonian_TB.o total_Heisenberg_Ham.o setup_DM.o convert.o zeeman.o anisotropy_heisenberg.o exchange_Heisenberg.o constant.o symmetry_operators.o get-position.o couplage_ME.o operator_pointer_utils.o derived_types.o io-utils.o io_files_utils.o
ext_param_rw.o: derived_types.o io-utils.o io_files_utils.o
basic_types.o:
pointgrp.o: io_files_utils.o constant.o vector.o derived_types.o basic_types.o user_information.o
dyna-utils.o:
kmesh.o: io_files_utils.o get-position.o

convert.o:
simu-parameters.o: basic_types.o
get-random.o: mtprng.o
IO-modes.o: derived_types.o
get-position.o: io_files_utils.o derived_types.o
init-DW.o:
init-skyrmions.o: io-utils.o vector.o constant.o create_lattice.o
init-Sk-lattice.o: vector.o init-skyrmions.o
init-heavyside.o:
init-spiral.o: vector.o get-position.o io-utils.o
initial-configuration.o: init_constant.o init_random_config.o init-spiral.o init-DW.o init-skyrmions.o init-Sk-lattice.o init-heavyside.o io-utils.o
simu-parameters.o: basic_types.o
io_files_utils.o: derived_types.o
io-utils.o: constant.o convert.o derived_types.o simu-parameters.o vector.o
rw_motif.o: derived_types.o vector.o
gauge_util.o:
entropic.o: cal_local_energy.o total_Energy.o relaxationtype.o sampling.o
mpi_prop.o:
constant.o:
Beff_sd.o: external_fields.o internal_fields.o
user_information.o:
CreateSpinFile.o: get-position.o constant.o io-utils.o
topocharge_local.o: vector.o create_lattice.o inp_rw.o rw_lattice.o
fit.o:
derived_types.o: basic_types.o
rw_zdir.o:
vector.o:
rw_efield.o: create_lattice.o vector.o
topoplot.o: topo_commons.o rw_lattice.o inp_rw.o constant.o create_lattice.o
WriteSpinAndCorrFile.o: CreateSpinFile.o rw_lattice.o create_lattice.o constant.o
init_rand_seed.o: main.o
MakeLattice.o: main.o
relaxationtype.o:
fftw.o: kmesh.o vector.o constant.o io-utils.o basic_types.o derived_types.o energy_commons.o
randist.o: rw_dyna.o mtprng.o
rw_dyna.o: main.o
DM_vector.o: main.o
topocharge_sd.o: vector.o rw_lattice.o create_lattice.o constant.o inp_rw.o create_lattice.o
coordinatesv4.o: rw_lattice.o create_lattice.o gauge_util.o
sampling.o: vector.o
InitSpin.o: get-random.o initial-configuration.o IO-modes.o get-position.o create_lattice.o derived_types.o rw_lattice.o order_lattice.o mtprng.o
createSpinFile.o: rw_lattice.o constant.o inp_rw.o
computeDM.o: main.o
forces.o: derivative.o
rw_lattice.o: io_files_utils.o io-utils.o vector.o derived_types.o create_lattice.o constant.o
Table_dist.o: convert.o user_information.o vector.o derived_types.o
spstm.o: rw_lattice.o create_lattice.o vector.o
setup_DM.o: sym_utils.o vector.o derived_types.o Table_dist.o pointgrp.o
create_lattice.o: derived_types.o
error.o:
mapping.o: rotation.o create_lattice.o inp_rw.o vector.o derived_types.o
indexation.o: derived_types.o vector.o
order_lattice.o: derived_types.o
pathImc.o:
qorienplot.o: rw_lattice.o inp_rw.o constant.o create_lattice.o
SkX_utils.o: rw_lattice.o vector.o
stdtypes.o:
mtprng.o: stdtypes.o
Energy_distrib.o: energy_commons.o operator_pointer_utils.o convert.o io_files_utils.o derived_types.o matrix.o

minimize.o: Beff_sd.o energy_commons.o internal_fields.o vector.o dyna-utils.o WriteSpinAndCorrFile.o solver.o

gneb.o: energy_commons.o minimize.o io-gneb.o gneb_utils.o gneb_parameters.o initialize_path.o
gneb_parameters.o: io_files_utils.o 
gneb_utils.o: spline.o rotation.o projection.o tangent.o WriteSpinAndCorrFile.o io-gneb.o Beff_sd.o gneb_parameters.o
initialize_path.o: dyna-utils.o io-gneb.o gneb_parameters.o gneb_utils.o
io-gneb.o: WriteSpinAndCorrFile.o CreateSpinFile.o vector.o convert.o io-utils.o io_files_utils.o get-random.o

MonteCarlo.o: rw_MC.o initialize_temperature.o inp_rw.o mtprng.o check_restart.o Correlation.o vector.o store_relaxation.o parallel_tempering.o check_restart.o
Relaxation.o: store_relaxation.o mtprng.o Correlation.o MCstep.o
check_restart.o:
store_relaxation.o: constant.o Correlation.o
parallel_tempering.o: paratemp.o MCstep.o check_restart.o store_relaxation.o initialize_temperature.o CreateSpinFile.o
paratemp.o: fit.o mtprng.o constant.o inp_rw.o
initialize_temperature.o:
cal_local_energy.o: energy_commons.o dipolar-energy.o
arrange_neihg.o: rotation.o matrix_power.o create_lattice.o inp_rw.o sym_utils.o constant.o
MCstep.o: choosespin.o cal_local_energy.o main.o create_lattice.o inp_rw.o rw_lattice.o  sampling.o relaxationtype.o mtprng.o
solver.o: torques.o constant.o create_lattice.o vector.o mtprng.o randist.o rw_dyna.o
Relaxation1000.o: main.o  inp_rw.o create_lattice.o mtprng.o topocharge_all.o
inp_rw.o: rw_zdir.o rw_lattice.o constant.o mpi_prop.o
spindynamics.o: precision.o print_Beff.o tracker.o solver_order.o plot_FFT.o forces.o solver_commun.o info-dynamics.o torques.o BTeff_sd.o update_time.o excitations.o measure_temp.o dyna-utils.o Energy_distrib.o createSpinFile.o topocharge_sd.o rw_lattice.o create_lattice.o inp_rw.o rw_dyna.o vector.o randist.o solver.o Beff_sd.o  mtprng.o cal_local_energy.o
total_Energy.o: energy_commons.o vector.o
qorien.o: vector.o inp_rw.o rw_lattice.o
CalculateAverages.o: topo_commons.o vector.o rw_lattice.o create_lattice.o 
main.o: get-random.o minimize.o tightbinding.o io-utils.o init-variables.o user_information.o WriteSpinAndCorrFile.o qorien.o vector.o Correlation.o constant.o create_lattice.o inp_rw.o rw_lattice.o CalculateAverages.o mtprng.o check_restart.o
Correlation.o: vector.o constant.o derived_types.o create_lattice.o
sym_utils.o: vector.o
topocharge_all.o: rw_lattice.o topocharge_local.o

setup-simu.o: rw_TB.o derivative.o null.o dipolar-energy.o io-utils.o io_files_utils.o topo_commons.o excitations.o external_fields.o get-position.o internal_fields.o energy_commons.o user_information.o rw_motif.o arrange_neihg.o indexation.o mapping.o derived_types.o create_lattice.o inp_rw.o rw_lattice.o setup_DM.o main.o sym_utils.o Table_dist.o mtprng.o
cal_local_energy.o:
choosespin.o: mtprng.o
