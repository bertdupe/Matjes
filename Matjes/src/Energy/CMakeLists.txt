set(main_f90 ${main_f90}
src/Energy/anisotropy_heisenberg.f90
src/Energy/chemical_potential_TB.f90
src/Energy/couplage_ME.f90
#src/Energy/demag-matrix.f90   #requires m_fft (not used?)
src/Energy/demag-tensor.f90
src/Energy/dipolar-energy.f90
src/Energy/energy_commons.f90
src/Energy/Energy_distrib.f90
src/Energy/exchange_Heisenberg.f90
src/Energy/exchange_tight_binding.f90
src/Energy/four-spin.f90
src/Energy/J_sd_exchange.f90
src/Energy/summerfeld_exp.f90
src/Energy/total_Hamiltonian_TB.f90
src/Energy/total_Heisenberg_Ham.f90
src/Energy/zeeman.f90

src/Energy/ham_public.f90
src/Energy/ham_base.f90
src/Energy/ham_coo.f90
src/Energy/ham_dense.f90
src/Energy/ham_manual.f90
src/Energy/ham_sparse_mkl.f90
src/Energy/ham_set_util.f90
)
