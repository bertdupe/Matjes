set(main_f90 ${main_f90}
src/Energy/anisotropy_heisenberg.f90
src/Energy/coupling_ME_J.f90
src/Energy/coupling_ME_D.f90
src/Energy/dipolar_direct.f90
src/Energy/dipolar_fft.f90
src/Energy/dipolar_util.f90
src/Energy/demag-tensor.f90
src/Energy/energy_distrib.f90
src/Energy/energy_output_contribution.f90
src/Energy/exchange_Heisenberg_J.f90
src/Energy/exchange_Heisenberg_D.f90
src/Energy/exchange_TJ.f90
src/Energy/setH_util.f90
src/Energy/set_Hamiltonians.f90
src/Energy/set_fft_Hamiltonians.f90
src/Energy/zeeman.f90
src/Energy/harmonic_phonon.f90
src/Energy/stark.f90
src/Energy/ASR_phonon.f90
src/Energy/mag_biquadratic.f90
src/Energy/forces_from_file.f90
src/Energy/4spin.f90
#src/Energy/demag-matrix.f90   #requires m_fft (not used?)
)
