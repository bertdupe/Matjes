module m_simu_parameters
use m_basic_types
!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  A module that contains the different parameters to run the simulations
!  Each paramater is defined as a character string in a specific variable
!  The variables are defined depending on the type of options that they trigger
!
!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Basic type to store strings with different length
!
!!!!!!!!!!!!!!!!!!!!!!!!!


type(var_name),dimension(*),parameter :: type_simu= [ var_name('magnet-dynamics'),           &
                                                      var_name('spin-dynamics'),           &
                                                      var_name('metropolis'),        &
                                                      var_name('GNEB'),              &
                                                      var_name('parallel-tempering'),          &
                                                      var_name('minimization'),      &
                                                      var_name('entropic'),          &
                                                      var_name('tight-binding'),      &
                                                      var_name('molecular_dynamics'),      &
                                                      var_name('minimize_infdamp') ,	&
                                                      var_name('wavefunc_eval') ,	&
                                                      var_name('llg_diag')]

type(var_name),dimension(4),parameter :: type_excitations= [ var_name('rampe'),           &
                                                      var_name('heavyside'),          &
                                                      var_name('TPulse'),             &
                                                      var_name('EMwave') ]

end module m_simu_parameters
