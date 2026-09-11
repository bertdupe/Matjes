module m_proba_commun
use m_probability
use m_io_utils
use m_io_files_utils
use m_proba_plot
use m_proba_base
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none

private
public :: select_sampling,alloc_P_distrib


contains

subroutine select_sampling(probability_distrib)
class(proba_data), intent(inout)  :: probability_distrib

integer :: io_input,N_bins
character(len=30) :: Psampling_name,Pplot_name
logical :: use_proba,plot_distrib


N_bins=50
plot_distrib=.false.
use_proba=.false.
Psampling_name='Boltzmann'
Pplot_name='histogram'


io_input=open_file_read('input')

call get_parameter(io_input,'input','Psampling_name',Psampling_name)
call get_parameter(io_input,'input','Pplot_name',Pplot_name)
call get_parameter(io_input,'input','N_bins',N_bins)
call get_parameter(io_input,'input','plot_distrib',plot_distrib)
call get_parameter(io_input,'input','calculate_Pdistrib',use_proba)

call close_file('input',io_input)

if (.not.use_proba) return


probability_distrib%N_bins=N_bins
probability_distrib%io_Pdist=plot_distrib

select case(Psampling_name)
   case('Boltzmann')
     write(output_unit,'(a)') 'Boltzmann sampling selected'
     probability_distrib%sampling_type => Boltzmann
   case('Fermi_dirac')
     write(output_unit,'(a)') 'Fermi_dirac sampling selected'
     probability_distrib%sampling_type => Fermi_dirac
   case('bose_einstein')
     write(output_unit,'(a)') 'bose_einstein sampling selected'
     probability_distrib%sampling_type => bose_einstein
   case default
     STOP 'sampling not correct: Boltzmann Fermi_dirac bose_einstein'
end select

select case(Pplot_name)
   case('histogram')
     write(output_unit,'(a)') 'histogram plot selected'
     probability_distrib%plot_type => histogram
   case default
     STOP 'plot name not correct: histogram'
end select

probability_distrib%is_set=.true.

end subroutine

end module
