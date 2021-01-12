module m_tightbinding_r
use m_tb_params, only : TB_params
use m_tb_types
use m_derived_types, only : lattice
use m_occupation, only: calc_occupation
use m_occupation_mult, only: occupation_mult
use m_fermi, only: calc_fermi 
use m_dos, only: calc_dos 
use m_dos_sc, only: calc_dos_sc
use m_distribution, only: int_distrib,fermi_distrib,dE_fermi_distrib
use m_save_state_r,only: TB_write_states_r, TB_read_states_r
use m_init_Hr
use m_H_tb_public
use, intrinsic :: iso_fortran_env, only : output_unit
implicit none
private
public :: tightbinding_r
contains

subroutine tightbinding_r(lat,h_par)
    type(lattice), intent(in) :: lat
    type(parameters_TB_Hsolve),intent(in)     ::  h_par

    real(8),allocatable         ::  eigval(:)
    complex(8),allocatable      ::  eigvec(:,:)

    real(8)                     ::  E_f
    logical                     :: calc_eigval,calc_eigvec

    logical                     :: read_success
    class(H_tb),allocatable     :: H

    procedure(int_distrib),pointer  :: dist_ptr => null()

    !check what to calculate (eigenvalue/eigenvector)
    calc_eigval=TB_params%flow%dos_r.or.TB_params%flow%occ_r.or.TB_params%flow%spec_r.or.TB_params%flow%fermi_r
    calc_eigvec=(TB_params%flow%dos_r.and.TB_params%is_sc).or.TB_params%flow%occ_r

    !possibly read previous solution
    if(TB_params%flow%read_solution_r)then
        Call TB_read_states_r(eigval,eigvec,read_success)
        calc_eigvec=calc_eigvec.and..not.allocated(eigvec)
        calc_eigval=calc_eigval.and..not.allocated(eigval)
    else
        read_success=.false.
    endif

    !initialize Hamiltonian
    if(calc_eigvec.or.calc_eigval) Call get_all_Hr(lat,TB_params%io_H,H)

    !calculate eigenvalue/eigenvector
    if(calc_eigvec)then
        write(output_unit,'(//A)') 'Start calculating eigenvectors'
        Call H%get_evec(eigval,eigvec)
        read_success=.false.
    elseif(calc_eigval)then
        write(output_unit,'(//A)') 'Start calculating eigenvalues'
        Call H%get_eval(eigval)
        read_success=.false.
    endif

    !write eigenvalue/eigenvector to file for later reuse
    if(.not.read_success.and.TB_params%flow%write_solution_r) Call TB_write_states_r(eigval,eigvec)

    !write spectrum
    if(TB_params%flow%spec_r)then
        write(output_unit,'(/A)') 'start write spectrum'
        Call write_realarray('eigval.dat',eigval)
    endif

    !Calculate Fermi energy (only useful without SC)
    E_f=TB_params%io_ef%E_F_in
    if(TB_params%flow%fermi_r)then
        write(output_unit,'(/A)') 'start calculate Fermi energy'
        if(TB_params%is_sc)then
            STOP "calculation of Fermi energy doesn't work when using superconductivity"
        else
            Call calc_fermi(eigval, TB_params%io_EF%N_electrons*h_par%ncell, TB_params%io_ef%kt, E_f)
        endif
    endif

    if(TB_params%flow%dos_r)then
        write(output_unit,'(/A)') 'start calculate DOS'
        if(TB_params%is_sc)then
            Call calc_dos_sc(eigval,eigvec,TB_params%io_dos,'dos_r_sc.dat')
        else
            Call calc_dos(eigval,TB_params%io_dos,'dos_r.dat')
        endif
    endif

    if(TB_params%flow%occ_r)then
        write(output_unit,'(/A)') 'start calculate occupation'
         !maybe use different smearing than EF input
        dist_ptr=>fermi_distrib
        Call calc_occupation(h_par,eigvec,eigval,E_f,TB_params%io_Ef%kt,'occ.dat',dist_ptr)
        dist_ptr=>dE_fermi_distrib
        Call calc_occupation(h_par,eigvec,eigval,E_f,TB_params%io_Ef%kt,'occ_dE.dat',dist_ptr)
    endif

    if(TB_params%flow%occ_mult_r)then
        write(output_unit,'(/A)') 'start calculate multiple occupations'
        Call occupation_mult(h_par,TB_params%io_occ_mult,eigval,eigvec)
    endif

end subroutine 


subroutine write_realarray(fname,realarr)
    use m_io_files_utils, only: close_file,open_file_write
    real(8),intent(in)          :: realarr(:)
    character(len=*),intent(in) ::  fname
    integer                     :: i,io

    io=open_file_write(fname)
    write(io,'(E16.8)') realarr
    call close_file(fname,io)
end subroutine 

end module
