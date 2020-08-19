module m_tightbinding_r
use m_tb_params, only : TB_params
use m_tb_types
use m_basic_types, only : vec_point
use m_energy_r, only: get_eigenval_r,get_eigenvec_r
use m_occupation, only: calc_occupation,calc_occupation_sc
use m_fermi, only: calc_fermi 
use m_dos, only: calc_dos 
use m_dos_sc, only: calc_dos_sc
use m_distribution, only: int_distrib,fermi_distrib,dE_fermi_distrib
implicit none
private
public :: tightbinding_r
contains

subroutine tightbinding_r(h_par,mode_mag)
    type(parameters_TB_Hsolve),intent(in)     ::  h_par
    type(vec_point),intent(in)                ::  mode_mag(:)

    real(8),allocatable         ::  eigval(:)
    complex(8),allocatable      ::  eigvec(:,:)

    real(8)                     ::  E_f
    logical                     :: calc_eigval,calc_eigvec

	procedure(int_distrib),pointer	:: dist_ptr => null()


    calc_eigval=TB_params%flow%dos_r.or.TB_params%flow%occ_r.or.TB_params%flow%spec_r.or.TB_params%flow%fermi_r
    calc_eigvec=(TB_params%flow%dos_r.and.TB_params%is_sc).or.TB_params%flow%occ_r

    if(calc_eigvec)then
        write(*,*) 'start eigenvec_r'
        Call get_eigenvec_r(h_par,eigval,eigvec,mode_mag)
    elseif(calc_eigval)then
        write(*,*) 'start eigenval_r'
        Call get_eigenval_r(h_par,eigval,mode_mag)
    endif

    !write spectrum
    if(TB_params%flow%spec_r)then
        write(*,*) 'start write spectrum'
        Call write_realarray('eigval.dat',eigval)
    endif

    !Calculate Fermi energy (only useful without SC)
    E_f=TB_params%io_ef%E_F_in
    if(TB_params%flow%fermi_r)then
        write(*,*) 'start calculate Fermi energy'
        if(TB_params%is_sc)then
            STOP "calculation of Fermi energy doesn't work when using superconductivity"
        else
            Call calc_fermi(eigval, TB_params%io_EF%N_electrons*h_par%ncell, TB_params%io_ef%kt, E_f)
        endif
    endif
    if(TB_params%flow%dos_r)then
        write(*,*) 'start calculate DOS'
        if(TB_params%is_sc)then
            Call calc_dos_sc(eigval,eigvec,TB_params%io_dos,'dos_r_sc.dat')
        else
            Call calc_dos(eigval,TB_params%io_dos,'dos_r.dat')
        endif
    endif

    if(TB_params%flow%occ_r)then
        write(*,*) 'start calculate occupation'
         !maybe use different smearing than EF input
        dist_ptr=>fermi_distrib
        Call calc_occupation(h_par,eigvec,eigval,E_f,TB_params%io_Ef%kt,'occ.dat',dist_ptr)
        dist_ptr=>dE_fermi_distrib
        Call calc_occupation(h_par,eigvec,eigval,E_f,TB_params%io_Ef%kt,'occ_dE.dat',dist_ptr)
    endif

end subroutine 


subroutine write_realarray(fname,realarr)
    use m_io_files_utils, only: close_file,open_file_write
    real(8),intent(in)          :: realarr(:)
    character(len=*),intent(in) ::  fname
    integer                     :: i,io

    io=open_file_write(fname)
    do i=1,size(realarr)
       write(io,'(E16.8)') realarr(i)
    enddo
    call close_file(fname,io)
end subroutine 

end module
