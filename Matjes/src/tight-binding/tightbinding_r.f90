module m_tightbinding_r
use m_tb_params, only : TB_params
use m_basic_types, only : vec_point
use m_energy_r, only: get_eigenval_r,calc_occupation,get_eigenvec_r
use m_fermi, only: calc_fermi 
use m_dos, only: calc_dos 
use m_dos_sc, only: calc_dos_sc
implicit none
private
public :: tightbinding_r
contains

subroutine tightbinding_r(dimH,TB_pos_ext,mode_mag)
    integer,intent(in)          ::  dimH
    integer,intent(in)          ::  TB_pos_ext(2)
    type(vec_point),intent(in)  ::  mode_mag(:)
    real(8),allocatable         ::  eigval(:)
    complex(8),allocatable      ::  eigvec(:,:)
    real(8)                     ::  E_f

    real(8),allocatable         ::  occupation(:)
    integer                     ::  n_cell

    integer                     :: io_input
    logical                     :: calc_eigval,calc_eigvec


    calc_eigval=TB_params%flow%dos_r.or.TB_params%flow%occ_r.or.TB_params%flow%spec_r.or.TB_params%flow%fermi_r
    calc_eigvec=(TB_params%flow%dos_r.and.TB_params%io_H%is_sc).or.TB_params%flow%occ_r

    n_cell=size(mode_mag)
    if(calc_eigvec)then
        allocate(eigval(dimH),source=0.0d0)
        allocate(eigvec(dimH,dimH),source=cmplx(0.0d0,0.0d0,8))
        write(*,*) 'get eigenvec_r'
        Call get_eigenvec_r(dimH,TB_pos_ext,eigval,eigvec,mode_mag)
    elseif(calc_eigval)then
        allocate(eigval(dimH),source=0.0d0)
        write(*,*) 'get eigenval_r'
        Call get_eigenval_r(dimH,TB_pos_ext,eigval,mode_mag)
    endif

    !diagonalize hamiltonian in real spac
    if(TB_params%flow%spec_r)then
        Call write_realarray('eigval.dat',eigval)
    endif
    E_f=TB_params%io_ef%E_F_in
    if(TB_params%flow%fermi_r)then
        Call calc_fermi(eigval, TB_params%io_EF%N_electrons*n_cell, TB_params%io_ef%kt, E_f)
    endif
    if(TB_params%flow%dos_r)then
        if(TB_params%io_H%is_sc)then
            Call calc_dos_sc(eigval,eigvec,TB_params%io_dos,'dos_r_sc.dat')
        else
            Call calc_dos(eigval,TB_params%io_dos,'dos_r.dat')
        endif
    endif

    if(TB_params%flow%occ_r)then
        !get the occupation
        write(*,*) 'get occupation'
        allocate(occupation(dimH))
        Call calc_occupation(dimH,eigvec,eigval,E_f,TB_params%io_Ef%kt,occupation) !maybe use different smearing than EF input
        Call write_realarray('occupation.dat',occupation)
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
