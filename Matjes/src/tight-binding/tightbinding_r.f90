module m_tightbinding_r
use m_rw_TB, only : TB_params
use m_basic_types, only : vec_point
use m_energy_r, only: get_eigenval_r,get_occupation
use m_fermi, only: calc_fermi 
use m_io_files_utils, only: close_file,open_file_write
implicit none
private
public :: tightbinding_r
contains

subroutine tightbinding_r(dimH,TB_pos_ext,mode_mag)

    integer,intent(in)          ::  dimH
    integer,intent(in)          ::  TB_pos_ext(2)
    type(vec_point),intent(in)  ::  mode_mag(:)
    real(8),allocatable         ::  eigval(:)
    real(8)                     ::  E_f

    real(8),allocatable         ::  occupation(:)
    real(8)                     ::  kt
    integer                     ::  n_cell

    kt=0.01d0
    n_cell=size(mode_mag)
    !diagonalize hamiltonian in real spac
    allocate(eigval(dimH),source=0.0d0)
    Call get_eigenval_r(dimH,TB_pos_ext,eigval,mode_mag)
    Call calc_fermi(eigval, TB_params%N_electrons*n_cell, TB_params%kt, E_f)
    Call write_realarray('eigval.dat',eigval)

    !get the occupation
    allocate(occupation(dimH))
    Call get_occupation(dimH,tb_pos_ext,occupation,mode_mag,E_f) 
    Call write_realarray('occupation.dat',occupation)

end subroutine 


subroutine write_realarray(fname,realarr)
    implicit none
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
