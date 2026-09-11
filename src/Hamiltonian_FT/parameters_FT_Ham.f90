module m_parameters_FT_Ham
implicit none

type parameters_FT_HAM_IO
    integer             ::  dimH=-1         !final size of Hamiltonian including all modifications
    !solving parameters
    integer             ::  i_diag=2  !different diagonalization methods
    logical             ::  sparse=.false.  !do calculation sparse
    real(8)             ::  Ebnd(2)=[-1.0d+99,1.0d+99]     !minimal and maximal energy values to consider in restricted eigensolver routines
    integer             ::  estNe=0                       !estimated number of eigenvalues in interval
    real(8)             ::  diag_acc=1.0d-12    ! accuracy of iterative eigenvalue solution (so far only fpm input)
contains
    procedure   :: read_file
end type

private
public parameters_FT_HAM_IO

contains

subroutine read_file(this,io,fname)
    use m_io_read_util
    use m_io_utils
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    class(parameters_FT_HAM_IO),intent(inout)  :: this
    integer,intent(in)                         :: io
    character(len=*), intent(in)               :: fname
    ! Internal variables
    logical :: test

    call get_parameter(io,fname,'FT_sparse',        this%sparse)
    call get_parameter(io,fname,'FT_diag',          this%i_diag)
    call get_parameter(io,fname,'FT_diag_acc',      this%diag_acc)
    call get_parameter(io,fname,'FT_diag_Ebnd',     this%Ebnd)
    call get_parameter(io,fname,'FT_diag_Emin',     this%Ebnd(1))
    call get_parameter(io,fname,'FT_diag_Emax',     this%Ebnd(2))
    call get_parameter(io,fname,'FT_diag_estNe',    this%estNe)
    if(this%Ebnd(1)>=this%Ebnd(2))then
        write(error_unit,'(2/A/2(E16.8/))') "WARNING, FT diagonalization minimal energy bound is smaller than maximal energy bound:", this%Ebnd
        STOP "Fix input"
    endif
end subroutine

end module
