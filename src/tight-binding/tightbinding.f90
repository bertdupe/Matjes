module m_tightbinding
use m_TB_types, only: parameters_TB
use m_derived_types, only : lattice
use m_tightbinding_r, only: tightbinding_r
use m_tightbinding_k, only: tightbinding_k
use mpi_basic, only: mpi_type
use,intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public tightbinding
contains
subroutine tightbinding(lat,comm)
    type(lattice), intent(inout)    :: lat
    type(mpi_type),intent(in)       :: comm
    ! internal parameter
    type(parameters_TB)     :: tb_par

    if(comm%ismas) write(output_unit,'(a)') 'entering into the tight-binding routines'
    Call lat%bcast(comm)

    !read tight-binding io parameter from input and set TB_params(m_tb_params)
    if(comm%ismas)then
        call tb_par%read_file('input')
        Call tb_par%init(lat)
    endif
    Call tb_par%bcast(comm)

    !do real-space tight-binding stuff
    if(comm%ismas)then
        if(tb_par%flow%do_r)then
            if(comm%Np>1.and. comm%ismas) write(error_unit,"(//A)") "WARNING, using real-space tight-binding with MPI, this is not implemented and will result in no speed-up"
            Call tightbinding_r(lat,tb_par)   
        endif
    endif

    !do reciprocal-space tight-binding stuff
    if(tb_par%flow%do_k) Call tightbinding_k(lat,tb_par,comm)

end subroutine tightbinding
end module
