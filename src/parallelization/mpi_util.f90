module mpi_util
use mpi_basic
implicit none
contains

subroutine distrib_int_index(com,ind_in,ind_out)
    !subroutine which takes a mpi_communicator and the start and end index (ind_in) of a contiguous integer loop variable .
    !return the loop start and end index splitting up according to the number of processes.
    !If it is not possible to distribute the loop indices equally it increases the range from the beginning( lower mpi ranks)
    class(mpi_type),intent(in)      :: com
    integer,intent(in)              :: ind_in(2)
    integer,intent(out)             :: ind_out(2)

    integer         :: size_int
    integer         :: remain, per, add(2)

    size_int=ind_in(2)-ind_in(1)+1
    per=size_int/com%Np
    remain=mod(size_int,com%Np)
    add=[min(com%id,remain),min(com%id+1,remain)]
    ind_out(1)=ind_in(1)+com%id*per+add(1)
    ind_out(2)=ind_in(1)+(com%id+1)*per-1+add(2)
end subroutine



end module
