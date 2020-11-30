module mpi_util
use mpi_basic
implicit none
contains

subroutine distrib_int_index(com,ind_in,ind_out,displ,cnt)
    !subroutine which takes a mpi_communicator and the start and end index (ind_in) of a contiguous integer loop variable .
    !return the loop start and end index splitting up according to the number of processes.
    !If it is not possible to distribute the loop indices equally it increases the range from the beginning( lower mpi ranks)
    !also return the displacement(displ) and count (cnt) array necessary to Scatterv/GatherV
    class(mpi_type),intent(in)      :: com
    integer,intent(in)              :: ind_in(2)
    integer,intent(out)             :: ind_out(2)
    integer,intent(out)             :: displ(com%Np)
    integer,intent(out)             :: cnt(com%Np)

    integer         :: ind_out_id(2,com%Np)
    integer         :: id_N(com%Np),add_N(2,com%np)
    integer         :: size_int
    integer         :: remain, per, add(2)
    integer         :: i

    do i =1,com%Np
        id_N(i)=i-1
    enddo
    size_int=ind_in(2)-ind_in(1)+1
    per=size_int/com%Np
    remain=mod(size_int,com%Np)
    add_N(1,:)=min(id_N,remain)
    add_N(2,:)=min(id_N+1,remain)
    ind_out_id(1,:)=ind_in(1)+id_N*per+add_N(1,:)
    ind_out_id(2,:)=ind_in(1)+(id_N+1)*per-1+add_N(2,:)

    ind_out=ind_out_id(:,com%id+1)
    displ=ind_out_id(1,:)-1
    cnt=ind_out_id(2,:)-ind_out_id(1,:)+1
end subroutine



end module
