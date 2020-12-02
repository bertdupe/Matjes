module mpi_util
use mpi_basic
use m_mc_exp_val, only: exp_val,measure_gather,measure_scatter
implicit none

interface gatherv
    module procedure gatherv_real
    module procedure measure_gather
end interface gatherv

interface scatterv 
    module procedure measure_scatter
end interface scatterv 

interface scatter 
    module procedure scatter_int
end interface
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

subroutine scatter_int(arr,loc,com)
    use mpi_basic
    integer,intent(in)             :: arr(:)
    class(mpi_type),intent(in)     :: com
    integer,intent(out)            :: loc
#ifdef CPP_MPI    
    integer                        :: ierr

    Call MPI_SCATTER(arr(1),1,MPI_INTEGER,loc,1,MPI_INTEGER,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine gatherv_real(arr,com,displ,cnt)
    use mpi_basic
    real(8),intent(in)             :: arr(:)
    class(mpi_type),intent(in)     :: com
    integer,intent(in)             :: displ(com%Np)
    integer,intent(in)             :: cnt(com%Np)
#ifdef CPP_MPI    
    integer                        :: ierr

    Call MPI_Gatherv(arr(1),cnt(com%id+1),MPI_DOUBLE_PRECISION,arr(1),cnt,displ,MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

end module
