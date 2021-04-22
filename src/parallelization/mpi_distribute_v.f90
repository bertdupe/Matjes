module mpi_distrib_v
use mpi_basic
use m_mc_exp_val, only: exp_val,measure_gatherv,measure_scatterv
use m_mc_therm_val, only: thermo_gatherv
implicit none

interface gatherv
    module procedure gatherv_real
    module procedure gatherv_real2mult
    module procedure measure_gatherv
    module procedure thermo_gatherv
end interface gatherv

interface scatterv 
    module procedure measure_scatterv
end interface scatterv 

contains 

!NOT FINISHED
subroutine get_two_level_comm(com_in,N,com_outer,com_inner)
    class(mpi_type),intent(in)      :: com_in
    integer                         :: N !number of entries that have to be distributed
    type(mpi_distv),intent(out)     :: com_outer
    type(mpi_type),intent(out)      :: com_inner

#ifdef CPP_MPI
    integer     :: div,color
    integer     :: mpi_comm_tmp !temporary MPI_comm
    integer     :: ierr,i

    if(N>=com_in%NP)then
        Call com_outer%init(com_in)
        com_outer%cnt=N/com_in%NP
        com_outer%cnt(1:modulo(N,com_in%NP))=com_outer%cnt(1:modulo(N,com_in%NP))+1
        com_outer%displ=0
        do i=2,com_in%NP
            com_outer%displ(i)=sum(com_outer%cnt(1:i-1))
        enddo
        Call MPI_COMM_SPLIT(com_in%com,com_in%id,com_in%id,mpi_comm_tmp,ierr)    !put everybody into own inner communicator to skip inner parallelization
        Call com_inner%set(mpi_comm_tmp)
    else
        div=(com_in%Np-1)/N+1
        color=com_in%id/div
        Call MPI_COMM_SPLIT(com_in%com,color,com_in%id,mpi_comm_tmp,ierr)
        Call com_inner%set(mpi_comm_tmp)
        color=1
        if(com_inner%ismas) color=0
        Call MPI_COMM_SPLIT(com_in%com,color,com_in%id,mpi_comm_tmp,ierr)
        Call com_outer%set(mpi_comm_tmp)
        allocate(com_outer%cnt(com_outer%Np),source=0)
        allocate(com_outer%displ(com_outer%Np),source=0)
        com_outer%cnt=[(1,i=1,com_outer%Np)]
        com_outer%displ=[(i-1,i=1,com_outer%Np)]
    endif
#else
    com_inner=com_in
    Call com_outer%init(com_in)
    com_outer%cnt=N
#endif


end subroutine



!subroutine set_distrib_com(com_in,N,com_out)
!    class(mpi_type),intent(in)      :: com_in
!    integer                         :: N !number of entries that have to be distributed
!    type(mpi_distv),intent(out)     :: com_out
!
!    integer     :: i
!    integer     :: id_task(N) 
!    integer     :: id_thread(com_in%Np) 
!    Call com_out%init(com_in)
!
!    if(N>=com_out%NP)then
!        do i=1,N
!            id_task(i)=modulo(i-1,com_out%NP)+1
!        enddo
!        do i=1,com_out%NP
!            com_out%cnt(i)=count(id_task==i)
!        enddo
!        do i=2,com_out%NP
!            com_out%displ(i)=sum(com_out%cnt(:i-1))
!        enddo
!    else
!    endif
!    CALL MPI_BARRIER(com_in%com,i)
!    STOP "DERPDERPDERP"
!
!end subroutine


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

subroutine gatherv_real(arr,com)
    use mpi_basic
    real(8),intent(in)             :: arr(:)
    class(mpi_distv),intent(in)    :: com
#ifdef CPP_MPI    
    integer                        :: ierr

    Call MPI_Gatherv(arr(1),com%cnt(com%id+1),MPI_DOUBLE_PRECISION,arr(1),com%cnt,com%displ,MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine gatherv_real2mult(arr,mult,com)
    use mpi_basic
    real(8),intent(in)             :: arr(:,:)
    class(mpi_distv),intent(in)    :: com
    integer,intent(in)             :: mult(com%Np)
#ifdef CPP_MPI    
    integer    :: ierr
    integer    :: cnt(size(com%cnt))
    integer    :: displ(size(com%displ)) 

    cnt=com%cnt*mult
    displ=com%displ*mult

    Call MPI_Gatherv(arr(1,1),cnt(com%id+1),MPI_DOUBLE_PRECISION,arr(1,1),cnt,displ,MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

end module

