module m_random_init
    Implicit none
    private
    public :: random_init

contains

subroutine random_init()
    use mpi_basic, only: mpi_world
    use m_get_random, only: init_mtprng

#ifdef CPP_MRG
    call init_mtprng(mpi_world%id+1)
#else
    Call init_rand_seed(mpi_world%id+1)
#endif
end subroutine

subroutine init_rand_seed(id)
    !initialize internal random number generatal
    integer,intent(in)  :: id
    INTEGER :: i, n, clock
    INTEGER, DIMENSION(:), ALLOCATABLE :: seed
    
    CALL RANDOM_SEED(size = n)
    allocate(seed(n))
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed=clock+37*(id)*[ (i-1, i = 1, n) ]
    CALL RANDOM_SEED(PUT = seed)
    deallocate(seed)
end subroutine init_rand_seed

end module
