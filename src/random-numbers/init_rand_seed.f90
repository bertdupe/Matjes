      SUBROUTINE init_rand_seed()
#ifdef CPP_MPI
      use m_mpi_prop, only : irank
#endif
      Implicit none
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      
      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

#ifdef CPP_MPI
      seed=clock+37*(irank+1)*(/ (i - 1, i = 1, n) /)
#else
      seed=clock+37*(/ (i - 1, i = 1, n) /)
#endif
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
      END SUBROUTINE init_rand_seed
