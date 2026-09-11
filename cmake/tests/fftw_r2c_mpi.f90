function my_function(i,j) result(res)
  integer(8), intent(in) :: i,j
  real(8) :: res

  res=exp(-1.0d0/(real(i,8)**2+real(j,8)**2)**2)

end function

program toto

  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3-mpi.f03'
  include 'mpif.h'
  integer(C_INTPTR_T), parameter :: L = 10000
  integer(C_INTPTR_T), parameter :: M = 10000
  type(C_PTR) :: plan, cdata, rdata
  complex(C_DOUBLE_COMPLEX), pointer :: data(:,:)
  real(C_DOUBLE), pointer :: in(:,:)
  integer(C_INTPTR_T) :: i, j, alloc_local, local_M, local_j_offset
  real(8) :: my_function
  integer process_Rank, size_Of_Cluster, ierror

  Call mpi_init(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

  call fftw_mpi_init()

!   get local data size and allocate (note dimension reversal)
  alloc_local = fftw_mpi_local_size_2d(M, L/2+1, MPI_COMM_WORLD, &
                                       local_M, local_j_offset)
  cdata = fftw_alloc_complex(alloc_local)
  call c_f_pointer(cdata, data, [L/2,local_M])

! split the real space in data
    rdata = fftw_alloc_real(2*alloc_local)
    call c_f_pointer(rdata, in, [2*(L/2+1),local_M])

  write(*,*) process_Rank, 'local offset', local_M, local_j_offset, 'Before plan'
!   create MPI plan for in-place forward DFT (note dimension reversal)
  plan = fftw_mpi_plan_dft_r2c_2d(M, L, in, data, MPI_COMM_WORLD, &
                              FFTW_MEASURE)

! initialize data to some function my_function(i,j)
  in=0.0d0
  do j = 1, local_M
    do i = 1, L
      in(i, j) = my_function(i, j + local_j_offset)
    end do
  end do

  write(*,*) process_Rank, 'before mpi_execute'
! compute transform (as many times as desired)
  call fftw_mpi_execute_dft_r2c(plan, in, data)
  write(*,*) process_Rank, 'after mpi_execute'

  call fftw_destroy_plan(plan)
  call fftw_free(cdata)

  call MPI_FINALIZE(ierror)

end program
