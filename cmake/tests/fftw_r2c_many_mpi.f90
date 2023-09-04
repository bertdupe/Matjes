function my_function(i,j) result(res)
  integer(8), intent(in) :: i,j
  real(8) :: res

  res=exp(-50.0d0/(real(i,8)**2+real(j,8)**2)**2)

end function

program toto

  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3-mpi.f03'
  include 'mpif.h'
  integer(C_INTPTR_T), parameter :: L = 10
  integer(C_INTPTR_T), parameter :: M = 10
  integer(C_INTPTR_T), parameter :: N = 9
  type(C_PTR) :: plan, cdata, rdata
  complex(C_DOUBLE_COMPLEX), pointer :: data(:,:,:)
  real(C_DOUBLE), pointer :: in(:,:,:)
  integer(C_INTPTR_T) :: k, i, j, alloc_local, local_M, M_offset
  real(8) :: my_function
  integer process_Rank, size_Of_Cluster, ierror

  Call mpi_init(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

  call fftw_mpi_init()

!   get local data size and allocate (note dimension reversal)
   alloc_local = fftw_mpi_local_size_many(2,[M, L],N,FFTW_MPI_DEFAULT_BLOCK,MPI_COMM_WORLD,local_M,M_offset)

   cdata = fftw_alloc_complex(N*alloc_local)
   call c_f_pointer(cdata, data, [L,M,N])

! split the real space in data
    rdata = fftw_alloc_real(2*N*alloc_local)
    call c_f_pointer(rdata, in, [L,M,N])


  !   create MPI plan for in-place forward DFT (note dimension reversal)
  plan = fftw_mpi_plan_many_dft_r2c(3,[N,M,L],N, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, in, data,&
           &  MPI_COMM_WORLD,FFTW_ESTIMATE)

! initialize data to some function my_function(i,j)
  in=0.0d0
    do j = 1, local_M
      do i = 1, L
        in(j+M_offset,i,:) = my_function(i, j+M_offset)
      end do
    end do




! compute transform (as many times as desired)
  call fftw_mpi_execute_dft_r2c(plan, in, data)

  call fftw_destroy_plan(plan)
  call fftw_free(cdata)

  call MPI_FINALIZE(ierror)

end program
