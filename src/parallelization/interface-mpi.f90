      module m_mpi
      interface allreduce
       module procedure allreduce_scalar
       module procedure allreduce_1D
       module procedure allreduce_2D
       module procedure allreduce_4D
       module procedure allreduce_3D
       module procedure allreduce_3D_bis
      end interface allreduce

      interface reduce
       module procedure reduce_vector_5d
       module procedure reduce_vector_2D
       module procedure reduce_vector_1D
       module procedure reduce_vector_I_1D
       module procedure reduce_scalar_r
       module procedure reduce_scalar_i
      end interface reduce

      interface bcast
       module procedure bcast_5d
      end interface bcast

      interface gather
       module procedure gather_1d
       module procedure gather_I_1d
       module procedure gather_2d
      end interface gather

      interface scatter
       module procedure scatter_1d
       module procedure scatter_I_1d
       module procedure scatter_2d
      end interface scatter

      contains
!----------------------------------------------------------------------------------------------------------
      function gather_1d(v,length,size_table,root,MPI_COMM)
      Implicit none
      integer, intent(in) :: length,root,MPI_COMM,size_table
      real(kind=8),dimension(size_table) :: gather_1d
      real(kind=8), intent(in) :: v(length)
      integer :: ierr

      include 'mpif.h'

      gather_1d=0.0d0

      call mpi_gather(v,length,MPI_REAL8,gather_1d,length,MPI_REAL8,root,MPI_COMM,ierr)

      end function gather_1d
!----------------------------------------------------------------------------------------------------------
      function gather_I_1d(v,length,isize,root,MPI_COMM)
      Implicit none
      integer, intent(in) :: length,root,MPI_COMM,isize
      integer,dimension(length*isize) :: gather_I_1d
      integer, intent(in) :: v(length)
      integer :: ierr

      include 'mpif.h'

      gather_I_1d=0

      call mpi_gather(v,length,MPI_INT,gather_I_1d,length,MPI_INT,root,MPI_COMM,ierr)

      end function gather_I_1d
!----------------------------------------------------------------------------------------------------------
      function gather_2d(v,width,length,size_table,root,MPI_COMM)
      Implicit none
      integer, intent(in) :: width,length,root,MPI_COMM,size_table
      real(kind=8),dimension(width,size_table) :: gather_2d
      real(kind=8), intent(in) :: v(width,length)
      integer :: ierr,n_comm

      include 'mpif.h'

      gather_2d=0.0d0
      n_comm=width*length

      call mpi_gather(v,n_comm,MPI_REAL8,gather_2d,n_comm,MPI_REAL8,root,MPI_COMM,ierr)

      end function gather_2d
!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------

      function scatter_1d(v,length,isize,root,MPI_COMM)
      Implicit none
      integer, intent(in) :: length,root,MPI_COMM,isize
      real(kind=8),dimension(length) :: scatter_1d
      real(kind=8), intent(in) :: v(isize*length)
      integer :: ierr

      include 'mpif.h'

      scatter_1d=0.0d0

      call mpi_scatter(v,length,MPI_REAL8,scatter_1d,length,MPI_REAL8,root,MPI_COMM,ierr)

      end function scatter_1d
!----------------------------------------------------------------------------------------------------------
      function scatter_I_1d(v,length,isize,root,MPI_COMM)
      Implicit none
      integer, intent(in) :: length,root,MPI_COMM,isize
      integer,dimension(length) :: scatter_I_1d
      integer, intent(in) :: v(isize*length)
      integer :: ierr

      include 'mpif.h'

      scatter_I_1d=0

      call mpi_scatter(v,length,MPI_INT,scatter_I_1d,length,MPI_INT,root,MPI_COMM,ierr)

      end function scatter_I_1d
!----------------------------------------------------------------------------------------------------------
      function scatter_2d(v,width,length,isize,root,MPI_COMM)
      Implicit none
      integer, intent(in) :: width,length,root,MPI_COMM,isize
      real(kind=8),dimension(width,length) :: scatter_2d
      real(kind=8), intent(in) :: v(width,isize*length)
      integer :: ierr,n_comm

      include 'mpif.h'

      scatter_2d=0.0d0
      n_comm=width*length

      call mpi_scatter(v,n_comm,MPI_REAL8,scatter_2d,n_comm,MPI_REAL8,root,MPI_COMM,ierr)

      end function scatter_2d

!----------------------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------
      function bcast_5d(v,n1,n2,n3,n4,n5,mpi_type,root,MPI_COMM)
      Implicit none
      real(kind=8),dimension(n1,n2,n3,n4,n5) :: bcast_5d
      integer, intent(in) :: n1,n2,n3,n4,n5,mpi_type,root,MPI_COMM
      real(kind=8), intent(in) :: v(n1,n2,n3,n4,n5)
      integer :: ierr,n_comm

      include 'mpif.h'

      n_comm=n1*n2*n3*n4*n5

      call MPI_BCAST(v,n_comm,mpi_type,root,MPI_COMM,ierr)

      bcast_5d=v

      end function bcast_5d

      function allreduce_2D(v,n,MPI_COMM)
      Implicit none
  ! intent(in)
      integer, intent(in) :: n(2),MPI_COMM
      real(kind=8), intent(in) :: v(:,:)
! value of the function
      real(kind=8),dimension(n(1),n(2)) :: allreduce_2D
      integer :: N_comm,ierr

      include 'mpif.h'

      N_comm=product(n)
      allreduce_2D=0.0d0

      call mpi_allreduce(v,allreduce_2D,N_comm,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)

      end function allreduce_2D

      function allreduce_3D_bis(v,n,p,MPI_COMM)
      Implicit none
  ! intent(in)
      integer, intent(in) :: n(2),p,MPI_COMM
      real(kind=8), intent(in) :: v(:,:,:)
! value of the function
      real(kind=8),dimension(p,n(1),n(2)) :: allreduce_3D_bis
      integer :: N_comm,ierr(3)

      include 'mpif.h'

      N_comm=product(n)*p
      allreduce_3D_bis=0.0d0

      call mpi_allreduce(v,allreduce_3D_bis,N_comm,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)

      end function allreduce_3D_bis

      function allreduce_3D(v,n,MPI_COMM)
      Implicit none
  ! intent(in)
      integer, intent(in) :: n(3),MPI_COMM
      real(kind=8), intent(in) :: v(:,:,:)
! value of the function
      real(kind=8),dimension(n(1),n(2),n(3)) :: allreduce_3D
      integer :: N_comm,ierr

      include 'mpif.h'

      N_comm=product(n)
      allreduce_3D=0.0d0

      call mpi_allreduce(v,allreduce_3D,N_comm,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)

      end function allreduce_3D

      function allreduce_4D(v,n,length,MPI_COMM)
      Implicit none
  ! intent(in)
      integer, intent(in) :: n(3),length,MPI_COMM
      real(kind=8), intent(in) :: v(:,:,:,:)
! value of the function
      real(kind=8),dimension(length,n(1),n(2),n(3)) :: allreduce_4D
      integer :: N_comm,ierr

      include 'mpif.h'

      N_comm=product(n)*length
      allreduce_4D=0.0d0

      call mpi_allreduce(v,allreduce_4D,N_comm,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)

      end function allreduce_4D

      function allreduce_scalar(v,MPI_COMM)
      Implicit none
  ! intent(in)
      integer, intent(in) :: MPI_COMM
      real(kind=8), intent(in) :: v
! value of the function
      real(kind=8) :: allreduce_scalar
      integer :: ierr

      include 'mpif.h'

      allreduce_scalar=0.0d0

      call mpi_allreduce(v,allreduce_scalar,1,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)

      end function allreduce_scalar

      function allreduce_1D(v,length,MPI_COMM)
      Implicit none
  ! intent(in)
      integer, intent(in) :: length,MPI_COMM
      real(kind=8), intent(in) :: v(length)
! value of the function
      real(kind=8), dimension(length) :: allreduce_1D
      integer :: ierr

      include 'mpif.h'

      allreduce_1D=0.0d0

      call mpi_allreduce(v,allreduce_1D,length,MPI_REAL8,MPI_SUM,MPI_COMM,ierr)

      end function allreduce_1D

!-------------- reduce part
!!!!!!!!!!!!!!!!!!!!!!!!!!!
      function reduce_vector_5d(v,N1,N2,N3,N4,N5,MPI_COMM)
      Implicit none
! intent(in)
      integer, intent(in) :: N1,N2,N3,N4,N5,MPI_COMM
      real(kind=8), intent(in) :: v(:,:,:,:,:)
! value of the function
      real(kind=8), dimension(N1,N2,N3,N4,N5) :: reduce_vector_5d
      integer :: ierr
! dummy

      include 'mpif.h'

      reduce_vector_5d=0.0d0

      call mpi_reduce(v,reduce_vector_5d,N1*N2*N3*N4*N5,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)

      end function reduce_vector_5d
!-----------------------
      function reduce_vector_2D(v,width,length,MPI_COMM)
      Implicit none
! intent(in)
      integer, intent(in) :: width,length,MPI_COMM
      real(kind=8), intent(in) :: v(:,:)
! value of the function
      real(kind=8), dimension(width,length) :: reduce_vector_2D
      integer :: ierr
! dummy

      include 'mpif.h'

      reduce_vector_2D=0.0d0

      call mpi_reduce(v,reduce_vector_2D,width*length,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)

      end function reduce_vector_2D
!-----------------------
      function reduce_vector_1D(v,length,MPI_COMM)
      Implicit none
! intent(in)
      real(kind=8), intent(in) :: v(:)
      integer, intent(in) :: length,MPI_COMM
! value of the function
      real(kind=8), dimension(length) :: reduce_vector_1D
! dummy
      integer :: ierr

      include 'mpif.h'

      reduce_vector_1D=0.0d0

      call mpi_reduce(v,reduce_vector_1D,length,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)

      end function reduce_vector_1D

!-----------------------
      function reduce_vector_I_1D(v,length,MPI_COMM)
      Implicit none
! intent(in)
      integer, intent(in) :: v(:)
      integer, intent(in) :: length,MPI_COMM
! value of the function
      integer, dimension(length) :: reduce_vector_I_1D
! dummy
      integer :: ierr

      include 'mpif.h'

      reduce_vector_I_1D=0

      call mpi_reduce(v,reduce_vector_I_1D,length,MPI_INT,MPI_SUM,0,MPI_COMM,ierr)

      end function reduce_vector_I_1D

!-----------------------
      function reduce_scalar_r(v,MPI_COMM)
      Implicit none
! intent(in)
      real(kind=8), intent(in) :: v
      integer, intent(in) :: MPI_COMM
! value of the function
      real(kind=8) :: reduce_scalar_r
! dummy
      integer :: ierr

      include 'mpif.h'

      reduce_scalar_r=0.0d0

      call mpi_reduce(v,reduce_scalar_r,1,MPI_REAL8,MPI_SUM,0,MPI_COMM,ierr)

      end function reduce_scalar_r

!-----------------------
      function reduce_scalar_i(v,isize,MPI_COMM)
      Implicit none
! intent(in)
     integer, intent(in) :: v
      integer, intent(in) :: isize,MPI_COMM
! value of the function
      integer :: reduce_scalar_i
! dummy
      integer :: ierr

      include 'mpif.h'

      reduce_scalar_i=0

      call mpi_reduce(v,reduce_scalar_i,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM,ierr)

      end function reduce_scalar_i
      end module m_mpi
