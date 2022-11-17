module m_correlation_int
use m_correlation_base
use m_io_files_utils
use m_io_utils
use iso_fortran_env, only: int64,output_unit
implicit none

private
public :: corre_int,Correlation

type,extends(corre_base) :: corre_int

    contains

    procedure :: read_options

end type





interface Correlation
   module procedure Correlations
   module procedure autocorre_1d
   module procedure autocorre_2d
end interface Correlation

contains


subroutine read_options(this)
   class(corre_int), intent(inout)  :: this

   integer :: io_in,correlation_number
   character(len=100) :: correlation_descriptor

   io_in=open_file_read('input')
   call get_parameter(io_in,'input','correlation_number',correlation_number)
   call get_parameter(io_in,'input','correlation_descriptor',correlation_descriptor)
   call close_file('input',io_in)

   this%get_correlation => auto_correlation

   write(*,*) correlation_number,trim(correlation_descriptor)

   stop



end subroutine read_options



pure subroutine auto_correlation(data,t,tprime,chunk,r)
   real(8),intent(in)            :: data(:)
   integer,intent(in)            :: t,tprime,chunk(:)
   real(8),intent(out)           :: r(:)


   integer :: i,N_data
   real(8) :: ave,norm

   r=0.0d0

   N_data=size(data)

   ave=sum(data(1:tprime))/tprime
   norm=sum((data(1:tprime)-ave)**2)

   if (norm.lt.1.0d-8) return


   do i=1,tprime
      r=r+(data(i)-ave)*(data(i+t)-ave)
   enddo
   r=r/norm


end subroutine















! ===============================================================
! calculate the correlation length
! written by Lukas Deuchler
! date 11/4/2013
! email deuchlerl@gmail.com
      function Correlations(spin_sum,spin,shape_spin,n_MC,N_cell)
      use m_vector, only : norm
      Implicit none
      integer, intent(in) :: n_MC,shape_spin(:)
      real(kind=8), intent(in) :: spin(:,:,:,:,:),spin_sum(:,:,:,:,:),N_cell
! value of the function
      real(kind=8) :: Correlations(3)
! internal variable
      real(kind=8) :: a_cor(3),b_cor(3)
      real(kind=8) :: disti,E_int
      integer :: k,j1,j2,j3,j4,i1,i2,i3,i4

!!! use spin(1:6,N_site)
! Spin(1:3,N_site) are the positions on the lattice
! Spin(4:6,N_site) are the magnetic moments

! dummy variable
      Correlations=0.0d0
      disti=0

#ifdef CPP_OPENMP
       do k=1,3
       E_int=0.0d0
!$OMP parallel do reduction(+:E_int) private(i4,i3,i2,i1,j4,j3,j2,j1,a_cor,b_cor) default(shared)
#else
       do k=1,3
       E_int=0.0d0
#endif
        do i4=1,shape_spin(5)
         do i3=1,shape_spin(4)
          do i2=1,shape_spin(3)
           do i1=1,shape_spin(2)

          a_cor(k)=(spin(k,i1,i2,i3,i4) - spin_sum(k,i1,i2,i3,i4)/n_MC)

         Do j4=1,shape_spin(5)
          do j3=1,shape_spin(4)
           do j2=1,shape_spin(3)
            do j1=1,shape_spin(2)

          b_cor(k)=(spin(k,j1,j2,j3,j4) - spin_sum(k,j1,j2,j3,j4)/n_MC)

          E_int=E_int+a_cor(k)*b_cor(k)*sqrt(dble((j3-i3)**2+(j2-i2)**2+(j1-i1)**2))

           enddo
          enddo
         enddo
        enddo
           enddo
          enddo
         enddo
        enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
       Correlations(k)=E_int
       enddo

        Correlations=Correlations/N_cell**2

        end function Correlations
! ===============================================================

! calculate the autocorrelation of several multidimensional input
       function autocorre_1d(input,duration)
       implicit none
       real(kind=8), intent(in) :: input(:)
       integer, intent(in) :: duration
       real(kind=8) :: autocorre_1d(duration)
!internal
       integer :: i,j
       real(kind=8) :: norm,average

       autocorre_1d=0.0d0

       average=sum(input)/dble(duration)
       norm=sum((input-average)**2)
       if (norm.eq.0.0d0) return

       do i=1,duration
        do j=1,duration-i+1
         autocorre_1d(i)=autocorre_1d(i)+(input(j)-average)*(input(j+i-1)-average)
        enddo
         autocorre_1d(i)=autocorre_1d(i)/norm
       enddo

       end function autocorre_1d

       function autocorre_2d(input,N,duration)
       implicit none
       real(kind=8), intent(in) :: input(:,:)
       integer, intent(in) :: duration,N
       real(kind=8) :: autocorre_2d(N,duration)
!internal
       integer :: i,j,k
       real(kind=8) :: norm(N),average(N)

       autocorre_2d=0.0d0

       do i=1,N
        average(i)=sum(input(i,:))/dble(duration)
        norm(i)=sum((input(i,:)-average)**2)
        if (norm(i).eq.0.0d0) return
       enddo

       do i=1,duration
        do k=1,N
         do j=1,duration-i+1
          autocorre_2d(k,i)=autocorre_2d(k,i)+(input(k,j)-average(k))*(input(k,j+i-1)-average(k))
         enddo
         autocorre_2d(k,i)=autocorre_2d(k,i)/norm(k)
        enddo
       enddo

       end function autocorre_2d

end module m_correlation_int
