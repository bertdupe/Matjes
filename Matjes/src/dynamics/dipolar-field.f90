module m_dipolar_field
use m_derived_types, only : vec_point

!
! contains all the routines to calculate the dipolar field
!
!

! matrix that contains all the distances between the spins
real(kind=8), public, protected, allocatable, dimension(:,:,:) :: distances
! pointer that points to the spins of interests
type(vec_point),public, protected, allocatable, dimension(:) :: mode_dipole
! logical to turn on or off the dipole
logical, public, protected :: i_dip
! magnetic constant in H/nm
real(kind=8), Parameter :: mu_0=1.256637d-15
! conversion parameter (I do not remember what it is)
real(kind=8), parameter :: alpha=6.74582d-7
! saturation magnetization
real(kind=8), public, protected :: Ms

private
public :: get_dipole_B,get_ham_dipole

contains

subroutine get_dipole_B(B,iomp)
use m_constants, only : pi
use m_vector, only : norm
use m_constants, only : mu_B
implicit none
integer, intent(in) :: iomp
real(kind=8), intent(inout) :: B(:)
! internal
integer :: i,N
real(kind=8) :: rc(3),ss,B0(3),B_int(3)

N=size(mode_dipole)
B0=1.0d0/3.0d0
B_int=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel do reduction(+:B) private(i,j,rc,ss) default(shared)
#endif

do i=1,N

   if (iomp.eq.i) then
     B_int=B_int+B0
   else
     rc=distances(:,i,iomp)
     ss=norm(rc)
     B_int=B_int+(mode_dipole(i)%w*ss**2-3.0d0*dot_product(rc,mode_dipole(i)%w)*rc)/ss**5
   endif

enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

B(1:3)=B(1:3)+B_int/pi(2.0d0)*mu_B**2*Ms**2

end subroutine get_dipole_B


!
! check if we should calculate the dipole dipole interaction
!
subroutine get_ham_dipole(fname,my_lattice,motif)
use m_io_files_utils
use m_io_utils
use m_derived_types, only : cell,lattice
use m_get_position
use m_vector, only : norm
use m_operator_pointer_utils
use m_constants, only : mu_B
implicit none
type(cell), intent(in) :: motif
type(lattice), intent(in) :: my_lattice
character(len=*), intent(in) :: fname
! internal
integer :: i,j,io,N,k,l
real(kind=8), allocatable, dimension(:,:) :: positions
type(vec_point), allocatable, dimension(:) :: all_mode
real(kind=8) :: test(3),shorter_distance(3)
logical :: i_test

N=product(my_lattice%dim_lat)

io=open_file_read(fname)
call get_parameter(io,fname,'dipdip',i_dip)
call close_file(fname,io)

if (.not.i_dip) return

! get Ms from the motif in CGS units
!V=my_lattice%areal(:,1),
Ms=motif%atomic(1)%moment



allocate(distances(3,N,N))
distances=0.0d0
allocate(positions(3,N))
positions=0.0d0

io=open_file_read('positions.dat')
call get_position(positions,'positions.dat')
call close_file('positions.dat',io)

do i=1,N
   do j=1,N

      test=positions(:,j)-positions(:,i)

      do l=1,3
         if (my_lattice%boundary(l)) then
            do k=-1,1,2
               shorter_distance=test+real(k)*my_lattice%areal(:,l)*real(my_lattice%dim_lat(l))
               if (norm(shorter_distance).lt.norm(test)) test=shorter_distance
            enddo
         endif
      enddo

      distances(:,j,i)=positions(:,j)-positions(:,i)

   enddo
enddo

! make the pointer point to the spin matrix
allocate(mode_dipole(N),all_mode(N))
call dissociate(mode_dipole,N)
call dissociate(all_mode,N)
!! associate all_mode with all the lattice
call associate_pointer(all_mode,my_lattice)

call associate_pointer(mode_dipole,all_mode,'magnetic',i_test)

call dissociate(all_mode,N)

end subroutine get_ham_dipole

end module m_dipolar_field
