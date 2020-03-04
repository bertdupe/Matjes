module m_dipole_energy
use m_derived_types, only : vec_point
!
! contains all the routines to calculate the dipole dipole energy interaction
!
!

! matrix that contains all the distances between the spins
real(kind=8), protected, allocatable, dimension(:,:,:) :: distances
! pointer that points to the spins of interests
type(vec_point), protected, allocatable, dimension(:) :: mode_dipole
! logical to turn on or off the dipole
logical, public, protected :: i_dip
! conversion parameter (I do not remember what it is)
real(kind=8), parameter :: alpha=6.74582d-7
! saturation magnetization
real(kind=8) :: Ms

private
public :: get_dipole,get_ham_dipole

contains

!
! calculate the double sum of the dipole dipole energy
!
real(kind=8) function get_dipole(iomp)
use m_constants, only : pi
use m_vector, only : norm
use m_basic_types, only : vec
use m_constants, only : mu_B
implicit none
! input
integer, intent(in) :: iomp
! internal variable
integer :: i,j,N
! internal variable
real(kind=8) :: rc(3),ss

get_dipole=0.0d0
N=size(mode_dipole)

! the choosen spin is the spin i_s. It has nn nearest neighbours. The numbers of nearest
! neighbours are stored in n(:). for example n(1)=4 means 4 nearest neighbours
!!! first neighb

#ifdef CPP_OPENMP
!$OMP parallel do reduction(+:get_dipole) private(i,j,rc,ss) default(shared)
#endif

do i=1,N

   if (iomp.eq.i) then
     get_dipole=get_dipole-1.0d0/3.0d0
   else
     rc=distances(:,i,iomp)
     ss=norm(rc)
     get_dipole=get_dipole+(dot_product(mode_dipole(i)%w,mode_dipole(iomp)%w)*ss**2-3.0d0* &
          dot_product(rc,mode_dipole(i)%w)*dot_product(rc,mode_dipole(iomp)%w))/ss**5
   endif

enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif
get_dipole=get_dipole/pi(4.0d0)*mu_B**2*Ms**2

end function get_dipole

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

! get Ms from the motif
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






!
! get the dipole dipole position matrix elements
!
!function get_distances_dipole(r1,r2,periodic)
!implicit none
!integer, intent(in) :: N
!logical, intent(in) :: periodic(3)
!real(kind=8) :: get_distances_dipole(3)
!! internal
!
!write(*,*) 'get_positions_dipole'
!
!end function get_distances_dipole

!!!! part of the FFT dipole dipole interaction
#ifndef CPP_BRUTDIP
! prepare the demag tensor
!      if (i_dip) then
!       write(6,'(a)') "preparing spatial contribution to dipole convolution "
!       call setup_dipole(dim_lat,Periodic_log,net,count(my_motif%i_mom),size(world))
!       do k=1,dim_lat(3)*count(my_motif%i_mom)
!        do j=1,dim_lat(2)
!         do i=1,dim_lat(1)
!         mmatrix(1:3,i,j,k)=spin(4:6,i,j,k,mod(k-1,count(my_motif%i_mom))+1)*spin(7,i,j,k,1)
!         enddo
!        enddo
!       enddo
!       write(6,'(a)') 'FFT dipole dipole is set up'
!       else
!       allocate(ntensor(1,1,1,1),mmatrix(1,1,1,1),ctrans(1,1,1),rtrans(1,1,1), &
!     & hcomplex(1,1,1,1),mcomplex(1,1,1,1),hreal(1,1,1,1))
!
!       mmatrix=0.0d0
!       rtrans=0.0d0
!       hreal=0.0d0
!       ctrans=dcmplx(0.d0,0.d0)
!       mcomplex=dcmplx(0.d0,0.d0)
!       hcomplex=dcmplx(0.d0,0.d0)
!       ntensor=dcmplx(0.d0,0.d0)
!      endif
!#endif
#endif

end module m_dipole_energy
