module m_dipolar_field
use m_derived_types, only : vec_point
use m_fftw

!
! contains all the routines to calculate the dipolar field
!

! matrix that contains the FFT of the modes
complex(kind=8), allocatable :: FFT_mode(:,:)

! matrix that contains all the distances between the spins
real(kind=8), public, protected, allocatable, dimension(:,:) :: distances
! pointer that points to the spins of interests
type(vec_point),public, protected, allocatable, dimension(:) :: mode_dipole
! logical to turn on or off the dipole
logical, public, protected :: i_dip



!
! In case of the Ewald sum
!
integer :: N_shell_Ewald=10


!
!
! The Tesla are in H.A.m^-2
! magnetic constant in H/m
!real(kind=8), Parameter :: mu_0=12.5663706144d-7
! bohr magneton in A.m^2
!real(kind=8), Parameter :: mu_B=9.2740094980d-24
! to avoid going in too small numbers, we express the bohr magneton in A.m^2 * 10^27
!real(kind=8), Parameter :: mu_B=9.2740094980d3
! this constant contains 9.2740094980d-24*10^27*12.5663706144d-7=mu_b*mu_0*10^27
! this constant is in Tesla
real(kind=8), Parameter :: alpha=0.0116541d0
! saturation magnetization
real(kind=8), public, protected :: Ms

private
public :: get_dipole_B,get_ham_dipole,prepare_FFT_dipole,calculate_FFT_modes

contains

!
! if you want to calculate the Dipole Dipole with the internal FFTW libraries and the Ewald sum
!
#ifdef CPP_SUM
subroutine get_dipole_B(B,iomp)
use m_constants, only : pi,mu_b
use m_vector, only : norm
implicit none
integer, intent(in) :: iomp
real(kind=8), intent(inout) :: B(:)
! internal
integer :: i,N
real(kind=8) :: rc(3),ss,B0,B_int(3),r0(3)

N=size(mode_dipole)
B0=2.0d0/3.0d0
B_int=0.0d0
r0=distances(:,iomp)

#ifdef CPP_OPENMP
!$OMP parallel do reduction(+:B_int) private(i,rc,ss) default(shared)
#endif

do i=1,N

   if (iomp.ne.i) then
     rc=distances(:,i)-r0
     ss=norm(rc)
     B_int=B_int+(mode_dipole(i)%w*ss**2-3.0d0*dot_product(rc,mode_dipole(i)%w)*rc)/ss**5
   endif

enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

!
! the dipole field is in T and should be converted to eV. So it is multiplied by mu_B in eV/T and mu_0 which is one in the code
!

B(1:3)=B(1:3)-B_int/pi(4.0d0)*alpha*Ms**2*mu_b

end subroutine get_dipole_B
#endif











!
! just check if the FFT dipole should be calculated
!
subroutine prepare_FFT_dipole(N)
implicit none
integer, intent(in) :: N
! internal
logical :: check_dipole
integer :: N_k,dim_mode

check_dipole=.false.
#ifdef CPP_BRUTDIP
check_dipole=i_dip
#endif
#ifdef CPP_FFTW
check_dipole=i_dip
#endif
if (.not.check_dipole) return

dim_mode=size(mode_dipole(1)%w)
N_k=product(N_kpoint)
allocate(FFT_mode(dim_mode,N_k))
FFT_mode=0.0d0

end subroutine

!
! just check if the FFT dipole should be calculated
!
subroutine calculate_FFT_modes(iteration)
implicit none
integer, intent(in) :: iteration
! internal
logical :: check_dipole
integer :: dim_mode

check_dipole=.false.
#ifdef CPP_BRUTDIP
check_dipole=i_dip
#endif
#ifdef CPP_FFTW
check_dipole=i_dip
#endif

if (.not.check_dipole) return

write(6,'(a,2x,I10)') 'calculation of the FFT of the modes for the interation',iteration
!
! calculate the FFT of the modes. FFT comes out
!
dim_mode=size(mode_dipole(1)%w)
call calculate_fft(mode_dipole,distances,-1.0d0,dim_mode,FFT_mode)

write(6,'(a)') 'done'

end subroutine





!
! if you want to calculate the Dipole Dipole with the internal FFT
!
#ifdef CPP_BRUTDIP
subroutine get_dipole_B(B,iomp)
use m_constants, only : pi,mu_b
use m_vector, only : norm
implicit none
integer, intent(in) :: iomp
real(kind=8), intent(inout) :: B(:)
! internal
integer :: i,N_k
real(kind=8) :: phase
complex(kind=8) :: D_int(3,3),B_int(3)

N_k=product(N_kpoint)

B_int=0.0d0
#ifdef CPP_OPENMP
!$OMP parallel do reduction(+:B_int) private(D_int,phase) default(shared)
#endif

do i=1,N_k

    phase=dot_product(distances(:,iomp),kmesh(:,i))
    D_int=FFT_pos_D(:,:,i)*complex(cos(phase),sin(phase))
    B_int=B_int+matmul(FFT_mode(:,i),D_int)

enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

!
! the dipole field is in T and should be converted to eV. So it is multiplied by mu_B in eV/T and mu_0 which is one in the code
!
do i=1,3
   B(i)=B(i)-real(B_int(i))/pi(4.0d0)*alpha*Ms**2*mu_b
enddo

end subroutine get_dipole_B

#endif









!
! if you want to calculate the Dipole Dipole with the double sum
!
#ifdef CPP_FFTW
subroutine get_dipole_B(B,iomp)
use m_constants, only : pi,mu_b
use m_vector, only : norm
implicit none
integer, intent(in) :: iomp
real(kind=8), intent(inout) :: B(:)
! internal
integer :: i,N
real(kind=8) :: rc(3),ss,B0,B_int(3)

N=size(mode_dipole)
B0=2.0d0/3.0d0
B_int=0.0d0

#ifdef CPP_OPENMP
!$OMP parallel do reduction(+:B_int) private(i,rc,ss) default(shared)
#endif

do i=1,N

   if (iomp.ne.i) then
     rc=distances(:,i,iomp)
     ss=norm(rc)
     B_int=B_int+(mode_dipole(i)%w*ss**2-3.0d0*dot_product(rc,mode_dipole(i)%w)*rc)/ss**5
   endif

enddo
#ifdef CPP_OPENMP
!$OMP end parallel do
#endif

!
! the dipole field is in T and should be converted to eV. So it is multiplied by mu_B in eV/T and mu_0 which is one in the code
!

B(1:3)=B(1:3)-B_int/pi(4.0d0)*alpha*Ms**2*mu_b

end subroutine get_dipole_B
#endif

!
! check if we should calculate the dipole dipole interaction
!
subroutine get_ham_dipole(fname,my_lattice,motif)
use m_io_files_utils
use m_io_utils
use m_derived_types, only : cell,lattice
use m_get_position
use m_operator_pointer_utils
use m_constants, only : mu_B
implicit none
type(cell), intent(in) :: motif
type(lattice), intent(in) :: my_lattice
character(len=*), intent(in) :: fname
! internal
integer :: i,io,N
real(kind=8), allocatable, dimension(:,:) :: positions
type(vec_point), allocatable, dimension(:) :: all_mode
real(kind=8) :: r(3,3)
logical :: i_test

N=product(my_lattice%dim_lat)

io=open_file_read(fname)
call get_parameter(io,fname,'dipdip',i_dip)
call close_file(fname,io)

if (.not.i_dip) return

! get Ms from the motif in CGS units
!V=my_lattice%areal(:,1),
Ms=motif%atomic(1)%moment

! put the lattice in the correct order for Fortran
do i=1,3
   r(:,i)=my_lattice%areal(i,:)
enddo

allocate(distances(3,N))
distances=0.0d0
allocate(positions(3,N))
positions=0.0d0

call get_position(positions,'positions.dat')

!
! prepare the dipolar matrix (the matrix of the r to calculate the 1/r)
!

call calculate_distances(distances,positions,r,my_lattice%dim_lat,my_lattice%boundary)

!
!
#ifdef CPP_BRUTDIP
call set_k_mesh('input',my_lattice)
call calculate_FFT_Dr(distances,-1.0d0,my_lattice%dim_lat)

#endif
!
!

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
