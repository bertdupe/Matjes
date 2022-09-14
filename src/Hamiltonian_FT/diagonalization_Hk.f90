module m_diagonalization_Hk
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
use m_parameters_rw_high
use m_set_Hamiltonian_FT
use m_FT_Ham_public
use m_FT_Ham_base
use m_io_files_utils
use m_io_utils, only : get_parameter
use m_parameters_FT_Ham
! the following module is used for the TB. You will find the module in the directory tight-binding
use m_highsym, only : set_highs_path,mv_kpts
use m_kgrid_FT
use m_input_H_types
use m_derived_types
use m_FT_Ham_coo_rtok_base
implicit none

!
! 1. make the FT of each of the Ham individually
!  - create one type per Ham
!  - collect the real space differences for each of them
! 2. combine all the Hamiltonian and real space distances and the into a big Hamiltonian
! 3. for each of the k, get the phase factor and put it into the Hamiltonian
! 4. diagonalise Hk
!

private
public  :: diagonalize_Ham_FT

contains

subroutine diagonalize_Ham_FT(H_io,lat)
    type(io_h),intent(in)               :: H_io
    type(lattice), intent(in)           :: lat

    ! internal Hamiltonians
    type(H_inp_real_to_k),allocatable :: FT_Ham(:)
    class(FT_Ham_base),allocatable    :: Hk
    type(parameters_FT_HAM_IO)        :: io_H_diag         ! parameters for the diagonalization
    ! high symmetry lines
    type(parameters_IO_HIGHS)         :: high_lines
    ! full BZ DOS type
    type(k_grid_FT)                   :: fullBZ_kmesh

    ! dummy variables
    real(8)   :: k(3)
    integer   :: io_input,io_evec
    integer   :: i,n_kpts,n_eigen,shape_vec(2),j,N(3)
    real(8), allocatable :: kpts(:,:)
    complex(8),allocatable  :: eigenvalues(:)      ! array containing the eigenvalues
    complex(8),allocatable  :: eigenvectors(:,:)   ! array containing the eigenvectors
    character(len=100)      :: form_eig,form_vec
    logical                 :: do_highs_k = .True.,do_evec = .False.,do_dos_k = .False.

    ! initialization
    k=0.0d0
    shape_vec=0
    io_evec=0
    io_input=open_file_read('input')
    call io_H_diag%read_file(io_input,'input')

    call get_parameter(io_input,'input','do_highs_k',do_highs_k)
    call get_parameter(io_input,'input','do_dos_k',do_dos_k)
    if ((do_dos_k).and.(do_highs_k)) STOP 'set do_highs_k or do_dos_k to .False.'

    call get_parameter(io_input,'input','calc_evec',do_evec)

    ! read the high symmetry lines
    if (do_highs_k) then
       call high_lines%read_file('q',io_input,'input')
       call set_highs_path(lat,high_lines)
       call mv_kpts(kpts)
    endif

    ! get the full BZ kpoint
    if (do_dos_k) then
       call fullBZ_kmesh%kmesh_read('q',io_input,'input')
       N=fullBZ_kmesh%get_kgrid()
       call fullBZ_kmesh%set(lat%astar,N)
       n_kpts=fullBZ_kmesh%get_Nk()
       allocate(kpts(3,n_kpts),source=0.0d0)

       do i=1,n_kpts
          kpts(:,i)=fullBZ_kmesh%get_k(i)
       enddo
    endif

    call close_file('input',io_input)


    ! prepare the Hamiltonian based on the coo matrices for the FT
    call set_Hamiltonians_FT(FT_Ham,H_io,lat)     ! choose with which algoritm you want to work

    Call set_H(Hk,io_H_diag)   ! choose the Hamiltonian with which you would like to work (sparse, dense...)

    call Hk%init(FT_Ham,io_H_diag)    ! initialize the Hamiltonian matrix
    call Hk%set_work(eigenvalues,eigenvectors)

    io_input=open_file_write('dispersion.dat')
    write(form_eig,'( "(3(E20.12E3,x),", I10, "(x,E20.12E3,x,E20.12E3))" )') size(eigenvalues)
    if (do_evec) then
    ! write the eigen vectors
         shape_vec=shape(eigenvectors)
         io_evec=open_file_write('eigenvectors.dat')
         write(output_unit,'(/a)') "search for eigenvectors"
         write(output_unit,'(a/)') "list of eigenvectors - format kpts, (Real,Cmplx):"
         write(form_vec,'( "(3(E20.12E3,x),5x,",I10,"(E20.12E3,x))")') 2*shape_vec(1)*shape_vec(2)
    endif

    n_kpts=size(kpts,2)
    do i=1,n_kpts
       call Hk%set_k(FT_Ham,kpts(:,i))
       if (do_evec) then
           call Hk%calc_evec(size(eigenvalues),eigenvalues,eigenvectors,n_eigen)
           write(io_evec,form_vec) (eigenvectors(:,j),j=1,n_eigen)
        else
           call Hk%calc_eval(size(eigenvalues),eigenvalues,n_eigen)
       endif
       write(io_input,form_eig) kpts(:,i),eigenvalues
    enddo

    call close_file('dispersion.dat',io_input)

    if (do_evec) call close_file('eigenvectors.dat',io_evec)

    stop 'dispersion done'
end subroutine

end module m_diagonalization_Hk
