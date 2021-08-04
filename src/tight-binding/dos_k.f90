module m_dos_k
use m_tb_types
use m_derived_types, only : lattice
use m_kgrid_int, only: k_mesh_int, kmesh_t, get_kmesh
use m_tb_k_public       !mode that contains the more efficient TB k-space type
use m_Hk
use, intrinsic :: iso_fortran_env, only : output_unit

use m_dos_util
use m_derived_types, only: k_grid_t
use m_constants, only : pi
use mpi_basic, only: mpi_type


private
public calc_dos_k

interface init_dos_bnd
    module procedure init_dos_bnd_sc
    module procedure init_dos_bnd_nc
end interface

interface init_dos_orb
    module procedure init_dos_orb_sc
    module procedure init_dos_orb_nc
end interface

interface init_dos_sig
    module procedure init_dos_sig_sc
    module procedure init_dos_sig_nc
end interface

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! MPI Hk implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_dos_k(Hk,io_H,lat,io_dos,is_sc,work,comm)
    !subroutine that calls the super-conducting or normal-conducting routines for the DOS in k-space
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: io_H
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work
    logical,intent(in)                      :: is_sc
    type(mpi_type),intent(in)               :: comm     !mpi communicator

    if(is_sc)then
        Call calc_dos_sc(Hk,io_H,lat,io_dos,work,comm)
    else
        Call calc_dos_nc(Hk,io_H,lat,io_dos,work,comm)
    endif
end subroutine

subroutine calc_dos_nc(Hk,h_io,lat,io_dos,work,comm)
    class(H_k_base),intent(inout)           :: Hk       !k-space Hamiltonian class
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work
    type(mpi_type),intent(in)               :: comm

    type(dos_nc)                            :: dos  !main dos type to gather normal dos without any projections

    integer                                 :: Nin,Nout, dimH !maximal number of eigenvalues, found number of eigenvalues, dimension of Hamiltonian
    real(8),allocatable                     :: eval(:)      !eigenvalue array (Nin)
    complex(8),allocatable                  :: evec(:,:)    !eigenvector array (Nin,dimH)
    logical                                 :: use_evec     !logical if eigenvectors are needed for dos-calculation

    !parameters to describe integration k-grid
    class(kmesh_t),allocatable              :: k_grid       !type which contains the k-grid 
    integer                                 :: ik,Nk        !k-point loop variable and total number
    integer                                 :: bnd(2)       !index boundaries for k-loop (considering mpi)
    integer                                 :: idos         !dos loop variable
    real(8)                                 :: k(3)         !local k-point

    !parameters for projection on a contiguous set of orbitals somewhere in the unit-cell 
    type(dos_bnd_nc),allocatable            :: dos_bnd(:)
    logical                                 :: use_bnd
    integer                                 :: Ndos_bnd

    !parameters for projection on an orbital considering all supercell periodic images
    type(dos_orb_nc),allocatable            :: dos_orb(:)
    logical                                 :: use_orb
    integer                                 :: Ndos_orb

    !normal dos at several smearing
    type(dos_nc),allocatable                :: dos_sig(:)
    integer                                 :: Ndos_sig

    type(dos_all_nc)                        :: dos_all

    !initialize used k-grid
    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh,comm)
    Nk=k_grid%get_NK()

    !initialize the various dos data-types
    Call dos%init(io_dos)                                                   !normal dos without projections
    if(io_dos%all_states) Call dos_all%init_mult(io_dos,Hk%dimH)            !projection on each state separately
    Call init_dos_bnd(io_dos,dos_bnd,use_bnd,Ndos_bnd)                      !projection on a contiguous set of orbitals somewhere in the unit-cell
    Call init_dos_orb(io_dos,dos_orb,use_orb,Ndos_orb,h_io%norb*h_io%nspin) !projection on an orbital considering all supercell periodic images 
    Call init_dos_sig(io_dos,dos_sig,Ndos_sig)                             !normal dos, but with several gauss smearings 

    !prepare data
    Nin=Hk%get_size_eval()  !get maximal number of eigenvalues
    allocate(eval(Nin),source=0.0d0)
    use_evec=use_bnd.or.use_orb.or.io_dos%all_states
    if(use_evec)then
        dimH=Hk%get_dimH()
        allocate(evec(dimH,Nin),source=(0.0d0,0.0d0))
    endif

    !calculate eigenvalues for each k and add to dos
    Call comm%get_loop_bnd(Nk,bnd) !get k-loop bounds for parallelization
    do ik=bnd(1),bnd(2)
        if(io_dos%print_kint.and.comm%ismas) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',bnd(2)

        !get the diagonalization at a chosen k-point
        k=k_grid%get_K(ik)
        if(use_evec)then
            Call Hk%get_evec(k,Nin,eval,evec,Nout,work) 
        else
            Call Hk%get_eval(k,Nin,eval,Nout,work) 
        endif

        !add calculated states to the various dos flavors
        Call dos%add(eval)
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_orb
            Call dos_orb(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_sig
            Call dos_sig(idos)%add(eval)
        enddo
        if(io_dos%all_states) Call dos_all%add(eval,evec)
    enddo

    !gather all dos-data on the master thread
    if(comm%Np>1) Call dos_mpi_sum(dos,dos_bnd,dos_orb,dos_sig,dos_all, comm) 

    !print resulting dos to files
    if(comm%ismas) Call print_dos_output(dos,dos_bnd,dos_orb,dos_sig,dos_all,"_nc") 
end subroutine

subroutine calc_dos_sc(Hk,h_io,lat,io_dos,work,comm)
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work
    type(mpi_type),intent(in)               :: comm

    type(dos_sc)                            :: dos  !main dos type to gather normal dos without any projections

    integer                                 :: Nin,Nout, dimH   !maximal number of eigenvalues, found number of eigenvalues, dimension of Hamiltonian
    real(8),allocatable                     :: eval(:)      !eigenvalue array (Nin)
    complex(8),allocatable                  :: evec(:,:)    !eigenvector array (Nin,dimH)

    !parameters to describe integration k-grid
    class(kmesh_t),allocatable              :: k_grid       !type which contains the k-grid 
    integer                                 :: ik,Nk
    integer                                 :: bnd(2)   !indice boundaries for k-loop
    real(8)                                 :: k(3)

    integer                                 :: idos

    !parameters for projection on a contiguous set of orbitals somewhere in the unit-cell 
    type(dos_bnd_sc),allocatable            :: dos_bnd(:)   
    integer                                 :: Ndos_bnd
    logical                                 :: use_bnd

    !parameters for projection on an orbital considering all supercell periodic images
    type(dos_orb_sc),allocatable            :: dos_orb(:)
    integer                                 :: Ndos_orb
    logical                                 :: use_orb

    !normal dos at several smearing
    type(dos_sc),allocatable                :: dos_sig(:)
    integer                                 :: Ndos_sig

    type(dos_all_sc)                        :: dos_all

    !initialize used k-grid
    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh,comm)
    Nk=k_grid%get_NK()

    !initialize the various dos data-types
    Call dos%init(io_dos)                                                   !normal dos without projections
    if(io_dos%all_states) Call dos_all%init_mult(io_dos,Hk%dimH/2)          !projection on each state separately
    Call init_dos_bnd(io_dos,dos_bnd,use_bnd,Ndos_bnd)                      !projection on a contiguous set of orbitals somewhere in the unit-cell
    Call init_dos_orb(io_dos,dos_orb,use_orb,Ndos_orb,h_io%norb*h_io%nspin) !projection on an orbital considering all supercell periodic images 
    Call init_dos_sig(io_dos,dos_sig,Ndos_sig)                              !normal dos, but with several gauss smearings 

    !prepare data
    Nin=Hk%get_size_eval() !get maximal number of eigenvalues
    dimH=Hk%get_dimH()
    allocate(eval(Nin),source=0.0d0)
    allocate(evec(dimH,Nin),source=(0.0d0,0.0d0))

    Call comm%get_loop_bnd(Nk,bnd) !get k-loop bounds for parallelization
    !calculate eigenvalues for each k and add to dos
    do ik=bnd(1),bnd(2)
        if(io_dos%print_kint.and.comm%ismas) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk

        !get the diagonalization at a chosen k-point
        k=k_grid%get_K(ik)
        Call Hk%get_evec(k,Nin,eval,evec,Nout,work)

        !add calculated states to the various dos flavors
        Call dos%add(eval,evec)
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_orb
            Call dos_orb(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_sig
            Call dos_sig(idos)%add(eval,evec)
        enddo
        if(io_dos%all_states) Call dos_all%add(eval,evec)
    enddo

    !gather all dos-data on the master thread
    if(comm%Np>1) Call dos_mpi_sum(dos,dos_bnd,dos_orb, dos_sig, dos_all, comm) 

    !print resulting dos to files
    if(comm%ismas) Call print_dos_output(dos,dos_bnd,dos_orb,dos_sig, dos_all,"_sc") 
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   MINOR HELP ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_dos_output(dos,dos_bnd,dos_orb,dos_sig,dos_all,specifier)
    class(dos_t),intent(inout)              :: dos
    class(dos_t),intent(inout)              :: dos_bnd(:)
    class(dos_t),intent(inout)              :: dos_orb(:)
    class(dos_t),intent(inout)              :: dos_sig(:)
    class(dos_mult),intent(inout)           :: dos_all
    character(len=*),intent(in)             :: specifier
    logical             :: alloc
    integer             :: i
    character(len=3)    :: dos_id
    character(len=10)   :: sigma_char
   
    !print normal dos without projections
    Call dos%print('dos_k'//specifier//'.dat')


    !print dos for projection on a contiguous set of orbitals somewhere in the unit-cell 
    do i=1,size(dos_bnd)
        write(dos_id,'(I0.3)') i
        Call dos_bnd(i)%print("dos_k"//specifier//"_bnd_"//dos_id//".dat")
    enddo

    !print dos for projection on an orbital considering all supercell periodic images
    do i=1,size(dos_orb)
        write(dos_id,'(I0.3)') i
        Call dos_orb(i)%print("dos_k"//specifier//"_bnd_"//dos_id//".dat")
    enddo

    !print dos for different gauss smearing
    do i=1,size(dos_sig)
        write(sigma_char,'(F010.6)') dos_sig(i)%get_sigma()
        Call dos_sig(i)%print("dos_k_smearing_"//trim(adjustl(sigma_char))//".dat")
    enddo


    if(dos_all%is_set()) Call dos_all%print("dos_k_all"//specifier//".dat")
end subroutine

subroutine dos_mpi_sum(dos,dos_bnd,dos_orb,dos_sig,dos_all, comm)
    class(dos_t),   intent(inout)           :: dos
    class(dos_t),   intent(inout)           :: dos_bnd(:)
    class(dos_t),   intent(inout)           :: dos_orb(:)
    class(dos_t),   intent(inout)           :: dos_sig(:)
    class(dos_mult),intent(inout)           :: dos_all
    type(mpi_type),intent(in)               :: comm
    integer     ::  i

    Call dos%mpi_sum(comm)
    do i=1,size(dos_bnd)
        Call dos_bnd(i)%mpi_sum(comm)
    enddo
    do i=1,size(dos_orb)
        Call dos_orb(i)%mpi_sum(comm)
    enddo
    do i=1,size(dos_sig)
        Call dos_sig(i)%mpi_sum(comm)
    enddo
    if(dos_all%is_set()) Call dos_all%mpi_sum(comm)
end subroutine

subroutine init_dos_bnd_nc(io_dos,dos_bnd,use_bnd,Ndos_bnd)
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(dos_bnd_nc),allocatable,intent(inout)  :: dos_bnd(:)
    logical,intent(inout)                       :: use_bnd
    integer,intent(inout)                       :: Ndos_bnd
    integer :: idos

    use_bnd=allocated(io_dos%bnd)
    if(use_bnd)then
        Ndos_bnd=size(io_dos%bnd,2)
        allocate(dos_bnd(Ndos_bnd))
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%init_bnd(io_dos,io_dos%bnd(:,idos))
        enddo
    else
        Ndos_bnd=0
        allocate(dos_bnd(0))
    endif
end subroutine

subroutine init_dos_bnd_sc(io_dos,dos_bnd,use_bnd,Ndos_bnd)
    type(parameters_TB_IO_DOS),intent(in)       :: io_dos
    type(dos_bnd_sc),allocatable,intent(inout)  :: dos_bnd(:)
    logical,intent(inout)                       :: use_bnd
    integer,intent(inout)                       :: Ndos_bnd
    integer :: idos

    use_bnd=allocated(io_dos%bnd)
    if(use_bnd)then
        Ndos_bnd=size(io_dos%bnd,2)
        allocate(dos_bnd(Ndos_bnd))
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%init_bnd(io_dos,io_dos%bnd(:,idos))
        enddo
    else
        Ndos_bnd=0
        allocate(dos_bnd(0))
    endif
end subroutine

subroutine init_dos_orb_nc(io_dos,dos_orb,use_orb,Ndos_orb,orbital_per_unitcell)
    type(parameters_TB_IO_DOS),intent(in)       :: io_dos
    type(dos_orb_nc),allocatable,intent(inout)  :: dos_orb(:)
    logical,intent(inout)                       :: use_orb
    integer,intent(inout)                       :: Ndos_orb
    integer,intent(in)                          :: orbital_per_unitcell
    integer :: idos

    use_orb=allocated(io_dos%orb)
    if(use_orb)then
        Ndos_orb=size(io_dos%orb)
        allocate(dos_orb(Ndos_orb))
        do idos=1,Ndos_orb
            Call dos_orb(idos)%init_bnd(io_dos,io_dos%orb(idos),orbital_per_unitcell)
        enddo
    else
        Ndos_orb=0
        allocate(dos_orb(0))
    endif
end subroutine

subroutine init_dos_orb_sc(io_dos,dos_orb,use_orb,Ndos_orb,orbital_per_unitcell)
    type(parameters_TB_IO_DOS),intent(in)       :: io_dos
    type(dos_orb_sc),allocatable,intent(inout)  :: dos_orb(:)
    logical,intent(inout)                       :: use_orb
    integer,intent(inout)                       :: Ndos_orb
    integer,intent(in)                          :: orbital_per_unitcell
    integer :: idos

    use_orb=allocated(io_dos%orb)
    if(use_orb)then
        Ndos_orb=size(io_dos%orb)
        allocate(dos_orb(Ndos_orb))
        do idos=1,Ndos_orb
            Call dos_orb(idos)%init_bnd(io_dos,io_dos%orb(idos),orbital_per_unitcell)
        enddo
    else
        Ndos_orb=0
        allocate(dos_orb(0))
    endif
end subroutine

subroutine init_dos_sig_nc(io_dos,dos_sig,Ndos_sig)
    type(parameters_TB_IO_DOS),intent(in)       :: io_dos
    type(dos_nc),allocatable,intent(inout)      :: dos_sig(:)
    integer,intent(out)                         :: Ndos_sig
    integer :: idos, N

    if(allocated(io_dos%sigma_arr))then
        Ndos_sig=size(io_dos%sigma_arr)
        allocate(dos_sig(Ndos_sig))
        do idos=1,Ndos_sig
            Call dos_sig(idos)%init(io_dos)
            Call dos_sig(idos)%set_sigma(io_dos%sigma_arr(idos))
        enddo
    else
        Ndos_sig=0
        allocate(dos_sig(0))
    endif
end subroutine

subroutine init_dos_sig_sc(io_dos,dos_sig,Ndos_sig)
    type(parameters_TB_IO_DOS),intent(in)       :: io_dos
    type(dos_sc),allocatable,intent(inout)      :: dos_sig(:)
    integer,intent(out)                         :: Ndos_sig
    integer :: idos, N

    if(allocated(io_dos%sigma_arr))then
        Ndos_sig=size(io_dos%sigma_arr)
        allocate(dos_sig(Ndos_sig))
        do idos=1,Ndos_sig
            Call dos_sig(idos)%init(io_dos)
            Call dos_sig(idos)%set_sigma(io_dos%sigma_arr(idos))
        enddo
    else
        Ndos_sig=0
        allocate(dos_sig(0))
    endif
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! WORK IN PROGRESS
!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine get_kgrid_adapt(Hk_inp,h_io,lat,io_dos)
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    type(dos_sc)                            :: dos
    real(8),allocatable                     :: eval(:)

    !parameters to describe integration k-grid
!    type(k_grid_t)                          :: k_grid
    class(kmesh_t),allocatable              :: k_grid 
    integer                                 :: ik,Nk
    real(8)                                 :: k(3)
    integer                                 :: io

    integer                                 :: idos
    character(len=3)                        :: dos_id

    type(dos_bnd_sc),allocatable            :: dos_bnd(:)
    integer                                 :: Ndos_bnd
    logical                                 :: use_bnd

    type(dos_orb_sc),allocatable            :: dos_orb(:)
    integer                                 :: Ndos_orb
    logical                                 :: use_orb
    real(8),allocatable                     :: min_eval(:)

    Call dos%init(io_dos)
    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh)

    Nk=k_grid%get_NK()
    allocate(min_eval(Nk),source=0.0d0)


    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        Call Hk_eval(Hk_inp,k,h_io,eval) 
        min_eval(ik)=minval(abs(eval))
        deallocate(eval)
    enddo
    open(newunit=io,file='dos_emin.dat')
    do ik=1,Nk
        write(io,'(4E16.8)') k_grid%get_k(ik),min_eval(ik)
    enddo
    close(io)
    ERROR STOP "CONTINUE KGRID ADAPT"
end subroutine
end module
