module m_tightbinding_k
use m_tb_types
use m_derived_types, only : lattice
use m_highsym, only: calc_highsym
use m_TB_types
use m_fermi_dos, only: fermi_dos_nc , fermi_dos_proj_nc, fermi_dos_projall_nc
use m_kgrid_int, only: k_mesh_int, kmesh_t, get_kmesh
use m_tb_k_public       !mode that contains the more efficient TB k-space type
use m_dos_k, only: calc_dos_k
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use mpi_basic, only: mpi_type

use m_Hk
use m_init_Hk
implicit none
private
public :: tightbinding_k
contains


subroutine tightbinding_k(lat,tb_par,comm)
    type(lattice), intent(in)           :: lat
    type(parameters_TB),intent(in)      :: tb_par   !input type which contains all io-related settings
    type(mpi_type),intent(in)           :: comm     !mpi communicator

    type(Hk_inp_t)              :: Hk_inp   !older, less efficient implementation of the k-space Hamiltonian
    class(H_k_base),allocatable :: Hk  !newer , more efficient implementation of the k-space Hamiltonian
    type(work_ham)              :: work         !work arrays for H_k_base Hamiltonian evaluation
    integer                     :: Hk_mode      !replace by input from tb_par

    !Initialization of k-space Hamiltonian
    !!preparation of older k-space Hamiltonian input
    Call get_Hk_inp(lat,tb_par%io_H,Hk_inp)
    !!creates newer k-space Hamiltonian from older input
    Hk_mode=tb_par%io_H%i_diag
    Call set_Hk(Hk,Hk_mode,comm)
    Call Hk%init(Hk_inp,tb_par%io_H)
    Call Hk%set_work(work)
    Call Hk_inp%destroy()  !old Hk_inp type no longer needed, might want to change Hk%init  to move Hk_inp data instead of copying

    !k-space operations that do not support MPI-parallelization
    if(comm%ismas)then
        !calculate kgrid with bands within an energy window
        if(product(tb_par%io_kmesh%grid)/=0)then
            Call print_reduced_kmesh(Hk,tb_par%io_H,tb_par%io_kmesh,lat,work)
        endif

        !writes the projection onto each orbital summed over all occupied state (T=0)  (used for terrible Jsd estimation)
        if(tb_par%flow%proj_energy) Call proj_energy(Hk,tb_par%io_H,lat,tb_par%io_dos,work)
       
        !calculates k-resolved DOS at the fermi-energy in general and including projections
        if(tb_par%flow%fermi_dos_k)then !make some more sensible input logicals...
            if(tb_par%is_sc)then
                write(error_unit,'(/A)') "Trying to get k-resolved DOS at fermi energy, but system is superconductiving."
                write(error_unit,'(A/)') "This is not implemented and will be skipped."
            else
                if(.not.tb_par%io_dos%fermi_proj_all) Call fermi_dos_nc(Hk,tb_par%io_H,lat,tb_par%io_dos,work)
                if(allocated(tb_par%io_dos%fermi_orb)) Call fermi_dos_proj_nc(Hk,tb_par%io_H,lat,tb_par%io_dos,work)
                if(tb_par%io_dos%fermi_proj_all) Call fermi_dos_projall_nc(Hk,tb_par%io_H,lat,tb_par%io_dos,work)
            endif
        endif

        !calculated bandstructure
        if(tb_par%flow%highs_k) Call calc_highsym(lat,tb_par%io_highs,Hk,work,tb_par%is_sc)
    endif

    !k-space operations that do support MPI-parallelization
    !calculated the density of states
    if(tb_par%flow%dos_k) Call calc_dos_k(Hk,tb_par%io_H,lat,tb_par%io_dos,tb_par%is_sc,work,comm)   !new version
end subroutine 

subroutine proj_energy(Hk,h_io,lat,io_dos,work)
    use m_dos_util, only: dos_nc, dos_bnd_nc, dos_orb_nc
    use m_constants, only : pi
    use m_derived_types, only: k_grid_t
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work


    !parameters to describe integration k-grid
    class(kmesh_t),allocatable      :: k_grid 
    integer                         :: ik,Nk
    real(8)                         :: k(3)

    real(8),allocatable             :: eval(:)
    complex(8),allocatable          :: evec(:,:)


    real(8),allocatable             :: state_sum(:)
    integer                         :: N_state
    integer                         :: dimH
    integer                         :: Nin,Nout         !maximal and output number of eigenvalues

    integer                         :: i,iE

    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh)
    Nk=k_grid%get_NK()

    N_state=h_io%dimH
    Nin=Hk%get_size_eval()
    dimH=Hk%get_dimH()
    allocate(state_sum(N_state),source=0.0d0)
    allocate(eval(Nin),source=0.0d0)
    allocate(evec(dimH,Nin),source=(0.0d0,0.0d0))

    do ik=1,Nk
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        Call Hk%get_evec(k,Nin,eval,evec,Nout,work) 
        do iE=1,Nout
            if(eval(iE)>0.0d0) exit !T=0
            state_sum=state_sum+real(conjg(evec(:,ie))*evec(:,ie))*eval(ie)
        enddo
    enddo
    deallocate(eval,evec)
    state_sum=state_sum/Nk
    open(newunit=i,file='proj_energy.dat')
    write(i,'(F16.8)') state_sum
    close(i)
end subroutine

subroutine print_reduced_kmesh(Hk,h_io,kmesh_io,lat,work)
    use m_derived_types, only: k_grid_t
    class(H_k_base),intent(inout)           :: Hk
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(parameters_TB_IO_kint),intent(in)  :: kmesh_io
    type(work_ham),intent(inout)            :: work

    type(k_grid_t)                          :: k_grid
    type(k_mesh_int)                        :: kmesh_int
    integer         ::  io

    Call k_grid%set(lat%a_sc_inv,kmesh_io%grid)
    Call kmesh_int%set(k_grid,Hk,h_io,kmesh_io%Ecut,work)
    open(newunit=io,file='kmesh_int.dat')
    write(io,*) kmesh_int
    close(io)
end subroutine

end module
