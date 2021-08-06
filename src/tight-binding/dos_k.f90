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

interface calc_dos_k
    module procedure calc_dos_k_old
    module procedure calc_dos_k_new
    module procedure calc_dos_k_mpi
end interface

interface write_dos_nc
    module procedure write_dos_nc_old
    module procedure write_dos_nc_new
end interface

interface write_dos_sc
    module procedure write_dos_sc_old
    module procedure write_dos_sc_new
end interface

interface init_dos_bnd
    module procedure init_dos_bnd_sc
    module procedure init_dos_bnd_nc
end interface

interface init_dos_orb
    module procedure init_dos_orb_sc
    module procedure init_dos_orb_nc
end interface

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! NEW MPI Hk implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_dos_k_mpi(Hk,io_H,lat,io_dos,is_sc,work,comm)
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: io_H
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work
    logical,intent(in)                      :: is_sc
    type(mpi_type),intent(in)               :: comm     !mpi communicator

    if(is_sc)then
        STOP "IMPLEMENT FASFASFA"
!        Call write_dos_sc_new(Hk,io_H,lat,io_dos,work)
    else
        Call write_dos_nc_mpi(Hk,io_H,lat,io_dos,work,comm)
    endif
end subroutine


subroutine write_dos_nc_mpi(Hk,h_io,lat,io_dos,work,comm)
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work
    type(mpi_type),intent(in)               :: comm

    type(dos_nc)                            :: dos

    integer                                 :: Nin,Nout, dimH ,ierr
    real(8),allocatable                     :: eval(:)
    complex(8),allocatable                  :: evec(:,:)

    !parameters to describe integration k-grid
    class(kmesh_t),allocatable              :: k_grid 
    integer                                 :: ik,Nk, idos
    integer                                 :: bnd(2)   !boundaries for k=loop
    real(8)                                 :: k(3)

    character(len=3)                        :: dos_id

    !parameters for projection on a contiguous set of orbitals somewhere in the unit-cell 
    type(dos_bnd_nc),allocatable            :: dos_bnd(:)
    logical                                 :: use_bnd
    integer                                 :: Ndos_bnd

    !parameters for projection on an orbital considering all supercell periodic images
    type(dos_orb_nc),allocatable            :: dos_orb(:)
    logical                                 :: use_orb
    integer                                 :: Ndos_orb

    Call dos%init(io_dos)
    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh)
    Nk=k_grid%get_NK()

    !initialize data for projection on a contiguous set of orbitals somewhere in the unit-cell 
    Call init_dos_bnd(io_dos,dos_bnd,use_bnd,Ndos_bnd)
    
    !initialize data for projection on an orbital considering all supercell periodic images 
    Call init_dos_orb(io_dos,dos_orb,use_orb,Ndos_orb,h_io%norb*h_io%nspin)

    !prepare data
    Nin=Hk%get_size_eval()
    allocate(eval(Nin),source=0.0d0)
    if(use_bnd.or.use_orb)then
        dimH=Hk%get_dimH()
        allocate(evec(dimH,Nin),source=(0.0d0,0.0d0))
    endif

    !calculate eigenvalues for each k and add to dos
    Call comm%get_loop_bnd(Nk,bnd)
    do ik=bnd(1),bnd(2)
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        if(use_bnd.or.use_orb)then
            Call Hk%get_evec(k,Nin,eval,evec,Nout,work) 
        else
            Call Hk%get_eval(k,Nin,eval,Nout,work) 
        endif
        Call dos%add(eval)
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_orb
            Call dos_orb(idos)%add(eval,evec)
        enddo
    enddo

    Call dos_mpi_sum(dos,dos_bnd,dos_orb, comm)

    if(comm%ismas) Call print_dos_output(dos,dos_bnd,dos_orb,"_nc")
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! NEW Hk implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_dos_k_new(Hk,io_H,lat,io_dos,is_sc,work)
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: io_H
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work
    logical,intent(in)                      :: is_sc

    if(is_sc)then
        Call write_dos_sc_new(Hk,io_H,lat,io_dos,work)
    else
        Call write_dos_nc_new(Hk,io_H,lat,io_dos,work)
    endif
end subroutine

subroutine write_dos_nc_new(Hk,h_io,lat,io_dos,work)
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work

    type(dos_nc)                            :: dos

    integer                                 :: Nin,Nout, dimH
    real(8),allocatable                     :: eval(:)
    complex(8),allocatable                  :: evec(:,:)

    !parameters to describe integration k-grid
    class(kmesh_t),allocatable              :: k_grid 
    integer                                 :: ik,Nk, idos
    real(8)                                 :: k(3)

    character(len=3)                        :: dos_id

    !parameters for projection on a contiguous set of orbitals somewhere in the unit-cell 
    type(dos_bnd_nc),allocatable            :: dos_bnd(:)
    logical                                 :: use_bnd
    integer                                 :: Ndos_bnd

    !parameters for projection on an orbital considering all supercell periodic images
    type(dos_orb_nc),allocatable            :: dos_orb(:)
    logical                                 :: use_orb
    integer                                 :: Ndos_orb

    Call dos%init(io_dos)
    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh)
    Nk=k_grid%get_NK()

    !initialize data for projection on a contiguous set of orbitals somewhere in the unit-cell 
    Call init_dos_bnd(io_dos,dos_bnd,use_bnd,Ndos_bnd)
    
    !initialize data for projection on an orbital considering all supercell periodic images 
    Call init_dos_orb(io_dos,dos_orb,use_orb,Ndos_orb,h_io%norb*h_io%nspin)

    !prepare data
    Nin=Hk%get_size_eval()
    allocate(eval(Nin),source=0.0d0)
    if(use_bnd.or.use_orb)then
        dimH=Hk%get_dimH()
        allocate(evec(dimH,Nin),source=(0.0d0,0.0d0))
    endif

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        if(use_bnd.or.use_orb)then
            Call Hk%get_evec(k,Nin,eval,evec,Nout,work) 
        else
            Call Hk%get_eval(k,Nin,eval,Nout,work) 
        endif
        Call dos%add(eval)
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_orb
            Call dos_orb(idos)%add(eval,evec)
        enddo
    enddo

    Call print_dos_output(dos,dos_bnd,dos_orb,"_nc")
end subroutine

subroutine write_dos_sc_new(Hk,h_io,lat,io_dos,work)
    class(H_k_base),intent(inout)           :: Hk
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    type(work_ham),intent(inout)            :: work

    type(dos_sc)                            :: dos

    integer                                 :: Nin,Nout, dimH
    real(8),allocatable                     :: eval(:)
    complex(8),allocatable                  :: evec(:,:)

    !parameters to describe integration k-grid
    class(kmesh_t),allocatable              :: k_grid 
    integer                                 :: ik,Nk
    real(8)                                 :: k(3)

    integer                                 :: idos
    character(len=3)                        :: dos_id

    !parameters for projection on a contiguous set of orbitals somewhere in the unit-cell 
    type(dos_bnd_sc),allocatable            :: dos_bnd(:)
    integer                                 :: Ndos_bnd
    logical                                 :: use_bnd

    !parameters for projection on an orbital considering all supercell periodic images
    type(dos_orb_sc),allocatable            :: dos_orb(:)
    integer                                 :: Ndos_orb
    logical                                 :: use_orb

    Call dos%init(io_dos)
    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh)

    Nk=k_grid%get_NK()

    !initialize data for projection on a contiguous set of orbitals somewhere in the unit-cell 
    Call init_dos_bnd(io_dos,dos_bnd,use_bnd,Ndos_bnd)
    
    !initialize data for projection on an orbital considering all supercell periodic images 
    Call init_dos_orb(io_dos,dos_orb,use_orb,Ndos_orb,h_io%norb*h_io%nspin)

    !prepare data
    Nin=Hk%get_size_eval()
    dimH=Hk%get_dimH()
    allocate(eval(Nin),source=0.0d0)
    allocate(evec(dimH,Nin),source=(0.0d0,0.0d0))

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        Call Hk%get_evec(k,Nin,eval,evec,Nout,work) 
        Call dos%add(eval,evec)
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_orb
            Call dos_orb(idos)%add(eval,evec)
        enddo
    enddo

    Call print_dos_output(dos,dos_bnd,dos_orb,"_sc")
end subroutine



!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! OLD Hk implementation
!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine calc_dos_k_old(Hk_inp,io_H,lat,io_dos,is_sc)
    type(Hk_inp_t),intent(in)   :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: io_H
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    logical,intent(in)          :: is_sc

    !Call get_kgrid_adapt(Hk_inp,tb_par%io_H,lat,tb_par%io_dos)
    if(is_sc)then
        Call write_dos_sc(Hk_inp,io_H,lat,io_dos)
    else
        Call write_dos_nc(Hk_inp,io_H,lat,io_dos)
    endif
end subroutine

subroutine write_dos_nc_old(Hk_inp,h_io,lat,io_dos)
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    type(dos_nc)                            :: dos
    real(8),allocatable                     :: eval(:)
    complex(8),allocatable                  :: evec(:,:)

    !parameters to describe integration k-grid
    class(kmesh_t),allocatable              :: k_grid 
    integer                                 :: ik,Nk, idos
    real(8)                                 :: k(3)

    character(len=3)                        :: dos_id

    !parameters for projection on a contiguous set of orbitals somewhere in the unit-cell 
    type(dos_bnd_nc),allocatable            :: dos_bnd(:)
    logical                                 :: use_bnd
    integer                                 :: Ndos_bnd

    !parameters for projection on an orbital considering all supercell periodic images
    type(dos_orb_nc),allocatable            :: dos_orb(:)
    logical                                 :: use_orb
    integer                                 :: Ndos_orb

    Call dos%init(io_dos)
    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh)
    Nk=k_grid%get_NK()

    !initialize data for projection on a contiguous set of orbitals somewhere in the unit-cell 
    Call init_dos_bnd(io_dos,dos_bnd,use_bnd,Ndos_bnd)
    
    !initialize data for projection on an orbital considering all supercell periodic images 
    Call init_dos_orb(io_dos,dos_orb,use_orb,Ndos_orb,h_io%norb*h_io%nspin)

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        if(use_bnd.or.use_orb)then
            Call Hk_evec(Hk_inp,k,h_io,eval,evec) 
        else
            Call Hk_eval(Hk_inp,k,h_io,eval) 
        endif
        Call dos%add(eval)
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_orb
            Call dos_orb(idos)%add(eval,evec)
        enddo
        deallocate(eval)
        if(use_bnd) deallocate(evec)
    enddo

    Call print_dos_output(dos,dos_bnd,dos_orb,"_nc")
end subroutine

subroutine write_dos_sc_old(Hk_inp,h_io,lat,io_dos)
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    type(dos_sc)                            :: dos
    real(8),allocatable                     :: eval(:)
    complex(8),allocatable                  :: evec(:,:)

    !parameters to describe integration k-grid
    class(kmesh_t),allocatable              :: k_grid 
    integer                                 :: ik,Nk
    real(8)                                 :: k(3)

    integer                                 :: idos
    character(len=3)                        :: dos_id

    !parameters for projection on a contiguous set of orbitals somewhere in the unit-cell 
    type(dos_bnd_sc),allocatable            :: dos_bnd(:)
    integer                                 :: Ndos_bnd
    logical                                 :: use_bnd

    !parameters for projection on an orbital considering all supercell periodic images
    type(dos_orb_sc),allocatable            :: dos_orb(:)
    integer                                 :: Ndos_orb
    logical                                 :: use_orb

    Call dos%init(io_dos)
    Call get_kmesh(k_grid,lat,io_dos%kgrid,io_dos%fname_kmesh)

    Nk=k_grid%get_NK()

    !initialize data for projection on a contiguous set of orbitals somewhere in the unit-cell 
    Call init_dos_bnd(io_dos,dos_bnd,use_bnd,Ndos_bnd)
    
    !initialize data for projection on an orbital considering all supercell periodic images 
    Call init_dos_orb(io_dos,dos_orb,use_orb,Ndos_orb,h_io%norb*h_io%nspin)

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        Call Hk_evec(Hk_inp,k,h_io,eval,evec) 
        Call dos%add(eval,evec)
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        do idos=1,Ndos_orb
            Call dos_orb(idos)%add(eval,evec)
        enddo
        deallocate(eval,evec)
    enddo

    Call print_dos_output(dos,dos_bnd,dos_orb,"_sc")
end subroutine




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!   MINOR HELP ROUTINES
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine print_dos_output(dos,dos_bnd,dos_orb,specifier)
    class(dos_t),intent(inout)              :: dos
    class(dos_t),intent(inout)              :: dos_bnd(:)
    class(dos_t),intent(inout)              :: dos_orb(:)
    character(len=*),intent(in)             :: specifier
    logical             :: alloc
    integer             :: i
    character(len=3)    :: dos_id
   
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
end subroutine

subroutine dos_mpi_sum(dos,dos_bnd,dos_orb, comm)
    class(dos_t),intent(inout)              :: dos
    class(dos_t),intent(inout)              :: dos_bnd(:)
    class(dos_t),intent(inout)              :: dos_orb(:)
    type(mpi_type),intent(in)               :: comm
    integer     ::  i

    Call dos%mpi_sum(comm)
    do i=1,size(dos_bnd)
        Call dos_bnd(i)%mpi_sum(comm)
    enddo
    do i=1,size(dos_orb)
        Call dos_orb(i)%mpi_sum(comm)
    enddo
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
