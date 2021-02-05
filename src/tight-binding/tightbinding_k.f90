module m_tightbinding_k
use m_tb_types
use m_derived_types, only : lattice
use m_highsym, only: plot_highsym_kpts,set_highs_path
use m_TB_types
use m_fermi_dos, only: fermi_dos_nc 
use, intrinsic :: iso_fortran_env, only : output_unit

use m_init_Hk
implicit none
private
public :: tightbinding_k
contains


subroutine tightbinding_k(lat,tb_par)
    type(lattice), intent(in)       :: lat
    type(parameters_TB),intent(in)  :: tb_par

    ! N_cell is the variable that will contain the number of unit
    ! cells in the simulation
    integer :: N_cell
    real(8), allocatable :: eigval(:,:) !eigen values(N_state,N_k)
    real(8), allocatable :: dist_neigh(:,:) !neighbor distances for fourier transform

    type(Hk_inp_t)      :: Hk_inp

    Call get_Hk_inp(lat,tb_par%io_H,Hk_inp)

    if(tb_par%flow%dos_k)then
        if(tb_par%is_sc)then
            Call write_dos_sc(Hk_inp,tb_par%io_H,lat,tb_par%io_dos)
        else
            Call write_dos_nc(Hk_inp,tb_par%io_H,lat,tb_par%io_dos)
        endif
    endif

    if(tb_par%flow%fermi_dos_k) Call fermi_dos_nc(Hk_inp,tb_par%io_H,lat,tb_par%io_dos)

    if(tb_par%flow%highs_k)then
        Call set_highs_path(lat,tb_par%io_highs)
        Call plot_highsym_kpts(Hk_inp,tb_par%io_H) 
    endif
end subroutine 

subroutine write_dos_nc(Hk_inp,h_io,lat,io_dos)
    use m_dos, only: dos_nc, dos_bnd_nc, dos_orb_nc
    use m_constants, only : pi
    use m_derived_types, only: k_grid_t
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    type(dos_nc)                            :: dos
    real(8),allocatable                     :: eval(:)
    complex(8),allocatable                  :: evec(:,:)

    !parameters to describe integration k-grid
    type(k_grid_t)                          :: k_grid
    integer                                 :: ik,Nk
    real(8)                                 :: k(3)

    character(len=3)                        :: dos_id

    type(dos_bnd_nc),allocatable            :: dos_bnd(:)
    logical                                 :: use_bnd
    integer                                 :: idos,Ndos_bnd

    type(dos_orb_nc),allocatable            :: dos_orb(:)
    logical                                 :: use_orb
    integer                                 :: Ndos_orb


    Call dos%init(io_dos)
    Call k_grid%set(lat%a_sc_inv,io_dos%kgrid)
    Nk=k_grid%get_NK()

    Ndos_bnd=0
    use_bnd=allocated(io_dos%bnd)
    if(use_bnd)then
        Ndos_bnd=size(io_dos%bnd,2)
        allocate(dos_bnd(Ndos_bnd))
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%init_bnd(io_dos,io_dos%bnd(:,idos))
        enddo
    endif
    
    Ndos_orb=0
    use_orb=allocated(io_dos%orb)
    if(use_orb)then
        Ndos_orb=size(io_dos%orb)
        allocate(dos_orb(Ndos_orb))
        do idos=1,Ndos_orb
            Call dos_orb(idos)%init_bnd(io_dos,io_dos%orb(idos),h_io%norb*h_io%nspin)
        enddo
    endif

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

    Call dos%print('dos_k.dat')
    do idos=1,Ndos_bnd
        write(dos_id,'(I0.3)') idos
        Call dos_bnd(idos)%print("dos_k_nc_bnd_"//dos_id//".dat")
    enddo

    do idos=1,Ndos_orb
        write(dos_id,'(I0.3)') io_dos%orb(idos)
        Call dos_orb(idos)%print("dos_k_nc_orb_"//dos_id//".dat")
    enddo
end subroutine

subroutine write_dos_sc(Hk_inp,h_io,lat,io_dos)
    use m_dos, only: dos_sc,dos_bnd_sc, dos_orb_sc
    use m_constants, only : pi
    use m_derived_types, only: k_grid_t
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    type(dos_sc)                            :: dos
    real(8),allocatable                     :: eval(:)
    complex(8),allocatable                  :: evec(:,:)

    !parameters to describe integration k-grid
    type(k_grid_t)                          :: k_grid
    integer                                 :: ik,Nk
    real(8)                                 :: k(3)

    integer                                 :: idos
    character(len=3)                        :: dos_id

    type(dos_bnd_sc),allocatable            :: dos_bnd(:)
    integer                                 :: Ndos_bnd
    logical                                 :: use_bnd

    type(dos_orb_sc),allocatable            :: dos_orb(:)
    integer                                 :: Ndos_orb
    logical                                 :: use_orb
    Call dos%init(io_dos)
    Call k_grid%set(lat%a_sc_inv,io_dos%kgrid)
    Nk=k_grid%get_NK()

    Ndos_bnd=0
    use_bnd=allocated(io_dos%bnd)
    if(use_bnd)then
        Ndos_bnd=size(io_dos%bnd,2)
        allocate(dos_bnd(Ndos_bnd))
        do idos=1,Ndos_bnd
            Call dos_bnd(idos)%init_bnd(io_dos,io_dos%bnd(:,idos))
        enddo
    endif

    Ndos_orb=0
    use_orb=allocated(io_dos%orb)
    if(use_orb)then
        Ndos_orb=size(io_dos%orb)
        allocate(dos_orb(Ndos_orb))
        do idos=1,Ndos_orb
            Call dos_orb(idos)%init_bnd(io_dos,io_dos%orb(idos),h_io%norb*h_io%nspin)
        enddo
    endif

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

    Call dos%print('dos_k_sc.dat')
    do idos=1,Ndos_bnd
        write(dos_id,'(I0.3)') idos
        Call dos_bnd(idos)%print("dos_k_sc_bnd_"//dos_id//".dat")
    enddo
    do idos=1,Ndos_orb
        write(dos_id,'(I0.3)') io_dos%orb(idos)
        Call dos_orb(idos)%print("dos_k_sc_orb_"//dos_id//".dat")
    enddo
end subroutine

end module
