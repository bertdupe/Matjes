module m_tightbinding_k
use m_tb_types
use m_derived_types, only : t_cell,lattice
use m_get_position, only: calculate_distances,get_position
use m_fftw, only: set_k_mesh, kmesh, N_kpoint
use m_highsym, only: plot_highsym_kpts,set_highs_path
use m_fermi, only: get_fermi_k
use m_TB_types
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
            Call write_dos(Hk_inp,tb_par%io_H,lat,tb_par%io_dos)
        endif
    endif

    if(tb_par%flow%highs_k)then
        Call set_highs_path(lat,tb_par%io_highs)
        Call plot_highsym_kpts(Hk_inp,tb_par%io_H) 
    endif
end subroutine 

subroutine write_dos(Hk_inp,h_io,lat,io_dos)
    use m_dos, only: dos_nc, dos_bnd_nc
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

    type(dos_bnd_nc),allocatable            :: dos_bnd(:)
    logical                                 :: use_bnd
    integer                                 :: idos,Ndos
    character(len=3)                        :: bnd_id


    Call dos%init(io_dos)
    Call k_grid%set(lat%a_sc_inv,io_dos%kgrid)
    Nk=k_grid%get_NK()

    Ndos=0
    use_bnd=allocated(io_dos%bnd)
    if(use_bnd)then
        Ndos=size(io_dos%bnd,2)
        allocate(dos_bnd(Ndos))
        do idos=1,Ndos
            Call dos_bnd(idos)%init_bnd(io_dos,io_dos%bnd(:,idos))
        enddo
    endif

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        if(use_bnd)then
            Call Hk_evec(Hk_inp,k,h_io,eval,evec) 
        else
            Call Hk_eval(Hk_inp,k,h_io,eval) 
        endif
        Call dos%add(eval)
        do idos=1,Ndos
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        deallocate(eval)
        if(use_bnd) deallocate(evec)
    enddo

    Call dos%print('dos_k.dat')
    do idos=1,Ndos
        write(bnd_id,'(I0.3)') idos
        Call dos_bnd(idos)%print("dos_k_nc_bnd_"//bnd_id//".dat")
    enddo
end subroutine

subroutine write_dos_sc(Hk_inp,h_io,lat,io_dos)
    use m_dos, only: dos_sc,dos_bnd_sc
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

    type(dos_bnd_sc),allocatable            :: dos_bnd(:)
    integer                                 :: idos,Ndos
    character(len=3)                        :: bnd_id

    Call dos%init(io_dos)
    Call k_grid%set(lat%a_sc_inv,io_dos%kgrid)
    Nk=k_grid%get_NK()

    Ndos=0
    if(allocated(io_dos%bnd))then
        Ndos=size(io_dos%bnd,2)
        allocate(dos_bnd(Ndos))
        do idos=1,Ndos
            Call dos_bnd(idos)%init_bnd(io_dos,io_dos%bnd(:,idos))
        enddo
    endif

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        if(io_dos%print_kint) write(output_unit,'(2(AI6))') 'start dosk', ik,' of',Nk
        k=k_grid%get_K(ik)
        Call Hk_evec(Hk_inp,k,h_io,eval,evec) 
        Call dos%add(eval,evec)
        do idos=1,Ndos
            Call dos_bnd(idos)%add(eval,evec)
        enddo
        deallocate(eval,evec)
    enddo

    Call dos%print('dos_k_sc.dat')
    do idos=1,Ndos
        write(bnd_id,'(I0.3)') idos
        Call dos_bnd(idos)%print("dos_k_sc_bnd_"//bnd_id//".dat")
    enddo
end subroutine

end module
