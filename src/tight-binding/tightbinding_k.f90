module m_tightbinding_k
use m_tb_types
use m_derived_types, only : t_cell,lattice
use m_get_position, only: calculate_distances,get_position
use m_fftw, only: set_k_mesh, kmesh, N_kpoint
use m_highsym, only: plot_highsym_kpts,set_highs_path
use m_fermi, only: get_fermi_k
use m_TB_types

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
    use m_dos, only: dos_nc
    use m_constants, only : pi
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    character(len=*),parameter              :: fname='dos_k.dat'

    type(dos_nc)                            :: dos

    !parameters to describe integration k-grid
    integer                                 :: i
    integer                                 :: ik,Nk
    real(8)                                 :: k(3)
    integer                                 :: kint(3),k_offset(3)
    real(8)                                 :: kdiff(3,3)
    real(8),allocatable                     :: eval(:)

    Call dos%init(io_dos)

    !get k_grid parameter
    Nk=product(io_dos%kgrid)
    kdiff=lat%a_sc_inv/spread(real(io_dos%kgrid),2,3)*2.0d0*pi
    k_offset=[(product(io_dos%kgrid(1:i)),i=0,2)]

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        kint=modulo((ik-1)/k_offset,io_dos%kgrid)
        k=matmul(kint,kdiff)
        Call Hk_eval(Hk_inp,k,h_io,eval) 
        Call dos%add(eval)
        deallocate(eval)
    enddo

    Call dos%print(fname)
end subroutine

subroutine write_dos_sc(Hk_inp,h_io,lat,io_dos)
    use m_dos, only: dos_sc
    use m_constants, only : pi
    use m_derived_types, only: k_grid_t
    type(Hk_inp_t),intent(in)               :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    type(lattice), intent(in)               :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos
    character(len=*),parameter              :: fname='dos_k_sc.dat'

    type(dos_sc)                            :: dos

    !parameters to describe integration k-grid
    integer                                 :: i
    integer                                 :: ik,Nk
    real(8)                                 :: k(3)
    integer                                 :: kint(3),k_offset(3)
    real(8)                                 :: kdiff(3,3)
    real(8),allocatable                     :: eval(:)
    complex(8),allocatable                  :: evec(:,:)
    type(k_grid_t)                          :: k_grid

    Call dos%init(io_dos)
    Call k_grid%set(lat%a_sc_inv,io_dos%kgrid)
    Nk=k_grid%get_NK()

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        k=k_grid%get_K(ik)
        Call Hk_evec(Hk_inp,k,h_io,eval,evec) 
        Call dos%add(eval,evec)
        deallocate(eval,evec)
    enddo

    Call dos%print(fname)
end subroutine

end module
