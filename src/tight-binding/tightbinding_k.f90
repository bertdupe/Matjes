module m_tightbinding_k
use m_tb_types
use m_derived_types, only : t_cell,lattice
use m_get_position, only: calculate_distances,get_position
use m_tb_params, only : TB_params
use m_fftw, only: set_k_mesh, kmesh, N_kpoint
use m_highsym, only: plot_highsym_kpts,set_highs_path
use m_fermi, only: get_fermi_k
use m_dos, only: calc_dos
use m_TB_types

use m_init_Hk
implicit none
private
public :: tightbinding_k
contains


subroutine tightbinding_k(lat)
    type(lattice), intent(in) :: lat

    ! N_cell is the variable that will contain the number of unit
    ! cells in the simulation
    integer :: N_cell
    real(8), allocatable :: eigval(:,:) !eigen values(N_state,N_k)
    real(8), allocatable :: dist_neigh(:,:) !neighbor distances for fourier transform

    type(Hk_inp_t)      :: Hk_inp

    Call get_Hk_inp(lat,TB_params%io_H,Hk_inp)

    if(TB_params%flow%dos_k)then
        if(TB_params%is_sc)then
            Call write_dos_sc(Hk_inp,TB_params%io_H,lat,TB_params%io_dos)
        else
            Call write_dos(Hk_inp,TB_params%io_H,lat,TB_params%io_dos)
        endif
    endif

    if(TB_params%flow%highs_k)then
        Call set_highs_path(lat,TB_params%io_highs)
        Call plot_highsym_kpts(Hk_inp,TB_params%io_H) 
    endif
end subroutine 

subroutine write_dos(Hk_inp,h_io,lat,io_dos)
    use m_dos2, only: dos_nc
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
    use m_dos2, only: dos_sc
    use m_constants, only : pi
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

    Call dos%init(io_dos)
    !get k_grid parameter
    Nk=product(io_dos%kgrid)
    kdiff=lat%a_sc_inv/spread(real(io_dos%kgrid),2,3)*2.0d0*pi
    k_offset=[(product(io_dos%kgrid(1:i)),i=0,2)]

    !calculate eigenvalues for each k and add to dos
    do ik=1,Nk
        kint=modulo((ik-1)/k_offset,io_dos%kgrid)
        k=matmul(kint,kdiff)
        Call Hk_evec(Hk_inp,k,h_io,eval,evec) 
        Call dos%add(eval,evec)
        deallocate(eval,evec)
    enddo

    Call dos%print(fname)
end subroutine

end module
