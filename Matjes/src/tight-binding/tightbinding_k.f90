
module m_tightbinding_k
use m_basic_types, only : vec_point
use m_tb_types
use m_derived_types, only : t_cell,lattice
use m_get_position, only: calculate_distances,get_position
use m_tb_params, only : TB_params
use m_energy_k ,only: set_dist_neigh, get_energy_kpts
use m_fftw, only: set_k_mesh, kmesh, N_kpoint
use m_highsym, only: plot_highsym_kpts,set_highs_path
use m_fermi, only: get_fermi_k
use m_dos, only: calc_dos
use m_TB_types
implicit none
private
public :: tightbinding_k
contains


subroutine tightbinding_k(h_par,mode_mag,my_lattice,my_motif)
    type(parameters_TB_Hsolve),intent(in)     ::  h_par
    type(vec_point),intent(in)  :: mode_mag(:)
    type(lattice), intent(in)   :: my_lattice
    type(t_cell), intent(in)      :: my_motif

    ! N_cell is the variable that will contain the number of unit
    ! cells in the simulation
    integer :: N_cell
    real(kind=8) :: E_F !fermi energy
    real(kind=8), allocatable :: eigval(:,:) !eigen values(N_state,N_k)
    real(kind=8), allocatable :: dist_neigh(:,:) !neighbor distances for fourier transform

    N_cell=product(shape(my_lattice%ordpar%l_modes))

    !get dist_neigh
    Call get_dist_neigh(N_cell,my_lattice,my_motif,dist_neigh)

    if(TB_params%flow%dos_k.or.TB_params%flow%dos_k)then
        !get kpoints on the grid
        call set_k_mesh('input',my_lattice)
        Call get_energy_kpts(kmesh,h_par,dist_neigh,mode_mag,eigval)
    endif

    E_F=TB_params%io_ef%E_F_in
    if(tb_params%flow%fermi_k)then
        call get_fermi_k(eigval, TB_params%io_EF%N_electrons,N_cell, TB_params%io_EF%kt, E_F)
    endif
    if(TB_params%flow%dos_k)then
        Call write_dos(eigval,TB_params%io_dos,'dos_k.dat')
    endif
    if(TB_params%flow%highs_k)then
        Call set_highs_path(my_lattice,TB_params%io_highs)
        Call plot_highsym_kpts(h_par,dist_neigh,mode_mag,E_F) 
    endif
end subroutine 

subroutine get_dist_neigh(N_cell,my_lattice,my_motif,dist_neigh)
    !functions which gives the difference vectors dist_neigh for the different shell neighbors
    integer,intent(in)                    :: N_cell
    type(lattice), intent(in)             :: my_lattice
    type(t_cell), intent(in)                :: my_motif
    real(kind=8), allocatable,intent(out) :: dist_neigh(:,:)

    real(kind=8)              :: pos(3,N_cell)
    real(kind=8)              :: distances(3,N_cell)
    real(kind=8), allocatable :: start_positions(:,:,:,:,:)
    
    pos=0.0d0
    distances=0.0d0
    allocate( start_positions(3, my_lattice%dim_lat(1), my_lattice%dim_lat(2), my_lattice%dim_lat(3), size(my_motif%atomic)) )
    call get_position( start_positions, my_lattice%dim_lat, my_lattice%areal, my_motif )
    pos=reshape( start_positions, [3, N_cell] )
    deallocate(start_positions)
    call calculate_distances(distances,pos,my_lattice%areal,my_lattice%dim_lat,my_lattice%periodic)
    Call set_dist_neigh(dist_neigh,distances)

end subroutine

subroutine write_dos(eigval,io_dos,fname)
    use m_sort
    real(8),intent(in)  ::  eigval(:,:)
    type(parameters_TB_IO_DOS),intent(in)  :: io_dos
    character(len=*)    ::  fname

    real(8),allocatable         :: eigval_sort(:)
    integer,allocatable         :: indices(:)

    !sort because the calc_dos input has to be sorted
    allocate(indices(size(eigval)),source=0)
    allocate(eigval_sort(size(eigval)))
    eigval_sort=reshape(eigval,[size(eigval)])
    call sort(size(eigval_sort), eigval_sort, indices, 1.0d-5)
    Call calc_dos(eigval_sort,io_dos,fname)

end subroutine

end module
