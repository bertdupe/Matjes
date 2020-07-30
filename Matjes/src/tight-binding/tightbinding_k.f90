
module m_tightbinding_k
use m_rw_TB, only : TB_params
use m_basic_types, only : vec_point
use m_derived_types, only : cell,lattice
use m_get_position, only: calculate_distances,get_position
use m_rw_TB, only : TB_params
use m_energy_k ,only: set_dist_neigh, get_energy_kpts
use m_fftw, only: set_k_mesh, kmesh, N_kpoint
use m_highsym, only: plot_highsym_kpts
use m_fermi, only: get_fermi_k
use m_dos, only: calc_dos
use m_io_files_utils, only: close_file,open_file_write,open_file_read
use m_io_utils, only: get_parameter
implicit none
private
public :: tightbinding_k
contains


subroutine tightbinding_k(dimH,TB_pos_ext,mode_mag,my_lattice,my_motif)
    integer,intent(in)          :: dimH
    integer,intent(in)          :: TB_pos_ext(2)
    type(vec_point),intent(in)  :: mode_mag(:)
    type(lattice), intent(in)   :: my_lattice
    type(cell), intent(in)      :: my_motif

    ! N_cell is the variable that will contain the number of unit
    ! cells in the simulation
    integer :: N_cell
    real(kind=8) :: E_F !fermi energy
    real(kind=8), allocatable :: eigval(:,:) !eigen values(N_state,N_k)
    real(kind=8), allocatable :: dist_neigh(:,:) !neighbor distances for fourier transform

    integer :: nb_kpoints

    N_cell=product(shape(my_lattice%l_modes))

    !get dist_neigh
    Call get_dist_neigh(N_cell,my_lattice,my_motif,dist_neigh)

    call set_k_mesh('input',my_lattice) !this BZ should most probably be scaled by the inverse supercell size
    nb_kpoints = product(N_kpoint)
    !get kpoints on the grid
    Call get_energy_kpts(kmesh,dimH,TB_pos_ext,dist_neigh,mode_mag,eigval)
    E_F=0.0d0
    call get_fermi_k(eigval, TB_params%N_electrons,N_cell, TB_params%kt, E_F)
    Call write_dos(eigval,E_F,'dos_k.dat')
    Call plot_highsym_kpts(dimH,TB_pos_ext,dist_neigh,mode_mag,my_lattice,E_F) 
end subroutine 

subroutine get_dist_neigh(N_cell,my_lattice,my_motif,dist_neigh)
    !functions which gives the difference vectors dist_neigh for the different shell neighbors
    integer,intent(in)                    :: N_cell
    type(lattice), intent(in)             :: my_lattice
    type(cell), intent(in)                :: my_motif
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
    call calculate_distances(distances,pos,my_lattice%areal,my_lattice%dim_lat,my_lattice%boundary)
    Call set_dist_neigh(dist_neigh,distances)

end subroutine

subroutine write_dos(eigval,E_F,fname)
    use m_sort
    real(8),intent(in)  ::  eigval(:,:)
    real(8),intent(in)  ::  E_F
    character(len=*)    ::  fname

    real(8),allocatable         :: eigval_sort(:)
    integer,allocatable         :: indices(:)
    integer                     :: io_input
    logical                     :: do_dos
    real(8)                     :: E_ext(2),sigma,dE

    io_input=open_file_read('input')
    do_dos=.False.
    call get_parameter(io_input,'input','do_dos_k',do_dos)
    if(do_dos)then
        E_ext=[-1.0d0,1.0d0]
        dE=0.01d0
        sigma=0.01d0
        call get_parameter(io_input,'input','dos_sigma',sigma)
        call get_parameter(io_input,'input','dos_E_ext',2,E_ext)
        call get_parameter(io_input,'input','dos_dE',dE)
        call close_file('input',io_input)
        !sort because the calc_dos input has to be sorted
        allocate(indices(size(eigval)),source=0)
        allocate(eigval_sort(size(eigval)))
        eigval_sort=reshape(eigval,[size(eigval)])
        call sort(size(eigval_sort), eigval_sort, indices, 1.0d-5)
        Call calc_dos(eigval_sort,E_f,E_ext,dE,sigma,fname)
    endif

end subroutine

end module
