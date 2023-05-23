module m_GNEB
use m_H_public,only: t_H
use m_gneb_utils
use m_GNEB_io,only: GNEB_input  
use m_io_gneb,only: write_path,write_en
use m_write_spin, only: WriteSpinAndCorrFile
use m_createspinfile, only: CreateSpinFile
use m_derived_types, only: lattice,io_parameter
use m_initialize_path, only: path_initialization
use m_spline, only: hermite_fit,spline_hermite_val
use, intrinsic :: iso_fortran_env, only : error_unit 
use mpi_basic
use m_hamiltonian_collection, only: hamiltonian
implicit none
private
public GNEB

contains
subroutine GNEB(my_lattice,io_simu,H,com)
    type(lattice), intent(inout)        :: my_lattice
    type(io_parameter), intent(inout)   :: io_simu
    type(hamiltonian),intent(inout)     :: H
    type(mpi_type),intent(in)           :: com
    write(6,'(a)') 'entering into the GNEB routine'

    if(com%ismas)then
        if(com%Np>1) write(error_unit,'(//A/A//)') "WARNING, USING MPI-PARALLELIZATION WHICH IS NOT IMPLEMENTED FOR GNEB CALCULATION","this calculation will only run the the master thread" 
        Call GNEB_run(my_lattice,io_simu,H)
    endif

end subroutine

subroutine GNEB_run(my_lattice,io_simu,H)
    type(lattice), intent(in)       :: my_lattice
    type(io_parameter), intent(in)  :: io_simu
    type(hamiltonian),intent(inout) :: H
    !internal variable
    type(GNEB_input)            :: io_gneb
    type(lattice),allocatable   :: images(:)
    real(8), allocatable        :: spinsp(:,:)
    real(8), allocatable        :: rx(:),ene(:),dene(:),xx(:),yy(:),dyy(:),c(:,:)
    real(8) :: rx0(1),drx0,ene0(1),dene0(1)
    integer :: imax,i,ci,size_order,N_cell
    logical :: gra_log
    integer :: nim,io_frequency
    
    N_cell=my_lattice%ncell
    size_order=my_lattice%M%dim_mode
    gra_log=io_simu%io_Xstruct
    ci=1
    Call io_gneb%read_file()
    nim=io_gneb%nim

    allocate(images(nim))
    do i=1,nim
        Call my_lattice%copy(images(i))
    enddo
    allocate(rx(nim),ene(nim),dene(nim),source=0.0d0)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! initialization of the path
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    call path_initialization(images,io_simu,io_gneb,H)
    call write_path(images)
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Running the Gneb calculation
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    if (io_gneb%do_gneb) then
       write (6,'(a)') "GNEB calculation in progress..."
       call find_path(nim,N_cell,io_gneb%mepftol,rx,ene,dene,images,io_gneb,H)
       write (6,'(a)') "Done!"
    end if
    if (io_gneb%do_gneb_ci) then
       write (6,'(a)') "CI-GNEB calculation in progress..."
       call find_path(nim,N_cell,io_gneb%mepftol_ci,rx,ene,dene,images,io_gneb,H,ci)
       write(6,'(a,I3)') 'ci:',ci
       write (6,'(a)') "Done!"
    end if
    call write_en(nim,rx,ene,dene,rx(nim),'en_path.out',io_gneb%do_norm_rx)
    call write_path(images)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! Do further things with GNEB results
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          
    allocate(xx(io_gneb%sample_num),yy(io_gneb%sample_num),dyy(io_gneb%sample_num),c(4,nim),source=0.d0)
    
    call hermite_fit(nim,io_gneb%sample_num,rx,ene,dene,xx,yy,dyy,c)
    call write_en(io_gneb%sample_num,xx,yy,dyy,rx(nim),'enfit_path.out',io_gneb%do_norm_rx)
          
    allocate(spinsp(size_order,N_cell),source=0.0d0)
    if (io_gneb%do_gneb_ci) then
        rx0(1) = rx(ci)
        spinsp=images(ci)%M%modes_v
    else
        call find_SP(nim,rx,c,imax,drx0)
        rx0(1) = rx(imax)+drx0
        if (imax==nim) then
            spinsp=images(nim)%M%modes_v
        else
            call find_SP_conf(images(imax),images(imax+1),rx(imax),rx(imax+1),rx(imax)+drx0,spinsp)
        end if
    end if
    
    call spline_hermite_val(nim,rx,c,rx0(1),ene0(1),dene0(1))
    call write_en(1,rx0,ene0,dene0,rx(nim),'ensp_path.out',io_gneb%do_norm_rx)
    call WriteSpinAndCorrFile(spinsp,'image-GNEB-saddlepoint.dat')
    call CreateSpinFile(spinsp,'povray-GNEB-saddlepoint.dat')
    
    deallocate(rx,ene,dene,xx,yy,dyy,c)
end subroutine GNEB_run 
end module
