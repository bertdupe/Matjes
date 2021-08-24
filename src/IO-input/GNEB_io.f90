module m_GNEB_io
use m_io_minimize, only: min_input
implicit none
private
public GNEB_input

type GNEB_input  
    character(len=35)   :: momfile_i='momfile_i'            !< File name for initial magn. configuration
    character(len=35)   :: momfile_f='momfile_f'            !< File name for final magn. configuration
    character(len=35)   :: restartfile_if='nofile'          !< Name of the file with initial and final states needed to restart energy minimization calculations
    character(len=35)   :: restartfile_path='nofile'        !< Name of the file with the path needed to restart GNEB calculations
!   real(8)             :: ind_tol  !seems to be obsolete
    real(8)             :: amp_rnd=0.0                      !< Amplitude of random perturbation of the components of magnetic moments
    real(8)             :: amp_rnd_path=0.0
!   integer             :: minalgo  !seems to be obsolete   !< Minimization algorithm (1=VPO,...)
    integer             :: minitrmax=10000000               !< Maximum number of iterations in energy minimization procedure
    real(8)             :: minftol=1.0d-9                   !< Convergence criterion in mRy
    integer             :: mintraj_step=10                  !< Save configuration every 'mintraj_step' step during energy minimization
    real(8)             :: dt=1.0d-2                        !< VPO timestep
    real(8)             :: mass=1.0d0                       !< VPO mass
    real(8)             :: spring=0.5d0                     !< Magnitude of the spring constant in GNEB calculations
    integer(8)          :: itrmax=huge(int(1,8))            !< Maximum number of iterations in MEP finding procedure
    real(8)             :: mepftol=1.0d-3                   !< Convergence criterion in the GNEB method, in mRy
    real(8)             :: mepftol_ci=1.0d-5                !< Convergence criterion in the CI-GNEB method, in mRy
    integer             :: every=10                         !< Save configuration every 'meptraj_step' step during MEP finding procedure
    logical             :: do_gneb=.False.                  !< Do GNEB calculations (Y/N)
    logical             :: do_gneb_ci=.False.               !< Do CI-GNEB calculations (Y/N)
    character(len=1)    :: do_norm_rx ='N'                  !< Normalize reaction coordinate (Y/N)
    character(len=1)    :: en_zero ='N'                     !< Level of zero energy. 'I' - initial state; 'F' - final state; 'N' - 0.0
    integer             :: sample_num=500                   !< Number of samples in the interpolated curve
    integer             :: nim=10
!parameters for the initialization of the path
    integer             :: initpath=-1                        !< Path initialization method (1=geodesic, 2=read from file,...,-1 ignore)
    logical             :: read_path=.False.                  !< read the entire path and do not initialize with geodesic spins
    logical             :: read_outer=.True.                  !< read the outer images 
    integer             :: min_type=1                         !< minimalization routine to minimize outer images (0=none,1=minimize_infdamp, 2=minimize)
    type(min_input)     :: io_min                        !input for the minimize_infdamp/minimize routine
contains
    procedure   :: read_file => rw_gneb
end type

contains
subroutine rw_gneb(gneb_io,fname_in)
    use m_io_files_utils, only : open_file_read,close_file
    use m_io_utils,only: get_parameter
    class(GNEB_input),intent(out)        :: gneb_io
    character(*),intent(in),optional    :: fname_in
    !internal
    character(*),parameter              :: fname_default='input'
    character(:), allocatable           :: fname
    integer                             :: io_param

    if(present(fname_in))then
        fname=fname_in
    else
        fname=fname_default
    endif
    io_param=open_file_read(fname)
    call get_parameter(io_param,fname,'momfile_i',gneb_io%momfile_i)
    call get_parameter(io_param,fname,'momfile_f',gneb_io%momfile_f)
    call get_parameter(io_param,fname,'restartfile_if',gneb_io%restartfile_if)
    call get_parameter(io_param,fname,'restartfile_path',gneb_io%restartfile_path)
    call get_parameter(io_param,fname,'spring',gneb_io%spring)
    call get_parameter(io_param,fname,'initpath',gneb_io%initpath)
    call get_parameter(io_param,fname,'amp_rnd',gneb_io%amp_rnd)
    call get_parameter(io_param,fname,'amp_rnd_path',gneb_io%amp_rnd_path)
    call get_parameter(io_param,fname,'mintraj_step',gneb_io%mintraj_step)
    call get_parameter(io_param,fname,'min_ftol',gneb_io%minftol)
    call get_parameter(io_param,fname,'mep_itrmax',gneb_io%itrmax)
    call get_parameter(io_param,fname,'meptraj_step',gneb_io%every)
    call get_parameter(io_param,fname,'mep_ftol',gneb_io%mepftol)
    call get_parameter(io_param,fname,'mep_ftol_ci',gneb_io%mepftol_ci)
    call get_parameter(io_param,fname,'do_gneb_simple',gneb_io%do_gneb)
    call get_parameter(io_param,fname,'do_gneb_ci',gneb_io%do_gneb_ci)
    call get_parameter(io_param,fname,'do_norm_rx',gneb_io%do_norm_rx)
    call get_parameter(io_param,fname,'en_zero',gneb_io%en_zero)
    call get_parameter(io_param,fname,'vpo_dt',gneb_io%dt)
    call get_parameter(io_param,fname,'vpo_mass',gneb_io%mass)
    call get_parameter(io_param,fname,'sample_num',gneb_io%sample_num)
    call get_parameter(io_param,fname,'nim',gneb_io%nim)

    !gneb_io%init_path to support older input format
    select case(gneb_io%initpath)
    case(-1)
        continue !do nothing
    case(1)
        gneb_io%read_path=.False.      
        gneb_io%read_outer=.True.      
        gneb_io%min_type=0            
    case(2)
        gneb_io%read_path=.True.      
        gneb_io%read_outer=.True.      
        gneb_io%min_type=2            
    case(3)
        gneb_io%read_path=.True.      
        gneb_io%read_outer=.True.      
        gneb_io%min_type=1            
    case default
        gneb_io%read_path=.False.      
        gneb_io%read_outer=.True.      
        gneb_io%min_type=1            
    end select
    call get_parameter(io_param,fname,'gnebinit_read_path',gneb_io%read_path)
    call get_parameter(io_param,fname,'gnebinit_read_outer',gneb_io%read_outer)
    call get_parameter(io_param,fname,'gnebinit_min_type',gneb_io%min_type)
    Call gneb_io%io_min%read_file(io_param,fname)

    call close_file(fname,io_param)
end subroutine
end module
