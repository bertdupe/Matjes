module m_input_types
use m_constants, only : pi
implicit none
private :: pi
public

type extpar_input
    real(8) ::  H(3)=0.0d0
    real(8) ::  E(3)=0.0d0
    real(8) ::  T(2)=0.0d0
    logical :: enable_H=.false.
    logical :: enable_E=.false.
    logical :: enable_T=.false.
    logical :: enable_M=.false.
    logical :: enable_u=.false.
end type

type MC_input
    !parameters what is run how often
    integer     :: n_Tsteps=1
    integer     :: n_sizerelax=1
    integer     :: n_thousand=1000
    integer     :: restart_MC_steps=0
    integer     :: T_relax=1
    integer     :: T_auto=1
    integer     :: Total_MC_Steps=1000

    logical     :: i_restart=.false.
    logical     :: print_relax=.false.
    logical     :: Cor_log=.false.
    logical     :: do_fluct=.True.

    real(8)     :: cone=pi

    logical     :: ising=.false.
    logical     :: underrelax=.false.
    logical     :: overrelax=.false.
    logical     :: equi=.false.
    logical     :: sphere=.false.
end type

type min_input
    real(8)     :: mass=1.0d0      ! mass used for minimize
    real(8)     :: dt=0.1d0         ! timestep
    integer(8)  :: N_minimization=huge(int(1,8))     !number of minimization steps
    integer     :: Efreq=100            ! frequency of printing the current energy/torque during the mimization
    real(8)     :: conv_torque=1.0d-6   !torque convergence criterion
end type


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
end type



end module
