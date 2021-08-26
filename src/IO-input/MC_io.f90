module m_MC_io
use m_constants, only : pi
use mpi_basic,only: mpi_type
implicit none

public :: MC_input

type MC_input
    !parameters what is run how often
    !THE MC parameter are used of MonteCarlo and parallel-tempering, but the meaning ofter differs and definitely should be unified by someone
    integer     :: n_Tsteps=1               !number of temperatures considered between initial and final temperature
    integer     :: n_sizerelax=1            !number of times the state is printed during the initial relaxation (if print_relax)
    integer     :: N_relaxation=1000        !number of outer loops for relaxation
    integer     :: T_relax=1                !number of MCsteps in inner relaxation (times number of spins)
    integer     :: T_auto=1                 !number of autocorrelation steps (afterwards multiplied with number of spins)
    integer     :: Total_MC_Steps=1000      !number of averaging steps
    integer     :: N_Topt=1                 !number temperature optimization steps in parallel tempering

    logical     :: print_relax=.false.      !writes information about the states during the initial relaxation to judge success of intial relaxation
    logical     :: Cor_log=.false.          !calculate correlation? (probably rather obsolete)
    logical     :: expval_save=.false.      !save the expectation values
    logical     :: expval_read=.false.      !read expectation values from previous calculation
    real(8)     :: cone=pi                  !initialization variable for equi and sphere spin selection
    !method for sampling of test state 
    !(in which direction should the spin of the chosen site point?)
    logical     :: ising=.false.            !Ising samling (reverses spin)
    logical     :: underrelax=.false.       !underrelaxation (update magnetization according to M=M-Mx(MxB))
    logical     :: overrelax=.false.        !overrelaxation  (update magnetization according to M=M+Mx(MxB))
    logical     :: equi=.false.             !random spin on the unit-sphere within the cone angle (which gets modified)
    logical     :: sphere=.false.           !sphereft sampling (sampling.f90)

    !fluction parameters
    logical     :: do_fluct=.True.                      !calculate fluctuation parameters
    real(8)     :: fluct_dir(3)=[1.0d0,0.0d0,0.0d0]     !direction with respect to which the fluctuation parameters are calculated
contains
    procedure   :: bcast => bcast_MC
    procedure   :: read_file => rw_MC
end type


private

contains

subroutine bcast_MC(this,com)
    use mpi_basic
    class(MC_input),intent(inout)  :: this
    type(mpi_type),intent(in)      :: com
#ifdef CPP_MPI
    integer     ::  ierr
    Call MPI_BCAST( this%n_Tsteps        ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%n_sizerelax     ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%N_relaxation    ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%T_relax         ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%T_auto          ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%Total_MC_Steps  ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%N_Topt          ,1,MPI_INTEGER,com%mas,com%com,ierr)

    Call MPI_BCAST( this%print_relax,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%Cor_log    ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%do_fluct   ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%ising      ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%underrelax ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%overrelax  ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%equi       ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%sphere     ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%expval_save,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%expval_read,1,MPI_LOGICAL,com%mas,com%com,ierr)

    Call MPI_BCAST( this%cone      ,1,MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
    Call MPI_BCAST( this%fluct_dir ,3,MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine rw_MC(inp_MC,io_in)
    use m_io_utils
    use m_io_files_utils
    class(MC_input),intent(out)  ::  inp_MC
    integer,intent(in),optional ::  io_in
    ! internal
    logical :: ising,underrel,overrel,sphere,equi
    logical :: methods(5)
    integer :: io_input,i
    
    methods=.False.
    io_input=open_file_read('input')
    
    call get_parameter(io_input,'input','n_Tsteps',inp_MC%n_Tsteps)
    call get_parameter(io_input,'input','n_sizerelax',inp_MC%n_sizerelax)
    call get_parameter(io_input,'input','n_relaxation',inp_MC%n_relaxation)
    call get_parameter(io_input,'input','Total_MC_Steps',inp_MC%Total_MC_Steps)
    call get_parameter(io_input,'input','T_relax',inp_MC%T_relax)
    call get_parameter(io_input,'input','T_auto',inp_MC%T_auto)
    call get_parameter(io_input,'input','n_Topt',inp_MC%N_Topt)
    
    call get_parameter(io_input,'input','cone',inp_MC%cone)
    call get_parameter(io_input,'input','print_relax',inp_MC%print_relax)
    call get_parameter(io_input,'input','Cor_log',inp_MC%Cor_log)

    call get_parameter(io_input,'input','do_fluct',inp_MC%do_fluct)
    call get_parameter(io_input,'input','fluct_direction',inp_MC%fluct_dir)

    call get_parameter(io_input,'input','ising',inp_MC%ising)
    call get_parameter(io_input,'input','underrelaxation',inp_MC%underrelax)
    call get_parameter(io_input,'input','overrelaxation',inp_MC%overrelax)
    call get_parameter(io_input,'input','sphere_samp',inp_MC%sphere)
    call get_parameter(io_input,'input','equi_samp',inp_MC%equi)

    call get_parameter(io_input,'input','expval_save',inp_MC%expval_save)
    call get_parameter(io_input,'input','expval_read',inp_MC%expval_read)

    methods=[inp_mc%ising,inp_MC%underrelax,inp_MC%overrelax,inp_MC%sphere,inp_MC%equi]
    if(count(methods)>1)then
        write(*,*) "CANNOT SET MORE THAN ONE OF THE FOLLOWING SAMPLING METHODS AT ONCE:"
        write(*,*) "ising:",methods(1)
        write(*,*) "underrelaxation:",methods(2)
        write(*,*) "overrelaxation:",methods(3)
        write(*,*) "sphere_samp:",methods(4)
        write(*,*) "equi_samp:",methods(5)
        STOP
    elseif(count(methods)==0)then
        write(*,*) "NEED TO SET AT LEAST ONE OF THE MONTECARLO SAMPLING METHODS"
        STOP
    endif

    call close_file('input',io_input)
    
end subroutine


end module
