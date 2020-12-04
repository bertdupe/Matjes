module m_MC_io
use m_constants, only : pi
use mpi_basic,only: mpi_type
implicit none

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
    logical     :: ising=.false.
    logical     :: underrelax=.false.
    logical     :: overrelax=.false.
    logical     :: equi=.false.
    logical     :: sphere=.false.
    real(8)     :: cone=pi
end type


private
public :: rw_MC,MC_input,bcast_MC

contains

subroutine bcast_MC(this,com)
    use mpi_basic
    class(MC_input),intent(inout)  :: this
    type(mpi_type),intent(in)      :: com
#ifdef CPP_MPI
    integer     ::  ierr
    Call MPI_BCAST( this%n_Tsteps        ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%n_sizerelax     ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%n_thousand      ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%restart_MC_steps,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%T_relax         ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%T_auto          ,1,MPI_INTEGER,com%mas,com%com,ierr)
    Call MPI_BCAST( this%Total_MC_Steps  ,1,MPI_INTEGER,com%mas,com%com,ierr)

    Call MPI_BCAST( this%i_restart  ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%print_relax,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%Cor_log    ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%do_fluct   ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%ising      ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%underrelax ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%overrelax  ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%equi       ,1,MPI_LOGICAL,com%mas,com%com,ierr)
    Call MPI_BCAST( this%sphere     ,1,MPI_LOGICAL,com%mas,com%com,ierr)

    Call MPI_BCAST( this%cone  ,1,MPI_DOUBLE_PRECISION,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

!subroutine rw_MC(n_Tsteps,n_sizerelax,n_thousand,restart_MC_steps,Total_MC_Steps,T_relax,T_auto,cone,i_restart,ising,underrel,overrel,sphere,equi,print_relax,Cor_log)
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
    call get_parameter(io_input,'input','n_relaxation',inp_MC%n_thousand)
    call get_parameter(io_input,'input','restart_MC_steps',inp_MC%restart_MC_steps)
    call get_parameter(io_input,'input','Total_MC_Steps',inp_MC%Total_MC_Steps)
    call get_parameter(io_input,'input','T_relax',inp_MC%T_relax)
    call get_parameter(io_input,'input','T_auto',inp_MC%T_auto)
    
    call get_parameter(io_input,'input','cone',inp_MC%cone)
    call get_parameter(io_input,'input','print_relax',inp_MC%print_relax)
    call get_parameter(io_input,'input','Cor_log',inp_MC%Cor_log)
    call get_parameter(io_input,'input','restart',inp_MC%i_restart)
    call get_parameter(io_input,'input','do_fluct',inp_MC%do_fluct)

    call get_parameter(io_input,'input','ising',inp_MC%ising)
    call get_parameter(io_input,'input','underrelaxation',inp_MC%underrelax)
    call get_parameter(io_input,'input','overrelaxation',inp_MC%overrelax)
    call get_parameter(io_input,'input','sphere_samp',inp_MC%sphere)
    call get_parameter(io_input,'input','equi_samp',inp_MC%equi)

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
