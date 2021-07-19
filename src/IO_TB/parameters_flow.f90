module m_parameters_TB_IO_FLOW
!module of type which contains the IO-parameters that determine which parts of the tight binding calculation are done
implicit none
private
public parameters_TB_IO_FLOW
integer     :: Npar=14

type parameters_TB_IO_FLOW
    logical         ::  do_r=.False.    !do any calculations in real-space
    logical         ::  dos_r=.False.   !calculated density of states in real-space.
    logical         ::  spec_r=.False.  !writes spectrum, i.e. all energy eigenvalues
    logical         ::  fermi_r=.False. !calculates fermi energy in real-space ( rather obsolete)
    logical         ::  occ_r=.False.   !print the occupation in real-space and derivative at a given energy
    logical         ::  occ_mult_r=.False. !print the occupation in real-space and derivative at several energies

    logical         ::  read_solution_r=.True.  !read the solution of the real-space diagonalization (at least try)
    logical         ::  write_solution_r=.False.!save the solution of the real-space diagonalization 

    logical         ::  do_k=.False.        !do any calcuation in reciprocal space
    logical         ::  dos_k=.False.       !do density of state in k-space
    logical         ::  highs_k=.False.     !calculate bandstructure along high-symmetry lines
    logical         ::  fermi_k=.False.     !calculate fermi energy (relatively obsolete)
    logical         ::  fermi_dos_k=.False. !plot fermi surface at Fermi energy
    logical         ::  proj_energy=.false. !do the projection onto each orbital summed over all occupied state (T=0)  (used for terrible Jsd estimation) 
contains
    procedure   :: read_file
    procedure   :: bcast => bcast_local
end type

contains


subroutine read_file(this,io,fname)
    use m_io_utils
    class(parameters_TB_IO_flow),intent(out)    :: this
    integer,intent(in)                          :: io
    character(len=*), intent(in)                :: fname

    call get_parameter(io,fname,'do_TB_r',      this%do_r)
    call get_parameter(io,fname,'do_dos_r',     this%dos_r)
    call get_parameter(io,fname,'do_occ_r',     this%occ_r)
    call get_parameter(io,fname,'do_spec_r',    this%spec_r)
    call get_parameter(io,fname,'do_fermi_r',   this%fermi_r)
    call get_parameter(io,fname,'do_occ_mult_r',this%occ_mult_r)

    call get_parameter(io,fname,'TB_read_solution_r',   this%read_solution_r)
    call get_parameter(io,fname,'TB_write_solution_r',  this%write_solution_r)

    call get_parameter(io,fname,'do_TB_k',          this%do_k)
    call get_parameter(io,fname,'do_dos_k',         this%dos_k)
    call get_parameter(io,fname,'do_fermi_k',       this%fermi_k)
    call get_parameter(io,fname,'do_fermi_dos_k',   this%fermi_dos_k)
    call get_parameter(io,fname,'do_highs_k',       this%highs_k)
    call get_parameter(io,fname,'do_proj_energy.',  this%proj_energy)
end subroutine

subroutine bcast_local(this,comm)
    use mpi_basic
    use mpi_util
    class(parameters_TB_IO_FLOW),intent(inout)  ::  this
    type(mpi_type),intent(in)                   ::  comm
#ifdef CPP_MPI
    logical     :: arr(Npar)
    Call to_array(this,arr)
    Call bcast(arr,comm)
    Call from_array(this,arr)
#else
    continue
#endif
end subroutine 


subroutine to_array(this,arr)
    class(parameters_TB_IO_flow),intent(in) :: this
    logical,intent(out)                     :: arr(Npar)

    arr =[&
    this%do_r, &
    this%dos_r, &
    this%spec_r, &
    this%fermi_r,& 
    this%occ_r, &
    this%occ_mult_r, &
    this%read_solution_r, &
    this%write_solution_r, &
    this%do_k, &
    this%dos_k, &
    this%highs_k, & 
    this%fermi_k, &
    this%fermi_dos_k, &
    this%proj_energy ]
end subroutine

subroutine from_array(this,arr)
    class(parameters_TB_IO_flow),intent(inout)  :: this
    logical,intent(in)                          :: arr(Npar)

    this%do_r             =arr( 1)  
    this%dos_r            =arr( 2)  
    this%spec_r           =arr( 3)  
    this%fermi_r          =arr( 4)  
    this%occ_r            =arr( 5)  
    this%occ_mult_r       =arr( 6)  
    this%read_solution_r  =arr( 7)  
    this%write_solution_r =arr( 8)  
    this%do_k             =arr( 9)  
    this%dos_k            =arr(10)  
    this%highs_k          =arr(11)  
    this%fermi_k          =arr(12)  
    this%fermi_dos_k      =arr(13)  
    this%proj_energy      =arr(14)  
end subroutine


end module
