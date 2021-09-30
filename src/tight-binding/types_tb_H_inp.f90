module m_types_tb_h_inp
!module which contains the types which read the manual tight-binding Hamiltonian input
!it definitely makes sense to put the separate terms in separate files, this got out of hand
use m_derived_types, only: lattice
use mpi_basic
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public TB_H_par,TB_hopping, TB_Jsd, TB_delta, alloc_TB_H, TBio_delta_onsite_scf, TBio_defect
type,abstract :: TB_H_par
contains
    procedure(int_print)  ,deferred :: print_std
    procedure(int_is_zero),deferred :: is_zero
    procedure(int_read),deferred    :: local_read
    procedure(int_bcast),deferred   :: bcast
end type

abstract interface
    subroutine int_check(par, lat)
        import TB_H_par, lattice
        class(TB_H_par), intent(inout):: par
        type(lattice),intent(in)      :: lat
    end subroutine

    subroutine int_bcast(this, comm)
        import TB_H_par, mpi_type
        class(TB_H_par),intent(inout)   :: this
        type(mpi_type),intent(in)       :: comm
    end subroutine

    subroutine int_read(par, str, iostat)
        import TB_H_par
        class(TB_H_par), intent(inout):: par
        character(len=*),intent(in)   :: str
        integer, intent(out)          :: iostat
    end subroutine

    subroutine int_print(this,io_in)
        !prints the set information out to the standart output
        import TB_H_par
        class(TB_H_par),intent(in)     :: this
        integer,intent(in),optional    :: io_in
    end subroutine

    function int_is_zero(this)result(is_zero)
        !check if the value is zero
        import TB_H_par
        class(TB_H_par),intent(in)     :: this
        logical :: is_zero
    end function
end interface

type,extends(TB_H_par) :: TB_hopping
    integer     :: attype(2)=[0,0]      !atom types which are connected
    integer     :: orbital(2)=[0,0]     !orbitals of atom-type (excluding spin)
    integer     :: spin(2)=[0,0]        !spin of each orbital  (0=no spin, 1= spin-up, 2= spin-down)    !check after setup that this is sensible
    integer     :: dist=0               !nth distance (1=nearest neighbor)
    real(8)     :: val=0.0d0            !value for this pair of atom-types at distance(dist)
    contains
    procedure :: local_read => TB_hopping_read
    procedure :: TB_hopping_assign
    generic   :: assignment(=) => TB_hopping_assign
    procedure :: print_std => hop_print
    procedure :: is_zero => hopping_is_zero
    procedure :: check => TB_hopping_check
    procedure :: bcast => TB_hopping_bcast
end type

type,extends(TB_H_par) :: TB_delta
    integer     :: attype(2)=[0,0]          !atom types which are connected
    integer     :: orbital(2)=[0,0]         !orbitals of atom-type (excluding spin)
    integer     :: dist=0                   !nth distance (1=nearest neighbor)
    complex(8)  :: val=(0.0d0,0.0d0)   !value for this pair of atom-types at distance(dist)
    contains
    procedure :: local_read => TB_delta_read
    procedure :: TB_delta_assign
    generic :: assignment(=) => TB_delta_assign
    procedure :: print_std => del_print
    procedure :: is_zero => delta_is_zero
    procedure :: check => TB_delta_check
    procedure :: bcast => TB_delta_bcast
end type

type,extends(TB_H_par) :: TB_Jsd
    integer     :: attype=0     !atom type
    integer     :: orbital=0    !orbital 
    real(8)     :: val=0.0d0    !magnitude of Jsd coupling
    contains
    procedure :: TB_Jsd_read
    procedure :: local_read => TB_Jsd_read
    procedure :: TB_Jsd_assign
    generic :: assignment(=) => TB_Jsd_assign
    procedure :: print_std => Jsd_print
    procedure :: is_zero => Jsd_is_zero
    procedure :: check => TB_Jsd_check
    procedure :: bcast => TB_Jsd_bcast
end type

type,extends(TB_H_par) :: TBio_delta_onsite_scf
    integer     :: attype=0     !atom type
    integer     :: orbital(2)=0    !orbital of atom-type (excluding spin)
    real(8)     :: val=0.0d0    !attractive potential magnitude
    !add maximal energy?
    contains
    procedure :: local_read => delta_onsite_scf_read
    procedure :: delta_onsite_scf_assign
    generic :: assignment(=) => delta_onsite_scf_assign
    procedure :: print_std => delta_onsite_scf_print
    procedure :: is_zero => delta_onsite_scf_is_zero
    procedure :: check => delta_onsite_scf_check
    procedure :: bcast => delta_onsite_scf_bcast
end type

type,extends(TB_H_par) :: TBio_defect
    integer     :: atom=0           !atom id
    integer     :: orbital=0        !orbital of atom-type (excluding spin)
    integer     :: site(3)=[0,0,0]  !position in super-lattice
    real(8)     :: nonmag=0.0d0     !spin independent term
    real(8)     :: mag=0.0d0       !attractive potential magnitude
    !add maximal energy?
    contains
    procedure :: local_read => defect_read
    procedure :: defect_assign
    generic :: assignment(=) => defect_assign
    procedure :: print_std => defect_print
    procedure :: is_zero => defect_is_zero
    procedure :: check => defect_check
    procedure :: bcast => defect_bcast
end type

contains 

subroutine alloc_TB_H(par,var_name)
    class(TB_H_par),allocatable,intent(out) ::  par
    character(len=*) var_name

    select case ( var_name)
    case("TB_hopping")
        allocate(TB_hopping::par)
    case("TB_delta")
        allocate(TB_delta::par)
    case("TB_Jsd")
        allocate(TB_Jsd::par)
    case("TB_scfdelta")
        allocate(TBio_delta_onsite_scf::par)
    case("TB_defect")
        allocate(TBio_defect::par)
    case default
        write(error_unit,'(3A)') 'Could not allocate general TB_H_par input type based on name: "',var_name,'"'
        STOP "Implement in types_tb_H_inp.f90, probably new type has been defined"
    end select
end subroutine

subroutine defect_read(par, str, iostat)
    class(TBio_defect), intent(inout)   :: par
    character(len=*),intent(in)         :: str
    integer, intent(out)                :: iostat
    
    type(TBio_defect)             :: tmp
   
    read(str,*,iostat=iostat) tmp%atom,tmp%orbital,tmp%nonmag,tmp%mag, tmp%site
    if (iostat > 0)then    !try to read without z-site
        read(str,*,iostat=iostat) tmp%atom,tmp%orbital,tmp%nonmag,tmp%mag, tmp%site(1:2)
        tmp%site(3)=1
    endif
    if (iostat > 0)then    !try to read without yz-site
        read(str,*,iostat=iostat) tmp%atom,tmp%orbital,tmp%nonmag,tmp%mag, tmp%site(1)
        tmp%site(2:3)=1
    endif
    if(iostat==0) par=tmp
end subroutine

subroutine defect_assign(par,par_in)
    class(TBio_defect), intent(out):: par
    type(TBio_defect) , intent(in ):: par_in

    par%atom   =par_in%atom  
    par%orbital=par_in%orbital 
    par%nonmag =par_in%nonmag 
    par%mag    =par_in%mag       
    par%site   =par_in%site
end subroutine

subroutine defect_check(this,lat)
    class(TBio_defect), intent(in)  :: this
    type(lattice),intent(in)        :: lat
    logical     :: spin_error

    if(this%atom<1.or.this%atom>size(lat%cell%atomic))then
        write(error_unit,'(//A/A)') "Atom in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(this%orbital<1.or.this%orbital>lat%cell%atomic(this%atom)%orbitals)then
        write(error_unit,'(//A/A)') "Orbital in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(any(this%site<1).or.any(this%site>lat%dim_lat))then
        write(error_unit,'(//A/A)') "Site input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif
    !will not catch in Hamiltonian is not magnetic but this%spin=2... but that information is not in lat..., either check later,ignore, or include here
end subroutine

subroutine defect_bcast(this,comm)
    use mpi_basic
    use mpi_util
    class(TBio_defect), intent(inout)   :: this
    type(mpi_type),intent(in)           :: comm
#ifdef CPP_MPI
    integer     :: arr_i(5)
    real(8)     :: arr_r(2)
    arr_i=[this%atom,this%orbital,this%site(1),this%site(2),this%site(3)]
    arr_r=[this%nonmag, this%mag]

    Call bcast(arr_i,comm)
    Call bcast(arr_r,comm)

    this%atom =arr_i(1)
    this%orbital=arr_i(2)
    this%site   =arr_i(3:5)
    this%nonmag =arr_r(1)
    this%mag    =arr_r(2)
#else
    continue
#endif
end subroutine


subroutine TB_hopping_read(par, str, iostat)
    class(TB_hopping), intent(inout)    :: par
    character(len=*),intent(in)         :: str
    integer, intent(out)                :: iostat
    
    type(TB_hopping)    :: tmp
   
    read(str,*,iostat=iostat) tmp%attype,tmp%orbital,tmp%spin,tmp%dist,tmp%val
    if (iostat > 0)then    !try to read without spin
        read(str,*,iostat=iostat) tmp%attype,tmp%orbital,tmp%dist,tmp%val
    endif
    if(iostat==0) par=tmp
    if(par%attype(1)>par%attype(2))then
        par%attype=tmp%attype(2:1:-1)
        par%orbital=tmp%orbital(2:1:-1)
    endif
end subroutine

subroutine TB_hopping_assign(par,par_in)
    class(TB_hopping), intent(out):: par
    type(TB_hopping) , intent(in ):: par_in

    par%attype =par_in%attype  
    par%orbital=par_in%orbital 
    par%spin   =par_in%spin    
    par%dist   =par_in%dist       
    par%val    =par_in%val        
end subroutine

subroutine TB_hopping_check(this,lat)
    class(TB_hopping), intent(in)   :: this
    type(lattice),intent(in)        :: lat
    logical     :: spin_error

    if(any(this%attype<1).or.any(this%attype>lat%cell%n_attype))then
        write(error_unit,'(//A/A)') "Atom types in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(any(this%orbital<1).or.any(this%orbital>lat%cell%atomic(this%attype)%orbitals))then
        write(error_unit,'(//A/A)') "Orbital in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(this%dist<0.or.(this%dist==0.and.this%attype(1)/=this%attype(2)))then
        write(error_unit,'(//A/A)') "Distance in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    spin_error=any(this%spin<0)
    spin_error=spin_error.or.any(this%spin>2)
    spin_error=spin_error.or.any(this%spin==0).and.any(this%spin/=0)
    !will not catch in Hamiltonian is not magnetic but this%spin=2... but that information is not in lat..., either check later,ignore, or include here
    if(spin_error)then
        write(error_unit,'(//A/A)') "Spin in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif
end subroutine

subroutine TB_hopping_bcast(this,comm)
    use mpi_basic
    use mpi_util
    class(TB_hopping), intent(inout)    :: this
    type(mpi_type),intent(in)           :: comm
#ifdef CPP_MPI
    integer     :: arr_i(7)
    arr_i=[this%attype(1),this%attype(2),this%orbital(1),this%orbital(2),this%spin(1),this%spin(2),this%dist]

    Call bcast(arr_i,comm)
    Call bcast(this%val,comm)
    this%attype =arr_i(1:2)
    this%orbital=arr_i(3:4)
    this%spin   =arr_i(5:6)
    this%dist   =arr_i(7)
#else
    continue
#endif
end subroutine


subroutine TB_delta_read(par, str, iostat)
    class(TB_delta), intent(inout)  :: par
    character(len=*),intent(in)     :: str
    integer, intent(out)            :: iostat
    
    type(TB_delta)  :: tmp
    real(8)         :: tmp_val(2)
   
    read(str,*,iostat=iostat) tmp%attype, tmp%orbital ,tmp%dist, tmp%val
    if (iostat > 0)then    !try to read 2 real values format
        read(str,*,iostat=iostat) tmp%attype, tmp%orbital ,tmp%dist, tmp_val
        tmp%val=cmplx(tmp_val(1),tmp_val(2),8)
    endif
    if(iostat==0) par=tmp
    if(par%attype(1)>par%attype(2))then
        par%attype=tmp%attype(2:1:-1)
        par%orbital=tmp%orbital(2:1:-1)
    endif
end subroutine


subroutine TB_delta_assign(par,par_in)
    class(TB_delta), intent(out):: par
    type(TB_delta),  intent(in ):: par_in

    par%attype =par_in%attype  
    par%orbital=par_in%orbital 
    par%dist   =par_in%dist       
    par%val    =par_in%val        
end subroutine

subroutine TB_delta_check(this,lat)
    class(TB_delta), intent(in) :: this
    type(lattice),intent(in)    :: lat
    logical     :: spin_error

    if(any(this%attype<1).or.any(this%attype>lat%cell%n_attype))then
        write(error_unit,'(//A/A)') "Atom types in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(any(this%orbital<1).or.any(this%orbital>lat%cell%atomic(this%attype)%orbitals))then
        write(error_unit,'(//A/A)') "Orbital in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(this%dist<-1.or.(this%dist==0.and.this%attype(1)/=this%attype(2)))then
        write(error_unit,'(//A/A)') "Distance in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif
end subroutine

subroutine TB_delta_bcast(this,comm)
    use mpi_basic
    use mpi_util
    class(TB_delta), intent(inout)  :: this
    type(mpi_type),intent(in)       :: comm
#ifdef CPP_MPI
    integer     :: arr_i(5)
    arr_i=[this%attype(1),this%attype(2),this%orbital(1),this%orbital(2),this%dist]

    Call bcast(arr_i,comm)
    Call bcast(this%val,comm)
    this%attype =arr_i(1:2)
    this%orbital=arr_i(3:4)
    this%dist   =arr_i(5)
#else
    continue
#endif
end subroutine


subroutine TB_Jsd_read(par, str, iostat)
    class(TB_Jsd), intent(inout)  :: par
    character(len=*),intent(in)   :: str
    integer, intent(out)          :: iostat
    
    type(TB_Jsd)                  :: tmp
   
    read(str,*,iostat=iostat) tmp%attype, tmp%orbital , tmp%val
    if(iostat==0) par=tmp
end subroutine


subroutine TB_Jsd_assign(par,par_in)
    class(TB_Jsd), intent(out):: par
    type(TB_Jsd),  intent(in ):: par_in

    par%attype =par_in%attype  
    par%orbital=par_in%orbital 
    par%val    =par_in%val        
end subroutine


subroutine delta_onsite_scf_check(this,lat)
    class(TBio_delta_onsite_scf), intent(in)    :: this
    type(lattice),intent(in)                    :: lat

    if(this%attype<1.or.this%attype>lat%cell%n_attype)then
        write(error_unit,'(//A/A)') "Atom types in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(any(this%orbital<1).or.any(this%orbital>lat%cell%atomic(this%attype)%orbitals))then
        write(error_unit,'(//A/A)') "Orbital in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif
    !check negative value?
end subroutine

subroutine delta_onsite_scf_read(par, str, iostat)
    class(TBio_delta_onsite_scf), intent(inout) :: par
    character(len=*),intent(in)                 :: str
    integer, intent(out)                        :: iostat
    
    type(TBio_delta_onsite_scf)        :: tmp
    read(str,*,iostat=iostat) tmp%attype, tmp%orbital, tmp%val
    if(iostat==0) par=tmp
end subroutine


subroutine delta_onsite_scf_assign(par,par_in)
    class(TBio_delta_onsite_scf), intent(out):: par
    type(TBio_delta_onsite_scf),  intent(in ):: par_in

    par%attype =par_in%attype  
    par%orbital=par_in%orbital 
    par%val    =par_in%val        
end subroutine

subroutine delta_onsite_scf_bcast(this,comm)
    use mpi_basic
    use mpi_util
    class(TBio_delta_onsite_scf), intent(inout) :: this
    type(mpi_type),intent(in)                   :: comm
#ifdef CPP_MPI
    integer     :: arr_i(3)
    arr_i=[this%attype,this%orbital(1),this%orbital(2)]

    Call bcast(arr_i,comm)
    Call bcast(this%val,comm)
    this%attype =arr_i(1)
    this%orbital=arr_i(2:3)
#else
    continue
#endif
end subroutine



subroutine TB_Jsd_check(this,lat)
    class(TB_Jsd), intent(in) :: this
    type(lattice),intent(in)    :: lat
    logical     :: spin_error

    if(this%attype<1.or.this%attype>lat%cell%n_attype)then
        write(error_unit,'(//A/A)') "Atom types in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif

    if(this%orbital<1.or.this%orbital>lat%cell%atomic(this%attype)%orbitals)then
        write(error_unit,'(//A/A)') "Orbital in TB-input out of range","Error in entry:" 
        Call this%print_std(error_unit)
        Stop "Check input"
    endif
end subroutine

subroutine TB_Jsd_bcast(this,comm)
    use mpi_basic
    use mpi_util
    class(TB_Jsd), intent(inout)    :: this
    type(mpi_type),intent(in)       :: comm
#ifdef CPP_MPI

    Call bcast(this%attype,comm)
    Call bcast(this%orbital,comm)
    Call bcast(this%val,comm)
#else
    continue
#endif
end subroutine


subroutine del_print(this,io_in)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(TB_delta),intent(in)    :: this
    integer,intent(in),optional     :: io_in
    integer                         :: io_unit

    io_unit=output_unit
    if(present(io_in)) io_unit=io_in

    write(io_unit,'(A)')         'Delta Hamiltonian input data:'
    write(io_unit,'(A,2I6)')     '  atom types:', this%attype
    write(io_unit,'(A,2I6)')     '  orbitals  :', this%orbital
    write(io_unit,'(A,2I6)')     '  distance  :', this%dist
    write(io_unit,'(A,2E16.8/)') '  energy    :', this%val
end subroutine

subroutine delta_onsite_scf_print(this,io_in)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(TBio_delta_onsite_scf),intent(in)     :: this
    integer,intent(in),optional                 :: io_in
    integer                                     :: io_unit

    io_unit=output_unit
    if(present(io_in)) io_unit=io_in

    write(io_unit,'(A)')         'onsite scf delta Hamiltonian input data:'
    write(io_unit,'(A,2I6)')     '  atom type:', this%attype
    write(io_unit,'(A,2I6)')     '  orbital  :', this%orbital
    write(io_unit,'(A,2E16.8/)') '  potential:', this%val
end subroutine


subroutine hop_print(this,io_in)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(TB_hopping),intent(in)    :: this
    integer,intent(in),optional     :: io_in
    integer                         :: io_unit

    io_unit=output_unit
    if(present(io_in)) io_unit=io_in

    write(io_unit,'(A)')         'Hopping Hamiltonian input data:'
    write(io_unit,'(A,2I6)')     '  atom types:', this%attype
    write(io_unit,'(A,2I6)')     '  orbitals  :', this%orbital
    write(io_unit,'(A,2I6)')     '  spin      :', this%spin
    write(io_unit,'(A,2I6)')     '  distance  :', this%dist
    write(io_unit,'(A,E16.8/)')  '  energy    :', this%val
end subroutine

subroutine Jsd_print(this,io_in)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(TB_Jsd),intent(in)        :: this
    integer,intent(in),optional     :: io_in
    integer                         :: io_unit

    io_unit=output_unit
    if(present(io_in)) io_unit=io_in

    write(io_unit,'(A)')         'Jsd Hamiltonian input data:'
    write(io_unit,'(A,2I6)')     '  atom types:', this%attype
    write(io_unit,'(A,2I6)')     '  orbitals  :', this%orbital
    write(io_unit,'(A,E16.8/)')  '  energy    :', this%val
end subroutine

subroutine defect_print(this,io_in)
    use, intrinsic :: iso_fortran_env, only : output_unit
    class(TBio_defect),intent(in)   :: this
    integer,intent(in),optional     :: io_in
    integer                         :: io_unit

    io_unit=output_unit
    if(present(io_in)) io_unit=io_in

    write(io_unit,'(A)')         'defect Hamiltonian input data:'
    write(io_unit,'(A,2I6)')     '  atom id     :', this%atom
    write(io_unit,'(A,2I6)')     '  orbital     :', this%orbital
    write(io_unit,'(A,E16.8)')   '  non-magnetic:', this%nonmag
    write(io_unit,'(A,E16.8)')   '  magnetic    :', this%mag
    write(io_unit,'(A,3I6/)')    '  site        :', this%site
end subroutine


function Jsd_is_zero(this)result(is_zero)
    class(TB_Jsd),intent(in)   ::  this
    logical :: is_zero

    is_zero=this%val==(0.0d0,0.0d0)
end function

function hopping_is_zero(this)result(is_zero)
    class(TB_hopping),intent(in)   ::  this
    logical :: is_zero

    is_zero=this%val==0.0d0
end function

function delta_is_zero(this)result(is_zero)
    class(TB_delta),intent(in)   ::  this
    logical :: is_zero

    is_zero=this%val==0.0d0
end function

function delta_onsite_scf_is_zero(this)result(is_zero)
    class(TBio_delta_onsite_scf),intent(in)   ::  this
    logical :: is_zero

    is_zero=this%val==0.0d0
end function


function defect_is_zero(this)result(is_zero)
    class(TBio_defect),intent(in)   ::  this
    logical :: is_zero

    is_zero=this%mag==0.0d0.and.this%nonmag==0.0d0
end function

end module
