module m_types_tb_h_inp
use m_derived_types, only: lattice
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
implicit none
private
public TB_H_par,TB_hopping, TB_Jsd, TB_delta, alloc_TB_H, TBio_delta_onsite_scf, TBio_defect
type,abstract :: TB_H_par
contains
    procedure(int_print)  ,deferred :: print_std
    procedure(int_is_zero),deferred :: is_zero
    procedure(int_read),deferred    :: local_read
    generic :: read(formatted) => local_read
end type

abstract interface
    subroutine int_check(par, lat)
        import TB_H_par, lattice
        class(TB_H_par), intent(inout):: par
        type(lattice),intent(in)      :: lat
    end subroutine

    subroutine int_read(par, unit, iotype, v_list, iostat, iomsg)
        import TB_H_par
        class(TB_H_par), intent(inout):: par
        integer, intent(in)           :: unit
        character(*), intent(in)      :: iotype
        integer, intent(in)           :: v_list(:)
        integer, intent(out)          :: iostat
        character(*), intent(inout)   :: iomsg
    end subroutine

    subroutine int_print(this,io_in)
        !prints the set information out to the standart output
        import TB_H_par
        class(TB_H_par),intent(in)     :: this
        integer,intent(in),optional     :: io_in
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

subroutine defect_read(par, unit, iotype, v_list, iostat, iomsg)
    class(TBio_defect), intent(inout):: par
    integer, intent(in)           :: unit
    character(*), intent(in)      :: iotype
    integer, intent(in)           :: v_list(:)
    integer, intent(out)          :: iostat
    character(*), intent(inout)   :: iomsg
    
    type(TBio_defect)             :: tmp
   
    read(unit,*,iostat=iostat,iomsg=iomsg) tmp%atom,tmp%orbital,tmp%nonmag,tmp%mag, tmp%site
    backspace(unit)
    if (iostat > 0)then    !try to read without z-site
        read(unit,*,iostat=iostat,iomsg=iomsg) tmp%atom,tmp%orbital,tmp%nonmag,tmp%mag, tmp%site(1:2)
        tmp%site(3)=1
        backspace(unit)
    endif
    if (iostat > 0)then    !try to read without yz-site
        read(unit,*,iostat=iostat,iomsg=iomsg) tmp%atom,tmp%orbital,tmp%nonmag,tmp%mag, tmp%site(1)
        tmp%site(2:3)=1
        backspace(unit)
    endif
    read(unit,*)    !together with backspace make sure that is advances
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


subroutine TB_hopping_read(par, unit, iotype, v_list, iostat, iomsg)
    class(TB_hopping), intent(inout):: par
    integer, intent(in)           :: unit
    character(*), intent(in)      :: iotype
    integer, intent(in)           :: v_list(:)
    integer, intent(out)          :: iostat
    character(*), intent(inout)   :: iomsg
    
    type(TB_hopping)              :: tmp
   
    read(unit,*,iostat=iostat,iomsg=iomsg) tmp%attype,tmp%orbital,tmp%spin,tmp%dist,tmp%val
    backspace(unit)
    if (iostat > 0)then    !try to read without spin
        read(unit,*,iostat=iostat,iomsg=iomsg) tmp%attype,tmp%orbital,tmp%dist,tmp%val
        backspace(unit)
    endif
    read(unit,*)    !together with backspace make sure that is advances
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

subroutine TB_delta_read(par, unit, iotype, v_list, iostat, iomsg)
    class(TB_delta), intent(inout):: par
    integer, intent(in)           :: unit
    character(*), intent(in)      :: iotype
    integer, intent(in)           :: v_list(:)
    integer, intent(out)          :: iostat
    character(*), intent(inout)   :: iomsg
    
    type(TB_delta)                :: tmp
   
    read(unit,*,iostat=iostat,iomsg=iomsg) tmp%attype, tmp%orbital ,tmp%dist, tmp%val
    backspace(unit)
    if (iostat > 0)then    !try to read 2 real values format
        read(unit,*,iostat=iostat,iomsg=iomsg) tmp%attype, tmp%orbital ,tmp%dist, tmp%val%re, tmp%val%im
        backspace(unit)
    endif
    read(unit,*)    !together with backspace make sure that it advances
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

subroutine TB_Jsd_read(par, unit, iotype, v_list, iostat, iomsg)
    class(TB_Jsd), intent(inout)  :: par
    integer, intent(in)           :: unit
    character(*), intent(in)      :: iotype
    integer, intent(in)           :: v_list(:)
    integer, intent(out)          :: iostat
    character(*), intent(inout)   :: iomsg
    
    type(TB_Jsd)                  :: tmp
   
    read(unit,*,iostat=iostat,iomsg=iomsg) tmp%attype, tmp%orbital , tmp%val
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

subroutine delta_onsite_scf_read(par, unit, iotype, v_list, iostat, iomsg)
    class(TBio_delta_onsite_scf), intent(inout)  :: par
    integer, intent(in)           :: unit
    character(*), intent(in)      :: iotype
    integer, intent(in)           :: v_list(:)
    integer, intent(out)          :: iostat
    character(*), intent(inout)   :: iomsg
    
    type(TBio_delta_onsite_scf)        :: tmp
   
    read(unit,*,iostat=iostat,iomsg=iomsg) tmp%attype, tmp%orbital , tmp%val
    if(iostat==0) par=tmp
end subroutine


subroutine delta_onsite_scf_assign(par,par_in)
    class(TBio_delta_onsite_scf), intent(out):: par
    type(TBio_delta_onsite_scf),  intent(in ):: par_in

    par%attype =par_in%attype  
    par%orbital=par_in%orbital 
    par%val    =par_in%val        
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
