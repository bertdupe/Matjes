module m_kgrid_int
use m_kgrid, only: k_grid_t, kmesh_t
use m_Hk, only: Hk_inp_t,Hk_eval
use m_tb_types
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
private
public  :: get_kmesh, kmesh_t, k_mesh_int

type,extends(kmesh_t) :: k_mesh_int
    real(8),allocatable :: k(:,:)
    integer             :: Nk=0
    integer             :: Nk_init=0
contains
    procedure :: get_k 
    procedure :: get_Nk
    procedure :: get_normalize

    procedure :: set
    procedure :: kmesh_write
    procedure :: kmesh_read
    generic :: write(formatted) => kmesh_write
    generic :: read(formatted) => kmesh_read
end type

contains

subroutine get_kmesh(kgrid,lat,grid,fname)
    use m_derived_types, only: lattice
    class(kmesh_t),intent(inout),allocatable    ::  kgrid
    type(lattice), intent(in)                   ::  lat
    integer,intent(in)                          ::  grid(3)
    character(:),intent(in),allocatable         ::  fname
    logical         ::  fexist
    integer         ::  io,stat

    fexist=.false.
    if(allocated(fname))then
        inquire(file=fname,exist=fexist)
        if(.not.fexist)then
            write(error_unit,'(//3A)')  "Warning, input suggests to read kmesh from file '",fname,"', but that file does not exist"
            write(error_unit,'(A,3I6)') "         using instead grid with size :",grid
        endif
    endif
    if(fexist)then
        allocate(k_mesh_int::kgrid)
    else
        allocate(k_grid_t::kgrid)
    endif
    select type(kgrid)
    type is(k_mesh_int)
        write(output_unit,'(/2A)') "Reading kmesh from file: ", fname
        open(newunit=io,file=fname,action='read',status='old')
        read(io,*,iostat=stat) kgrid
        close(io)
        if(stat/=0)then
            write(error_unit,'(3/3A)') "Failed to read kmesh from file: '",fname,"'"
            STOP
        endif
    type is (k_grid_t)
        write(output_unit,'(/A3I6)') "Initializing kgrid with following repetitions:", grid
        Call kgrid%set(lat%a_sc_inv,grid)
    class default
        ERROR STOP "THIS SHOULD NEVER BE REACHED"
    end select
end subroutine
        

    



function get_k(this,i)result(k)
    use m_constants,    only : pi
    class(k_mesh_int),intent(in)    :: this
    integer,intent(in)              :: i
    real(8)                         :: k(3)
    integer                         :: kint(3)

    if(i<1.or.i>this%NK) ERROR STOP "wanted k-grid index not within kgrid-lattice"
    k=this%k(:,i)
end function

function get_Nk(this)result(Nk)
    use m_constants,    only : pi
    class(k_mesh_int),intent(in)    :: this
    integer                         :: Nk

    Nk=this%Nk
end function

function get_normalize(this)result(norm)
    class(k_mesh_int),intent(in)    :: this
    real(8)                         :: norm
    norm=1.0d0/real(this%Nk)
end function

subroutine set(this,kgrid,Hk_inp,h_io,Ecut)
    class(k_mesh_int)            :: this
    type(k_grid_t),intent(in)   :: kgrid
    type(Hk_inp_t),intent(in)   :: Hk_inp
    type(parameters_TB_IO_H),intent(in)     :: h_io
    real(8),intent(in)          :: Ecut(2)
    integer                     :: Nk_init

    real(8),allocatable         :: k_found(:,:)
    integer                     :: Nfound
    real(8)                     :: k(3)
    real(8),allocatable         :: eval(:)
    integer                     :: iE, ik

    if(Ecut(2)<=Ecut(1)) then
        write(error_unit,'(/2(/A,E16.8))') "Ecut(1)=",Ecut(1),"Ecut(2)=",Ecut(2)
        ERROR STOP "integration kmesh energy bound Ecut(2) must be larger than Ecut(1)"
    endif
    Nk_init=kgrid%get_Nk()
    if(Nk_init<1)then
        ERROR STOP "kgrid Nk<1"
    endif

    allocate(k_found(3,Nk_init),source=0.d0)
    Nfound=0
    do ik=1,Nk_init
        k=kgrid%get_K(ik)
        Call Hk_eval(Hk_inp,k,h_io,eval) 
        if(any(eval>Ecut(1).and.eval<Ecut(2)))then
            Nfound=Nfound+1
            k_found(:,Nfound)=k
        endif
        deallocate(eval)
    enddo
    this%Nk=Nfound
    this%Nk_init=Nk_init
    allocate(this%k,source=k_found(:,1:Nfound))
    deallocate(k_found)
end subroutine

subroutine kmesh_write(par, unit, iotype, v_list, iostat, iomsg)
    class(k_mesh_int), intent(in) :: par
    integer, intent(in)           :: unit
    character(*), intent(in)      :: iotype
    integer, intent(in)           :: v_list(:)
    integer, intent(out)          :: iostat
    character(*), intent(inout)   :: iomsg
   
    if(par%Nk<1)then
        write(error_unit,*) "cannot write kmesh since it has no k-points"
        iostat=1
    endif
    write(unit,'(2I12/)',iostat=iostat,iomsg=iomsg) par%Nk,par%Nk_init
    if(iostat/=0) return
    write(unit,'(3E16.8)',iostat=iostat,iomsg=iomsg) par%K
end subroutine

subroutine kmesh_read(par, unit, iotype, v_list, iostat, iomsg)
    class(k_mesh_int), intent(inout) :: par
    integer, intent(in)           :: unit
    character(*), intent(in)      :: iotype
    integer, intent(in)           :: v_list(:)
    integer, intent(out)          :: iostat
    character(*), intent(inout)   :: iomsg

    character(100)                  :: derp
   
    if(allocated(par%K)) deallocate(par%K)
    read(unit,'(A)',iostat=iostat,iomsg=iomsg) derp !no clue why this helps, but otherwise the reading fails
    backspace(unit)
    read(unit,'(2I12/)',iostat=iostat,iomsg=iomsg) par%Nk,par%Nk_init
    allocate(par%K(3,par%Nk),source=0.0d0)
    read(unit,'(3E16.8)',iostat=iostat,iomsg=iomsg) par%K
end subroutine

end module
