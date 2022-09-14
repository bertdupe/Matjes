module m_input_H_types
implicit none
public


type :: Hr_pair_tensor !tensor hamiltonian real pair
    integer     ::  attype(2)=[0,0] !atom types which are connected
    integer,allocatable     ::  dist(:)    !nth distance (1=nearest neighbor)
    real(8),allocatable     ::  val(:,:)     !value for this tensor pair of atom-types at distance(dist)
    real(8),allocatable     ::  bound(:,:)    !bound along which the tensor should apply
contains
    procedure   ::  prt=>prt_Hr_pair_tensor
end type

type :: Hr_pair !hamiltonian real pair
    integer     ::  attype(2)=[0,0] !atom types which are connected
    integer,allocatable     ::  dist(:)    !nth distance (1=nearest neighbor)
    real(8),allocatable     ::  val(:)     !value for this pair of atom-types at distance(dist)
contains 
    procedure   ::  prt=>prt_Hr_pair
end type

type :: Hr_pair_single !hamiltonian real pair
    integer     ::  attype(2)=[0,0] !atom types which are connected
    integer     ::  dist=0      !nth distance (1=nearest neighbor)
    real(8)     ::  val=0.0     !value for this pair of atom-types at distance(dist)
end type

type :: Hr_pair_single_tensor !hamiltonian tensor real pair
    integer     ::  attype(2)=[0,0] !atom types which are connected
    integer     ::  dist=0      !nth distance (1=nearest neighbor)
    real(8)     ::  val(9)=0.0     !value of the tensor for this pair of atom-types at distance(dist)
    real(8)     ::  bound(3)=0.0    !bound along which the tensor should apply
end type

type :: Hr_triple_single !hamiltonian real triple
    integer     ::  attype(3)=[0,0,0] !atom types which are connected
    integer     ::  dist=0      !nth distance (1=nearest neighbor)
    real(8)     ::  val=0.0     !value for this pair of atom-types at distance(dist)
end type

type :: Hr_triple !hamiltonian real pair
    integer     ::  attype(3)=[0,0,0] !atom types which are connected
    integer,allocatable     ::  dist(:)    !nth distance (1=nearest neighbor)
    real(8),allocatable     ::  val(:)     !value for this pair of atom-types at distance(dist)
contains 
    procedure   ::  prt=>prt_Hr_triple
end type

type :: Hr_quadruplet_single !hamiltonian real triple
    integer     ::  attype(4)=[0,0,0,0] !atom types which are connected
    integer     ::  dist=0      !nth distance (1=nearest neighbor)
    real(8)     ::  val=0.0     !value for this pair of atom-types at distance(dist)
end type

type :: Hr_quadruplet !hamiltonian real pair
    integer     ::  attype(4)=[0,0,0,0] !atom types which are connected
    integer,allocatable     ::  dist(:)    !nth distance (1=nearest neighbor)
    real(8),allocatable     ::  val(:)     !value for this pair of atom-types at distance(dist)
contains
!    procedure   ::  prt=>prt_Hr_quadruplet
end type

type :: io_H_base
    logical :: is_set=.false.   !gets set to true if sensible input is read
    logical :: fft=.false.      !whether fft is used if available
end type

type,extends(io_H_base) :: io_H_aniso
    !anisotropy in space of Cartesian coordinates
    integer,allocatable     ::  attype(:)       !integer of atom type
    real(8),allocatable     ::  val(:,:)        !(4,size(attype))  (1:3) direction, (4) magnitude
    !anisotropy in space of normalized real-space lattice
    integer,allocatable     ::  attype_lat(:)   !integer of atom type
    real(8),allocatable     ::  val_lat(:,:)    !(4,size(attype_lat))  (1:3) direction, (4) magnitude
    real(8)                 :: c_H_ani=1.0d0
end type

type,extends(io_H_base) :: io_H_zeeman
    real(8)             :: c_zeeman=-1.0d0 !constant factor to furthermore rescale zeeman energy
end type

type,extends(io_H_base) :: io_H_J
    type(Hr_pair),allocatable   :: pair(:)
    real(8)                     :: c_H_J=-1.0d0
end type

type,extends(io_H_base) :: io_H_sp3
    type(Hr_triple),allocatable   :: triplet(:)
    integer,allocatable           :: at_type(:)
    real(8)                       :: c_H_sp3=-1.0d0
end type

type,extends(io_H_base) :: io_H_sp4
    integer,allocatable     :: at_type(:)
    real(8),allocatable     :: val(:)
    real(8)                 :: c_H_sp4=-1.0d0
end type

type,extends(io_H_base) :: io_H_D
    type(Hr_triple),allocatable :: trip(:) 
    real(8)                     :: c_H_D=-1.0d0
end type

type,extends(io_H_base) :: io_H_TJ
    type(Hr_pair),allocatable   :: pair(:) 
end type

type,extends(io_H_base) :: io_H_ME_J
    type(Hr_pair),allocatable   :: pair(:) 
end type

type,extends(io_H_base) :: io_H_ME_D
    type(Hr_pair),allocatable   :: pair(:) 
end type

type,extends(io_H_base) :: io_H_Ph
    type(Hr_pair),allocatable   :: pair(:) 
    real(8)                     :: c_ph=-1.0d0
end type

type,extends(io_H_base) :: io_H_stark
    real(8)             :: c_stark=-1.0d0 !constant factor to furthermore rescale stark energy
end type

type,extends(io_H_base) :: io_U_ASR
    type(Hr_pair),allocatable          :: pair(:)
    type(Hr_pair_tensor),allocatable   :: pair_tensor(:)
    integer,allocatable                :: attype(:)
    real(8)                            :: c_ASR=1.0d0 ! to enforce the acoustic sum rule for the phonon energy
end type

type,extends(io_H_base) :: io_H_dipole
    integer     :: period_cutoff(3)=[1,1,1]  !how many periodic images are considered in each direction
!    integer     :: dist_cutoff=-1.0 !not really used yet  
end type

type,extends(io_H_base) :: io_H_Exchten
    type(Hr_pair_tensor),allocatable   :: pair(:)
    real(8)                            :: c_H_Exchten=-1.0d0
end type

type,extends(io_H_base) :: io_H_Mag_Biq
    type(Hr_pair),allocatable   :: pair(:)
    real(8)                     :: c_H_Mbiq=-1.0d0
end type

type,extends(io_H_base) :: io_H_SC_D
    type(Hr_triple),allocatable   :: triplet(:)
    integer,allocatable           :: at_type(:)
    real(8)                       :: c_SC=-1.0d0
end type

type,extends(io_H_base) :: io_H_Ph4
    integer,allocatable     :: at_type(:)
    real(8),allocatable     :: val(:)
end type

type,extends(io_H_base) :: io_H_Ph_Biq
    type(Hr_pair),allocatable   :: pair(:)
end type

type,extends(io_H_base) :: io_H_Force_tensor
    type(Hr_pair_tensor),allocatable   :: pair(:)
    real(8)                            :: c_phtensor=-1.0d0
end type

type :: io_H
    type(io_H_aniso)            :: aniso
    type(io_H_zeeman)           :: zeeman
    type(io_H_ME_J)             :: ME_J
    type(io_H_ME_D)             :: ME_D
    type(io_H_J)                :: J
    type(io_H_sp4)              :: sp4
    type(io_H_sp3)              :: sp3
    type(io_H_TJ)               :: TJ
    type(io_H_D)                :: D
    type(io_H_Ph)               :: F
    type(io_H_stark)            :: stark
    type(io_U_ASR)              :: ASR_Ph
    type(io_H_Mag_Biq)          :: M_biq
    type(io_H_dipole)           :: dip
    type(io_H_Exchten)          :: Exchten
    type(io_H_SC_D)             :: SC
    type(io_H_Ph4)              :: Ph4
    type(io_H_Ph_Biq)           :: U_biq
    type(io_H_Force_tensor)     :: U_foten
    type(io_H_dipole)           :: dip_ph
end type

contains


subroutine reduce_Hr_triple(Hr_triple_in,Hr_triple_out)
    !subroutine which combines the entries with same atom types from Hr_triple_in into Hr_triple_out
    type(Hr_triple_single),intent(in)         ::  Hr_triple_in(:)
    type(Hr_triple),intent(out),allocatable   ::  Hr_triple_out(:)

    integer         ::  at_triples       (3,size(Hr_triple_in))
    integer         ::  at_triples_unique(3,size(Hr_triple_in))
    integer         ::  i,j,ii
    integer         ::  N_triple(size(hr_triple_in))
    integer         ::  entry_unique(size(hr_triple_in))

    if(size(Hr_triple_in)==1)then 
        !do nothing if only one entry is present
        allocate(Hr_triple_out(1))
        Hr_triple_out(1)%attype= Hr_triple_in(1)%attype
        Hr_triple_out(1)%dist  =[Hr_triple_in(1)%dist  ]
        Hr_triple_out(1)%val   =[Hr_triple_in(1)%val   ]
        return
    endif
    !intializations and set first unique triples entry
    do i=1,size(hr_triple_in)
        at_triples(:,i)=Hr_triple_in(i)%attype
    enddo
    at_triples_unique(:,1)=at_triples(:,1)
    entry_unique(1)=1
    N_triple=0
    N_triple(1)=1
    ii=1
    !go through remaining triple entries and find out which are unique
    outer: do i=2,size(Hr_triple_in)
        do j=1,ii
            if(all(at_triples(:,i)==at_triples_unique(:,j)))then
                entry_unique(i)=j
                N_triple(j)=N_triple(j)+1
                cycle outer
            endif
        enddo
        ii=ii+1
        at_triples_unique(:,ii)=at_triples(:,i)
        entry_unique(i)=ii
        N_triple(ii)=N_triple(ii)+1
    enddo outer
    !fill Hr_triple_out
    allocate(Hr_triple_out(ii))
    do i=1,ii
        Hr_triple_out(i)%attype=at_triples_unique(:,i)
        allocate(Hr_triple_out(i)%dist(N_triple(i)))
        allocate(Hr_triple_out(i)%val(N_triple(i)))
    enddo
    N_triple=0
    do i=1,size(hr_triple_in)
        j=entry_unique(i)
        N_triple(j)=N_triple(j)+1
        Hr_triple_out(j)%val(N_triple(j))=Hr_triple_in(i)%val
        Hr_triple_out(j)%dist(N_triple(j))=Hr_triple_in(i)%dist
    enddo
end subroutine


subroutine reduce_Hr_pair(Hr_pair_in,Hr_pair_out)
    !subroutine which combines the entries with same atom types from Hr_pair_in into Hr_pair_out
    type(Hr_pair_single),intent(in)         ::  Hr_pair_in(:)
    type(Hr_pair),intent(out),allocatable   ::  Hr_pair_out(:)

    integer         ::  at_pairs       (2,size(Hr_pair_in))
    integer         ::  at_pairs_unique(2,size(Hr_pair_in))
    integer         ::  i,j,ii
    integer         ::  N_pair(size(hr_pair_in))
    integer         ::  entry_unique(size(hr_pair_in))

    if(size(Hr_pair_in)==1)then 
        !do nothing if only one entry is present
        allocate(Hr_pair_out(1))
        Hr_pair_out(1)%attype= Hr_pair_in(1)%attype
        Hr_pair_out(1)%dist  =[Hr_pair_in(1)%dist  ]
        Hr_pair_out(1)%val   =[Hr_pair_in(1)%val   ]
        return
    endif
    !intializations and set first unique pairs entry
    do i=1,size(hr_pair_in)
        at_pairs(:,i)=Hr_pair_in(i)%attype
    enddo
    at_pairs_unique(:,1)=at_pairs(:,1)
    entry_unique(1)=1
    N_pair=0
    N_pair(1)=1
    ii=1
    !go through remaining pair entries and find out which are unique
    outer: do i=2,size(Hr_pair_in)
        do j=1,ii
            if(all(at_pairs(:,i)==at_pairs_unique(:,j)))then
                entry_unique(i)=j
                N_pair(j)=N_pair(j)+1
                cycle outer
            endif
        enddo
        ii=ii+1
        at_pairs_unique(:,ii)=at_pairs(:,i)
        entry_unique(i)=ii
        N_pair(ii)=N_pair(ii)+1
    enddo outer
    !fill Hr_pair_out
    allocate(Hr_pair_out(ii))
    do i=1,ii
        Hr_pair_out(i)%attype=at_pairs_unique(:,i)
        allocate(Hr_pair_out(i)%dist(N_pair(i)))
        allocate(Hr_pair_out(i)%val(N_pair(i)))
    enddo
    N_pair=0
    do i=1,size(hr_pair_in)
        j=entry_unique(i)
        N_pair(j)=N_pair(j)+1
        Hr_pair_out(j)%val(N_pair(j))=Hr_pair_in(i)%val
        Hr_pair_out(j)%dist(N_pair(j))=Hr_pair_in(i)%dist
    enddo
end subroutine

subroutine reduce_Hr_tensor_pair(Hr_pair_in,Hr_pair_out)
    !subroutine which combines the entries with same atom types from Hr_pair_in into Hr_pair_out
    type(Hr_pair_single_tensor),intent(in)         ::  Hr_pair_in(:)
    type(Hr_pair_tensor),intent(out),allocatable   ::  Hr_pair_out(:)

    integer         ::  at_pairs       (2,size(Hr_pair_in))
    integer         ::  at_pairs_unique(2,size(Hr_pair_in))
    integer         ::  i,j,ii
    integer         ::  N_pair(size(hr_pair_in))
    integer         ::  entry_unique(size(hr_pair_in))

    if(size(Hr_pair_in)==1)then
        !do nothing if only one entry is present
        allocate(Hr_pair_out(1))
        Hr_pair_out(1)%attype= Hr_pair_in(1)%attype
        Hr_pair_out(1)%dist  =[ Hr_pair_in(1)%dist          ]
        Hr_pair_out(1)%val   =reshape([Hr_pair_in(1)%val,1.0d0],(/9,1/))
        Hr_pair_out(1)%bound =reshape([Hr_pair_in(1)%bound,1.0d0],(/3,1/))
        return
    endif
    !intializations and set first unique pairs entry
    do i=1,size(hr_pair_in)
        at_pairs(:,i)=Hr_pair_in(i)%attype
    enddo
    at_pairs_unique(:,1)=at_pairs(:,1)
    entry_unique(1)=1
    N_pair=0
    N_pair(1)=1
    ii=1
    !go through remaining pair entries and find out which are unique
    outer: do i=2,size(Hr_pair_in)
        do j=1,ii
            if(all(at_pairs(:,i)==at_pairs_unique(:,j)))then
                entry_unique(i)=j
                N_pair(j)=N_pair(j)+1
                cycle outer
            endif
        enddo
        ii=ii+1
        at_pairs_unique(:,ii)=at_pairs(:,i)
        entry_unique(i)=ii
        N_pair(ii)=N_pair(ii)+1
    enddo outer
    !fill Hr_pair_out
    allocate(Hr_pair_out(ii))
    do i=1,ii
        Hr_pair_out(i)%attype=at_pairs_unique(:,i)
        allocate(Hr_pair_out(i)%dist(N_pair(i)))
        allocate(Hr_pair_out(i)%val(9,N_pair(i)))
        allocate(Hr_pair_out(i)%bound(3,N_pair(i)))
    enddo
    N_pair=0
    do i=1,size(hr_pair_in)
        j=entry_unique(i)
        N_pair(j)=N_pair(j)+1
        Hr_pair_out(j)%val(:,N_pair(j))=Hr_pair_in(i)%val
        Hr_pair_out(j)%dist(N_pair(j))=Hr_pair_in(i)%dist
        Hr_pair_out(j)%bound(:,N_pair(j))=Hr_pair_in(i)%bound
    enddo
end subroutine

subroutine prt_Hr_pair(this,io,fmt_pre)
    class(Hr_pair),intent(in)               :: this
    integer,intent(in)                      :: io
    character(len=*),intent(in),optional    :: fmt_pre
    integer                                 :: i
    character(len=:),allocatable            :: pre

    if(present(fmt_pre))then    
        pre=fmt_pre
    else
        pre=""
    endif

    write(io,'('//pre//'A)') "Hamiltonian pair data:"
    write(io,'('//pre//'2X,A,2I6)') "atom types:",this%attype
    write(io,'('//pre//'2X,A,I6)') "number of interactions: ",size(this%dist)
    write(io,'('//pre//'2X,A)') "distance    value"
    do i=1,size(this%dist)
        write(io,'('//pre//'I6,E16.8)') this%dist(i),this%val(i)
    enddo
end subroutine

subroutine prt_Hr_pair_tensor(this,io,fmt_pre)
    class(Hr_pair_tensor),intent(in)        :: this
    integer,intent(in)                      :: io
    character(len=*),intent(in),optional    :: fmt_pre

    integer                                 :: i,j,shape_val(2)
    character(len=:),allocatable            :: pre
    character(len=100)                      :: form

    shape_val=shape(this%val)
    if(present(fmt_pre))then
        pre=fmt_pre
    else
        pre=""
    endif

    if (size(this%dist).ne.shape_val(2)) STOP "ERROR the size of the exchange matrix and the shell differs"

    write(form,'("(I6,",I10,"E16.8)")') shape_val(1)

    write(io,'('//pre//'A)') "Hamiltonian pair data:"
    write(io,'('//pre//'2X,A,2I6)') "atom types:",this%attype
    write(io,'('//pre//'2X,A,I6)') "number of interactions: ",size(this%dist)
    write(io,'('//pre//'2X,A)') "distance    value"
    do i=1,shape_val(2)
        write(io,form) this%dist(i),(this%val(j,i),j=1,shape_val(1))
    enddo
end subroutine

subroutine prt_Hr_triple(this,io,fmt_pre)
    class(Hr_triple),intent(in)             :: this
    integer,intent(in)                      :: io
    character(len=*),intent(in),optional    :: fmt_pre
    integer                                 :: i
    character(len=:),allocatable            :: pre

    if(present(fmt_pre))then    
        pre=fmt_pre
    else
        pre=""
    endif

    write(io,'('//pre//'A)') "Hamiltonian triple data:"
    write(io,'('//pre//'2X,A,2I6)') "atom types:",this%attype
    write(io,'('//pre//'2X,A,I6)') "number of interactions: ",size(this%dist)
    write(io,'('//pre//'2X,A)') "distance    value"
    do i=1,size(this%dist)
        write(io,'('//pre//'I6,E16.8)') this%dist(i),this%val(i)
    enddo
end subroutine
end module
