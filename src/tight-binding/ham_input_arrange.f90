module m_ham_arrange
use m_types_tb_h_inp 
use m_derived_types, only: lattice

private
public :: Htb_inp, Htb_atpair, Htb_dist

type Htb_inp
    type(Htb_atpair), allocatable :: at(:)
contains
    procedure ::    set=> set_ham
end type 

type Htb_atpair
    integer     :: atpair(2)
    integer     :: Nspin(2)
    integer     :: Norb(2)
    type(Htb_dist), allocatable :: dist(:)
end type

type Htb_dist
    integer ::  dist=-1
    complex(8),allocatable  ::  Hloc(:,:)
contains
    procedure ::     add => add_hampair_entry
end type

contains 

subroutine set_ham(ham,hop,lat,nspin)
    use m_sort
    class(Htb_inp),intent(out)       :: Ham
    type(TB_hopping),intent(in)         :: hop(:)
    integer,intent(in)                  :: nspin
    type(lattice),intent(in)            :: lat

    !data which is extracted from hop and shuffeled around
    integer     :: attype(2,size(hop))
    integer     :: orbital(2,size(hop))
    integer     :: spin(2,size(hop))
    integer     :: dist(size(hop))
    complex(8)  :: val(size(hop))

    integer,allocatable :: Norb_attype(:)   !number of orbitals for each atom type

    integer,allocatable :: atpair_N(:)      !number of entries of each unique atom pair [1:Natpair]->[1,size(hop)]
    integer,allocatable :: atpair_unique(:) !indices of unique atom pairs in attype array [1:Natpair]->[1,size(hop)] (temporary)

    !local parameters for sorting
    integer             :: N
    integer,allocatable :: ind(:)
    integer             :: atpair_unique_id(size(hop))  !index of unique atom pair in atpair_unique [1:size(hop)]->[1:Natpair]

    integer,allocatable :: Ndist_unique(:)  !number of unique distances for each unique atom type pair
    integer :: Natpair !number of unique atom pairs
    integer :: Ndim !dimension of local Hamiltonian

    integer :: i,j
    integer :: i1,i2
    integer :: i_dist, tmp_dist

    !get data from derived type
    do i=1,size(hop)
        attype(:,i)=hop(i)%attype
        orbital(:,i)=hop(i)%orbital
        spin(:,i)=hop(i)%spin
        dist(i)=hop(i)%dist
        val(i)=cmplx(hop(i)%val,0.0d0,8)
    enddo 
    allocate(Norb_attype(lat%cell%n_attype),source=0)
    do i=1,size(lat%cell%atomic)
        Norb_attype(lat%cell%atomic(i)%type_id)=lat%cell%atomic(i)%orbitals
    enddo

    !arrange by atom type pairs so that the same combinations are next to each other
    Call get_unique_pair(attype,atpair_N,atpair_unique,atpair_unique_id)
    !first save unique atpair sets
    Natpair=size(atpair_N)
    allocate(Ham%at(Natpair))
    do i=1,Natpair
        Ham%at(i)%atpair=attype(:,atpair_unique(i))
        Ham%at(i)%Nspin=nspin
        Ham%at(i)%Norb=Norb_attype(attype(:,atpair_unique(i)))
    enddo
    deallocate(atpair_unique)
    !actual sorting
    N=size(hop)
    allocate(ind(N),source=0)
    Call sort(N,atpair_unique_id,ind)
    Call rearange_arrays(ind,orbital,spin,val,dist,attype)
    deallocate(ind)

    !sort each atom type pair by distance
    allocate(Ndist_unique(Natpair))
    do i=1,Natpair
        i1=1+sum(atpair_N(1:i-1))
        i2=sum(atpair_N(1:i))
        N=atpair_N(i)
        allocate(ind(N),source=0)
        Call sort(N,dist(i1:i2),ind)
        Call rearange_arrays(ind,orbital(:,i1:i2),spin(:,i1:i2),val(i1:i2))
        Ndist_unique(i)=1+count([(dist(j-1)<dist(j),j=i1+1,i2)]) !number of unique distances
        deallocate(ind)
    enddo

    !allocate Ham local Hamiltonians 
    do i=1,Natpair
        allocate(Ham%at(i)%dist(Ndist_unique(i)))
        do j=1,Ndist_unique(i)
            Ndim=product(Ham%at(i)%Nspin)*product(Ham%at(i)%Norb)
            allocate(Ham%at(i)%dist(j)%Hloc(Ndim,Ndim),source=(0.d0,0.0d0))
        enddo
    enddo

    !all all data into the Ham type
    do i=1,Natpair
        i_dist=0; tmp_dist=-1
        i1=1+sum(atpair_N(1:i-1))
        i2=sum(atpair_N(1:i))
        do j=i1,i2
            if(dist(j)/=tmp_dist) i_dist=i_dist+1 
            tmp_dist=dist(j)
            ham%at(i)%dist(i_dist)%dist=dist(j)
            Call ham%at(i)%dist(i_dist)%add(orbital(:,j),spin(:,j),val(j),ham%at(i)%Norb,ham%at(i)%Nspin)
        enddo
    enddo

end subroutine

subroutine rearange_arrays(ind,orbital,spin,val,dist,attype)
    !it might make sense to implement some better rearranging routine
    integer,intent(in)              :: ind(:)
    integer,intent(inout)           :: orbital(2,size(ind))
    integer,intent(inout)           :: spin   (2,size(ind))
    complex(8),intent(inout)        :: val      (size(ind))
    integer,intent(inout),optional  :: attype (2,size(ind))
    integer,intent(inout),optional  :: dist     (size(ind))

    orbital=orbital(:,ind)
    spin=spin(:,ind)
    val=val(ind)
    if(present(attype)) attype=attype(:,ind)
    if(present(dist)) dist=dist(ind)
end subroutine
    
subroutine get_unique_pair(atpair,atpair_N,atpair_unique,atpair_unique_id)
    integer,intent(in)              :: atpair(:,:)  !first dimension should be 2
    integer,allocatable,intent(out) :: atpair_unique(:)
    integer,allocatable,intent(out) :: atpair_N(:)
    integer,intent(out)             :: atpair_unique_id(size(atpair,2))

    integer             :: atpair_unique_N (size(atpair,2))
    integer             :: Nimplicit_entry (size(atpair,2))
    integer,allocatable :: tmp(:)
    integer             :: N_unique, Npair
    integer             :: i,j

    if(size(atpair,1)/=2) STOP "get unique pair atpair has to be of shape 2,:"

    Npair=size(atpair,2)
    N_unique=1
    atpair_unique_id=1
    atpair_unique_N=0
    atpair_unique_N(1)=1
    allocate(atpair_unique(Npair),source=1)
    outer:do i=2,Npair
        do j=1,N_unique
            if(all(atpair(:,i)==atpair(:,atpair_unique(j))))then
               atpair_unique_id(i)=j
               atpair_unique_N(j)=atpair_unique_N(j)+1
               cycle outer
            endif
        enddo
        N_unique=N_unique+1
        atpair_unique(N_unique)=i
        atpair_unique_N(N_unique)=atpair_unique_N(N_unique)+1
        atpair_unique_id(i)=N_unique
    enddo outer
    allocate(atpair_N,source=atpair_unique_N(1:N_unique))
    Call move_alloc(atpair_unique,tmp)
    allocate(atpair_unique,source=tmp(1:N_unique))
    deallocate(tmp)
end subroutine

recursive subroutine add_hampair_entry(this,orbital,spin,val_in,Norb,Nspin)
    class(Htb_dist),intent(inout)   ::  this
    integer,intent(in)                      ::  orbital(2),spin(2)
    integer,intent(in)                      ::  Norb(2),Nspin(2)
    complex(8),intent(in)                   ::  val_in
    integer                                 ::  ind(2)

    complex(8)  :: val
    integer     :: i

    if(all(spin==0))then
        !if both spins are 0, add entry to both spin channels
        !probably inefficient, but I don't care right now
        do i=1,Nspin(1)
            Call this%add(orbital,[i,i],val_in,Norb,Nspin)
        enddo
        return
    endif
    if(any(orbital<1).or.any(orbital>Norb)) STOP "INVALID TB-HAMILTONIAN ORBITAL INPUT, check input"
    if(any(spin<1).or.any(spin>Nspin)) STOP "INVALID TB-HAMILTONIAN SPIN INPUT, check input"

    val=val_in
    ind=(orbital-1)*Nspin+spin
    if(this%dist==0.and.ind(1)>ind(2)) val=conjg(val_in)    !make sure hermitian ... might give problems with deltas...  (this%dist==0 -> onsite)

    this%Hloc(ind(1),ind(2))=this%Hloc(ind(1),ind(2))+val
end subroutine



end module
