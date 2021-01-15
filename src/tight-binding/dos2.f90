module m_dos2
!THIS ROUTINE REALLY NEEDS TO GET CLEANED AND UNITIFIED
!the nc and sc routines are very different in the way they add up the dos contributions

use m_TB_types, only: parameters_TB_IO_DOS
implicit none
private
public dos_nc, dos_sc

type dos_t
    integer             :: N_entry=0    !number of times dos has been added !obsolete?, manually norm anyways...
    real(8)             :: sigma        !smearing sigma
    real(8),allocatable :: Eval(:)      !energy values
    real(8),allocatable :: dos(:)       !summed dos
    real(8)             :: dE
    real(8)             :: E_ext(2)
contains
    procedure :: init => init_dos
    procedure :: print=> print_dos
end type

type,extends(dos_t) ::  dos_nc
contains
    procedure :: add  => add_dos_nc
end type


type,extends(dos_t) :: dos_sc
contains
    procedure :: add  => add_dos_sc
end type

contains

subroutine init_dos(this,io_dos)
    class(dos_t),intent(inout)              :: this
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    integer             ::  NE,iE

    Ne=int((io_dos%E_ext(2)-io_dos%E_ext(1))/io_dos%dE)+1
    allocate(this%dos(Ne),source=0.0d0)
    allocate(this%Eval(Ne))
    do iE=1,Ne
        this%Eval(iE)=(iE-1)*io_dos%dE+io_dos%E_ext(1)
    enddo
    this%N_entry=0
    this%sigma=io_dos%sigma
    this%dE=io_dos%dE
    this%E_ext=io_dos%E_ext
end subroutine

subroutine add_dos_nc(this,eigval)
    class(dos_nc),intent(inout)  :: this
    real(8),intent(in)          :: eigval(:)

    real(8)                     :: dos_loc(size(this%dos))
    integer                     :: iE

    if(.not.allocated(this%Eval)) STOP "Trying to add to dos, by type is not initialized"
    dos_loc=0.0d0
    do iE=1,size(this%Eval)
        Call get_dos(eigval,this%Eval(iE),dos_loc(iE),this%sigma)
    enddo
    dos_loc=dos_loc/real(size(eigval))
    this%dos=this%dos+dos_loc
    this%N_entry=this%N_entry+1
end subroutine

subroutine add_dos_sc(this,eigval,eigvec)
    class(dos_sc),intent(inout) :: this
    real(8),intent(in)          :: eigval(:)
    complex(8),intent(in)       :: eigvec(:,:)

    real(8)                     :: dos_loc(size(this%dos))
    integer                     :: i, i_start, n_state

    if(.not.allocated(this%Eval)) STOP "Trying to add to dos, by type is not initialized"
    dos_loc=0.0d0

    !only consider the states below the zero energy
    i_start=0
    do i=1,size(eigval)
        if(eigval(i)>0.0d0)then
            i_start=i 
            exit
        endif
    enddo
    if(i_start==0) STOP "No positive eigenvalues found ind add_dos_sc. For dos choose energies in positive branch."

    dos_loc=0.0d0
    do i=i_start,size(eigval)
        Call add_dos_1(this,eigval(i),eigvec(:,i),dos_loc)
    enddo
    this%dos=this%dos+dos_loc/real(size(eigval)-i_start+1)
    this%N_entry=this%N_entry+1
end subroutine

subroutine add_dos_1(this,val,eigvec,dos)
    class(dos_sc),intent(in)::  this
    real(8),intent(in)      ::  val
    complex(8),intent(in)   ::  eigvec(:)
    real(8),intent(inout)   ::  dos(:)
    
    real(8)                 ::  pref
    real(8),parameter   ::  dist_inc=5.0d0
    integer             ::  i_min,i_max
    integer             ::  dimH

    !need to add correct prefactors from eigvec
    pref=1.0d0
    dimH=size(eigvec)

    !u-part of BdG
    i_min=max(int(((val-this%sigma*dist_inc)-this%E_ext(1))/this%dE)+1,1)
    i_max=min(int(((val+this%sigma*dist_inc)-this%E_ext(1))/this%dE),size(this%Eval))
    pref=real(dot_product(eigvec(1:dimH/2),eigvec(1:dimH/2)),8)
    Call add_gauss(val,pref,this%Eval(i_min:i_max),dos(i_min:i_max),this%sigma)

    
    !v-part of BdG
    i_min=max(int(((-val-this%sigma*dist_inc)-this%E_ext(1))/this%dE)+1,1)
    i_max=min(int(((-val+this%sigma*dist_inc)-this%E_ext(1))/this%dE),size(this%Eval))
    pref=real(dot_product(eigvec(dimH/2+1:dimH),eigvec(dimH/2+1:dimH)),8)
    Call add_gauss(-val,pref,this%Eval(i_min:i_max),dos(i_min:i_max),this%sigma)

end subroutine

subroutine add_gauss(val,pref,E,dos,sigma)
    !add to dos:
    !gauss distribution from single energy point val into the spacing supplied by E with std. sigma and an additional prefactor pref 
    real(8),intent(in)     ::  val,pref,sigma
    real(8),intent(in)     ::  E(:)
    real(8),intent(inout)  ::  dos(:)
    real(8)                ::  tmp(size(dos))
    tmp=(val-E)**2
    tmp=-tmp*0.5d0/sigma/sigma
    tmp=exp(tmp)
    tmp=tmp/sqrt(2.0d0*3.14159265359d0)/sigma
    dos=dos+pref*tmp
end subroutine


subroutine print_dos(this,fname)
    class(dos_t),intent(inout)  :: this
    character(len=*),intent(in) :: fname
    real(8),allocatable         :: dos_loc(:)
    real(8)                     :: norm
    integer ::  io,i
    if(.not.allocated(this%Eval)) STOP "Trying to print dos, by type is not initialized"
    if(this%N_entry<1) STOP "Trying to print dos, no eigenvalue sets have been added"
    open(newunit=io,file=fname)
    dos_loc=this%dos
    norm=sum(dos_loc)
    dos_loc=dos_loc/norm
    do i=1,size(this%dos)
       write(io,'(2E16.8)') this%Eval(i),dos_loc(i)
    enddo
    close(io)
end subroutine


subroutine get_dos(val,E,res,sigma)
    real(8),intent(in)  ::  val(:),E,sigma
    real(8),intent(out) ::  res

    real(8)             ::  tmp(size(val,1))

    Call gauss_dist(val,E,sigma,tmp)
    res=sum(tmp)
end subroutine

subroutine gauss_dist(val,mu,sigma,res)
    use m_constants, only : pi
    !gauss distribution using 1/\sqrt{2*\pi*sigma^2}*e^{-(val-mu)^2/(2*sigma^2)} for all values in val respective to one given mu and sigma
    real(8),intent(in)  ::  val(:),mu,sigma
    real(8),intent(out) ::  res(size(val,1))
    integer             ::  i

    !only consider vals which are within [mu-dist_inc*sigma,mu+dist_inc*sigma], because all 
    !other results in small res  anyways and a too small exponent becomes numerically problematic
    real(8),parameter   ::  dist_inc=5.0d0

    integer             ::  i_min,i_max,Ne
    
    Ne=size(val,1)
    res=0.0d0
    if(val(Ne) <=mu-dist_inc*sigma.or.val(1) >=mu+dist_inc*sigma) return
    i_min=1
    do i=1,Ne
        if(val(i) >=mu-dist_inc*sigma)then
            i_min=i
            exit
        endif
    enddo
    i_max=i_min
    do i=Ne,i_min,-1
        if(val(i) <=mu+dist_inc*sigma)then
            i_max=i
            exit
        endif
    enddo
    if(i_min.eq.i_max) return
    res(i_min:i_max)=(val(i_min:i_max)-mu)*(val(i_min:i_max)-mu)
    res(i_min:i_max)=-res(i_min:i_max)*0.5d0/sigma/sigma
    res(i_min:i_max)=exp(res(i_min:i_max))
    res=res/sqrt(2.0d0*pi)/sigma
end subroutine
end module
