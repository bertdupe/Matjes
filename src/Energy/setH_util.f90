module m_setH_util
implicit none
interface get_coo
    module procedure get_coo_r
    module procedure get_coo_c
end interface
    real(8),parameter,private   ::  acc=1.0d-10 !exclude all entries smaller by factor acc than maximal size

contains
subroutine get_coo_r(M,val,ind)
    !small subroutine to get indices and values from 2D matrix M to represent it in coo format
    real(8),intent(in)                  :: M(:,:)
    integer,intent(inout),allocatable   :: ind(:,:)
    real(8),allocatable                 :: val(:)
    
    integer     :: N,i1,i2,ii
    real(8)     :: maxv,thres

    if(allocated(ind)) deallocate(ind)
    if(allocated(val)) deallocate(val)
    maxv=maxval(abs(M))
    thres=maxv*acc
    N=count(abs(M)>thres)
    allocate(ind(2,N),source=0)
    allocate(val(N),source=0.0d0)
    ii=0
    do i2=1,size(M,2)
        do i1=1,size(M,1)
            if(abs(M(i1,i2))>thres)then
                ii=ii+1
                val(ii)=M(i1,i2)
                ind(:,ii)=[i1,i2]
            endif
        enddo
    enddo
end subroutine

subroutine get_coo_c(M,val,ind)
    !small subroutine to get indices and values from 2D matrix M to represent it in coo format
    complex(8),intent(in)               :: M(:,:)
    integer,intent(inout),allocatable   :: ind(:,:)
    complex(8),allocatable              :: val(:)
    
    integer     :: N,i1,i2,ii
    real(8)     :: maxv,thres

    if(allocated(ind)) deallocate(ind)
    if(allocated(val)) deallocate(val)
    maxv=maxval(abs(M))
    thres=maxv*acc
    N=count(abs(M)>thres)
    allocate(ind(2,N),source=0)
    allocate(val(N),source=(0.0d0,0.0d0))
    ii=0
    do i2=1,size(M,2)
        do i1=1,size(M,1)
            if(abs(M(i1,i2))>thres)then
                ii=ii+1
                val(ii)=M(i1,i2)
                ind(:,ii)=[i1,i2]
            endif
        enddo
    enddo
end subroutine


function ind(dimmodes,i_entry)
    !function which gets get correct index for the Hamiltonian in the (1:product(dimmodes))-space
    integer :: dimmodes(:)  !dim_mode for each entry
    integer :: i_entry(:)   !index in space of each dim_mode
    integer :: ind
    !logical
    integer :: N,i

    N=size(dimmodes)
    if(size(i_entry)/=N) STOP "ind input variables of differing rank"
    if(any(i_entry>dimmodes).or.any(i_entry<1)) STOP "index dimmodes,i_entry input inconsistend"
    ind=1
    do i=1,N
        ind=ind+(i_entry(i)-1)*product(dimmodes(:i-1))
    enddo
end function


end module
