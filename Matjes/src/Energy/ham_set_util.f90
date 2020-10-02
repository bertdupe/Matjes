module m_setH_util
implicit none
contains
subroutine get_coo(M,val,ind)
    !small subroutine to get indices and values from 2D matrix M to represent it in coo format
    real(8),intent(in)                  ::  M(:,:)
    integer,intent(inout),allocatable   ::  ind(:,:)
    real(8),allocatable                 ::  val(:)
    
    integer     :: N,i1,i2,ii
    real(8)     :: maxv,thres
    real(8),parameter   ::  acc=1.0d-10 !exclude all entries smaller by factor acc than maximal size

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

end module
