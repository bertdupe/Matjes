module m_pos_op
implicit none

type    :: pos_op
    complex(8),allocatable  ::  pos(:,:)
contains
    procedure   ::  init
end type
contains
subroutine init(this,dimH)
    class(pos_op),intent(inout)     :: this
    integer,intent(in)              :: dimH
    integer ::  i

    allocate(this%pos(dimH,dimH),source=(0.0d0,0.0d0))
    do i=1,dimH
        this%pos(i,i)=(1.0d0,0.0d0)
    enddo
end subroutine

end module
