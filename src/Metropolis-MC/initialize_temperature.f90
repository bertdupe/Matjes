module m_set_temp
implicit none
interface ini_temp
    module procedure ini_temp_normal
end interface ini_temp
contains
!
! subroutine that initializes and update the temperature in serial calculations
! everything is defined as local variable
! kt is an out
! the rest is purely in
!

subroutine ini_temp_normal(kt,kTfin,kTini,isize,i_print_W)
    implicit none
    real(8),intent(inout),allocatable   :: kt(:)
    real(8),intent(in)                  :: kTfin,kTini
    integer,intent(in)                  :: isize
    logical,intent(in)                  :: i_print_W
    ! internal variables
    integer :: i
    logical :: i_geometric

    if(.not.allocated(kt))then
        allocate(kt(isize))
    else
        if(size(kt)/=isize) STOP "Allocated kt, but size is different, probably a programming mistake..."
    endif

    inquire (file='geometric',exist=i_geometric)

    if (isize.gt.1) then
       if (i_geometric) then
           do i=1,isize
              kt(i)=(kTfin/kTini)**(dble(i-1)/dble(isize-1))*kTini
           enddo
       else
           do i=1,isize
              kt(i)=(ktfin-kTini)/dble(isize-1)*dble(i-1)+ktini
           enddo
       endif
    else
       if (i_print_W) then
              write(6,'(/,a)') 'WARNING: only 1 temperature will be used in the MC'
              write(6,'(a,/)') 'T=Tinitial'
       endif
       kt(1)=kTini
    endif
end subroutine ini_temp_normal

end module m_set_temp
