      module m_set_temp
      interface ini_temp
       module procedure ini_temp_normal
!#ifdef CPP_MPI
!       module procedure ini_temp_mpi_ghost
!       module procedure ini_temp_mpi
!#endif
      end interface ini_temp
      contains
!
! subroutine that initializes and update the temperature in serial calculations
! everything is defined as local variable
! kt is an out
! the rest is purely in
!
!#ifdef CPP_MPI
!      subroutine ini_temp_mpi(kt,kTfin,kTini,isize,irank,nRepProc,i_print_W)
!      implicit none
!      real(kind=8), intent(out) :: kt(:)
!      real(kind=8), intent(in) :: kTfin,kTini
!      integer, intent(in) :: isize,irank,nRepProc
!      logical, intent(in) :: i_print_W
!! internal variables
!     integer :: i
!     logical :: i_geometric
!
!     i_geometric=.False.
!
!     inquire (file='geometric',exist=i_geometric)
!
!     if (isize.gt.1) then
!        if (i_geometric) then
!            do i=1,nRepProc
!               kt(i)=(kTfin/kTini)**(dble(i-1+irank*nRepProc)/dble(isize-1))*kTini
!            enddo
!        else
!            do i=1,nRepProc
!               kt(i)=(ktfin-kTini)/dble(isize-1)*dble(i-1+irank*nRepProc)+ktini
!            enddo
!        endif
!     else
!        if (i_print_W) then
!               write(6,'(/,a)') 'WARNING: only 1 temperature will be used in the MC'
!               write(6,'(a,/)') 'T=Tinitial'
!        endif
!        kt(1)=kTini
!     endif
!
!     end subroutine ini_temp_mpi
!#endif

      subroutine ini_temp_normal(kt,kTfin,kTini,isize,i_print_W)
      implicit none
      real(kind=8), intent(out) :: kt(:)
      real(kind=8), intent(in) :: kTfin,kTini
      integer, intent(in) :: isize
      logical, intent(in) :: i_print_W
! internal variables
     integer :: i
     logical :: i_geometric

     i_geometric=.False.

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

!#ifdef CPP_MPI
!     subroutine ini_temp_mpi_ghost(kt,kTfin,kTini,irank,isize,n_ghost,i_print_W)
!     implicit none
!     real(kind=8), intent(out) :: kt(:),kTfin,kTini
!     integer, intent(in) :: irank,n_ghost,isize
!     logical, intent(in) :: i_print_W
!! internal variables
!     integer :: i
!     logical :: i_geometric
!
!     i_geometric=.False.
!
!     inquire (file='geometric',exist=i_geometric)
!
!     if (isize.gt.1) then
!        if (i_geometric) then
!            do i=1,isize
!               kt(i)=(kTfin/kTini)**(dble(irank/n_ghost)/dble(isize-1))*kTini
!            enddo
!        else
!            do i=1,isize
!               kt(i)=(ktfin-kTini)/dble(isize-1)*dble(irank/n_ghost)+ktini
!            enddo
!        endif
!     else
!        if (i_print_W) then
!               write(6,'(/,a)') 'WARNING: only temperature will be used in the MC'
!               write(6,'(a,/)') 'T=Tinitial'
!        endif
!        kt=kTini
!     endif
!
!     end subroutine ini_temp_mpi_ghost
!#endif

     end module m_set_temp
