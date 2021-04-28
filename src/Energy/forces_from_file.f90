module m_forces_from_file
use m_io_files_utils

type :: local_ham
   integer     :: index_atom = 0
   integer     :: index_neigh = 0
   real(8)     :: vector(3) = 0.0d0
   real(8)     :: H(3,3) = 0.0d0
end type

private
public :: get_ASR_file,get_forces_file

contains

! get the ASR from file
subroutine get_ASR_file(H,fname,offset)
real(8), intent(inout)         :: H(:,:)
character(len=*), intent(in)   :: fname
integer, intent(in)            :: offset(2)

integer                        :: io,i,j
type(local_ham)                :: ham

io=open_file_read(fname)
call get_ham_in_file(ham,(/0.0d0,0.0d0,0.0d0/),io)
call close_file(fname,io)

do j=1,3
   do i=1,3
      H(offset(1)+i,offset(2)+j)=ham%H(i,j)
   enddo
enddo

end subroutine


! get the forces from file
subroutine get_forces_file(H,fname,pos,offset)
real(8), intent(inout)         :: H(:,:)
character(len=*), intent(in)   :: fname
real(8), intent(in)            :: pos(3)
integer, intent(in)            :: offset(2)

! internal
integer                        :: io,i,j
type(local_ham)                :: ham

io=open_file_read(fname)
call get_ham_in_file(ham,pos,io)
call close_file(fname,io)

do j=1,3
   do i=1,3
      H(offset(1)+i,offset(2)+j)=ham%H(i,j)
   enddo
enddo

end subroutine



!!!!!!!
!! some private routines
!!!!!!!
subroutine get_ham_in_file(ham,pos_at,io)
integer, intent(in)            :: io
real(8), intent(in)            :: pos_at(3)
type(local_ham), intent(inout) :: ham
! internal variable
integer                        :: i,j,stat
real(8)                        :: norm_int
character(len=100)             :: str

rewind(io)
do

    read (io,'(a)',iostat=stat) str
    if (stat /= 0) exit
    str= trim(adjustl(str))
    if (len_trim(str)==0) cycle
    if (str(1:1) == '#' ) cycle
    backspace(io)
    read(io,*) ham%index_atom,ham%index_neigh
    read(io,*) (ham%vector(i),i=1,3)
    do j=1,3
        read(io,*) (ham%H(i,j),i=1,3)
    enddo

    norm_int=sum((pos_at-ham%vector)**2)

    if (norm_int.lt.1.0d-8) return

enddo

write(6,'(a,3(x,F12.6),a)') "Hamiltonian at position",pos_at," is 0.0"

end subroutine

end module m_forces_from_file
