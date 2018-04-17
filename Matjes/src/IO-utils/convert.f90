module m_convert
private
public :: convert

interface convert
 module procedure convert_int_string
 module procedure convert_string_string
end interface convert
contains

! convert a integer in a string of the good length
function convert_int_string(i)
implicit none
integer, intent(in) :: i
character(len=10) :: convert_int_string
! internal
character(len=10) :: str

write(str,'(I10)') i

convert_int_string=trim(adjustl(str))

end function

! merge 2 strings of different length into a string
function convert_string_string(str1,str2)
implicit none
character(len=*), intent(in) :: str1,str2
character(len=10) :: convert_string_string
! internal
integer :: l1,l2,l_tot

l1=len_trim(str1)
l2=len_trim(str2)
l_tot=l1+l2

if (l_tot.gt.10) then
   write(6,'(a)') 'error in the convert_string_string subroutine'
   stop
endif

write(convert_string_string,'(10a)') str1(1:l1),str2(1:l2)

end function

end module m_convert
