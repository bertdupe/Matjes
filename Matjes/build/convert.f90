module m_convert

interface convert
 module procedure convert_int_to_string
 module procedure string_plus_string_to_string
 module procedure string_plus_int_to_string
 module procedure string_plus_int_plus_string_to_string
end interface convert

private
public :: convert
contains

! convert a integer in a string of the good length
function convert_int_to_string(i)
implicit none
integer, intent(in) :: i
character(len=10) :: convert_int_to_string
! internal

convert_int_to_string=int_to_str(i)

end function

! merge 2 strings of different length into a string
function string_plus_string_to_string(str1,str2)
implicit none
character(len=*), intent(in) :: str1,str2
character(len=120) :: string_plus_string_to_string
! internal
integer :: l1,l2,l_tot

l1=len_trim(str1)
l2=len_trim(str2)
l_tot=l1+l2

if (l_tot.gt.120) then
   write(6,'(a)') 'error in the string_plus_string_to_string subroutine'
   stop
endif

write(string_plus_string_to_string,'(120a)') str1(1:l1),str2(1:l2)

end function

! merge a string and an integer into a string
function string_plus_int_to_string(str1,i)
implicit none
character(len=*), intent(in) :: str1
integer, intent(in) :: i
character(len=30) :: string_plus_int_to_string
! internal
integer :: l1,l2,l_tot
character(len=10) :: str2

l1=len_trim(str1)
str2=int_to_str(i)
l2=len_trim(str2)
l_tot=l1+l2

if (l_tot.gt.30) then
   write(6,'(a)') 'error in the string_plus_int_to_string subroutine'
   stop
endif

write(string_plus_int_to_string,'(30a)') str1(1:l1),str2(1:l2)

end function

! merge a string and an integer into a string
function string_plus_int_plus_string_to_string(str1,i,str2)
implicit none
character(len=*), intent(in) :: str1,str2
integer, intent(in) :: i
character(len=30) :: string_plus_int_plus_string_to_string
! internal
integer :: l1,l2,l_tot,l3
character(len=10) :: str3

l1=len_trim(str1)
str3=int_to_str(i)
l3=len_trim(str3)
l2=len_trim(str2)
l_tot=l1+l2+l3

if (l_tot.gt.30) then
   write(6,'(a)') 'error in the string_plus_int_to_string subroutine'
   stop
endif

write(string_plus_int_plus_string_to_string,'(30a)') str1(1:l1),str3(1:l3),str2(1:l2)

end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! private part of the module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function int_to_str(i)
implicit none
integer, intent(in) :: i
character(len=10) :: int_to_str
! internal
character(len=10) :: str

write(str,'(I10)') i

int_to_str=trim(adjustl(str))

end function

end module m_convert
