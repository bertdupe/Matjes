module m_convert

interface convert
 module procedure convert_int_to_string
 module procedure string_plus_string_to_string
 module procedure string_plus_int_to_string
 module procedure string_plus_int_plus_string_to_string
 module procedure string_plus_real_plus_string_to_string
 module procedure fourstring_plus_int_to_string
 module procedure threestring_plus_int_to_string
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
character(len=len_trim(str1)+len_trim(str2)) :: string_plus_string_to_string
! internal
integer :: l1,l2,l_tot
character(len=30) :: form

l1=len_trim(str1)
l2=len_trim(str2)
l_tot=l1+l2

if (l_tot.gt.120) then
   write(6,'(a)') 'error in the string_plus_string_to_string subroutine'
   stop
endif

write(form,'(a,I4,a)') '(',l_tot,'a)'
write(string_plus_string_to_string,form) str1(1:l1),str2(1:l2)

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
character(len=30) :: form

l1=len_trim(str1)
str2=int_to_str(i)
l2=len_trim(str2)
l_tot=l1+l2

if (l_tot.gt.30) then
   write(6,'(a)') 'error in the string_plus_int_to_string subroutine'
   stop
endif

write(form,'(a,I4,a)') '(',l_tot,'a)'
write(string_plus_int_to_string,form) str1(1:l1),str2(1:l2)

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
character(len=30) :: form

l1=len_trim(str1)
str3=int_to_str(i)
l3=len_trim(str3)
l2=len_trim(str2)
l_tot=l1+l2+l3

if (l_tot.gt.30) then
   write(6,'(a)') 'error in the string_plus_int_to_string subroutine'
   stop
endif

write(form,'(a,I4,a)') '(',l_tot,'a)'
write(string_plus_int_plus_string_to_string,form) str1(1:l1),str3(1:l3),str2(1:l2)

end function

! merge a string and a real into a string
function string_plus_real_plus_string_to_string(str1,x,str2)
implicit none
character(len=*), intent(in) :: str1,str2
real(kind=8), intent(in) :: x
character(len=50) :: string_plus_real_plus_string_to_string
! internal
integer :: l1,l2,l_tot,l3
character(len=50) :: str3
character(len=50) :: form

l1=len_trim(str1)
str3=real_to_str(x)
l3=len_trim(str3)
l2=len_trim(str2)
l_tot=l1+l2+l3

if (l_tot.gt.50) then
   write(6,'(a)') 'error in the string_plus_int_to_string subroutine'
   stop
endif

write(form,'(a,f10.5,a)') '(',l_tot,'a)'
write(string_plus_real_plus_string_to_string,form) str1(1:l1),str3(1:l3),str2(1:l2)

end function

! merge a 4 strings and an integer into a string
function fourstring_plus_int_to_string(str1,str2,str3,str4,i)
implicit none
character(len=*), intent(in) :: str1,str2,str3,str4
integer, intent(in) :: i
character(len=30) :: fourstring_plus_int_to_string
! internal
integer :: l1,l2,l_tot,l3,l4,l5
character(len=10) :: str5
character(len=30) :: form

l1=len_trim(str1)
l2=len_trim(str2)
l3=len_trim(str3)
l4=len_trim(str4)
str5=int_to_str(i)
l5=len_trim(str5)

l_tot=l1+l2+l3+l4+l5

if (l_tot.gt.30) then
   write(6,'(a)') 'error in the string_plus_int_to_string subroutine'
   stop
endif

write(form,'(a,I4,a)') '(',l_tot,'a)'
write(fourstring_plus_int_to_string,form) str1(1:l1),str2(1:l2),str3(1:l3),str4(1:l4),str5(1:l5)

end function


! merge a 3 strings and an integer into a string
function threestring_plus_int_to_string(str1,str2,i,str3)
implicit none
character(len=*), intent(in) :: str1,str2,str3
integer, intent(in) :: i
character(len=30) :: threestring_plus_int_to_string
! internal
integer :: l1,l2,l_tot,l3,l4
character(len=30) :: str4
character(len=100) :: form

l1=len_trim(str1)
l2=len_trim(str2)
l3=len_trim(str3)
str4=int_to_str(i)
l4=len_trim(str4)

l_tot=l1+l2+l3+l4

if (l_tot.gt.30) then
   write(6,'(a)') 'error in the string_plus_int_to_string subroutine'
   stop
endif

write(form,'(a,I4,a)') '(',l_tot,'a)'
write(threestring_plus_int_to_string,form) str1(1:l1),str2(1:l2),str4(1:l4),str3(1:l3)

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

function real_to_str(x)
implicit none
real(kind=8), intent(in) :: x
character(len=10) :: real_to_str
! internal
character(len=10) :: str

write(str,'(f10.5)') x

real_to_str=trim(adjustl(str))

end function

end module m_convert
