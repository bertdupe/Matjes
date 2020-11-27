module m_user_info
interface welcome
       module procedure welcome_mpi
       module procedure welcome_serial
end interface welcome

interface user_info
       module procedure user_info_raw
end interface user_info

private
public :: welcome,user_info

contains
! simple module to write welcome messages for the different parts of the code

subroutine welcome_serial()
implicit none

write(6,'(/,a)') 'Bonjour!! you are now using the Kieler code'
write(6,'(a)') 'All the developpers hope that you will enjoy this moment.'
write(6,'(a)') 'If you have a problem, if you are happy with the code or if you wish to chat a bit'
write(6,'(a)') 'send an email to bertrand.dupe@gmail.com'
#ifdef CPP_VERSIONGIT
write(6,'(/2a)') "You are using the git-version: ",CPP_VERSIONGIT
#endif
write(6,'(/)') 


end subroutine welcome_serial

subroutine welcome_mpi(j)
implicit none
integer, intent(in) :: j

Call welcome_serial()
write(6,'(a,I10,2x,a)') 'Congratulations, you are now using',j,'processors'

end subroutine welcome_mpi

subroutine user_info_raw(io,time1,string,space)
implicit none
integer, intent(in) :: io
character(len=*), intent(in) :: string
logical, intent(in) :: space
real(kind=8), intent(inout) :: time1
! internal
character(len=30) :: forme
real(kind=8) :: time2

if (space) then
   forme='(/a/)'
else
   forme='(a)'
endif

if (time1.lt.1.0d-8) then
   call cpu_time(time1)
else
   call cpu_time(time2)
   write(io,'(a,E10.5,a)') 'duration: ',time2-time1,' seconds'
   time1=0.0d0
endif

#ifdef CPP_MPI
if (irank.eq.0)
#endif
   write(io,forme) string
#ifdef CPP_MPI
endif
#endif

end subroutine user_info_raw

end module m_user_info
