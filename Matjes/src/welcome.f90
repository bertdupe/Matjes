       module m_welcome
       interface welcome
        module procedure welcome_mpi
        module procedure welcome_serial
       end interface welcome

       contains
! simple module to write welcome messages for the different parts of the code

       subroutine welcome_serial()
       implicit none

       write(6,'(/,a)') 'Bonjour!! you are now using the Kieler code'
       write(6,'(a)') 'All the developpers hope that you will enjoy this moment.'
       write(6,'(a)') 'If you have a problem, if you are happy with the code or if you wish to chat a bit'
       write(6,'(a)') 'send an email to bertrand.dupe@gmail.com'

       end subroutine welcome_serial

       subroutine welcome_mpi(j)
       implicit none
       integer, intent(in) :: j

       write(6,'(/,a)') 'Bonjour!! you are now using the Kieler code'
       write(6,'(a)') 'All the developpers hope that you will enjoy this moment.'
       write(6,'(a)') 'If you have a problem, if you are happy with the code or if you wish to chat a bit'
       write(6,'(a,/)') 'send an email to bertrand.dupe@gmail.com'
       write(6,'(a,I10,2x,a)') 'Congratulations, you are now using',j,'processors'

       end subroutine welcome_mpi


       end module m_welcome
