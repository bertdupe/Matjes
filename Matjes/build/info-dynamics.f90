module m_info_dynamics

private
public :: choice_solver

contains

subroutine choice_solver(integtype)
implicit none
integer, intent(in) :: integtype
! internal variables


select case (integtype)
       case(3)
!-----------------------
! Heun
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'Heun integration method selected'
#else
        write(6,'(a)') 'Heun integration method selected'
#endif
       case(2)
!-----------------------
! SIA
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'SIA: Heun semi-implicit+projector integration method selected'
#else
        write(6,'(a)') 'SIA: Heun semi-implicit+projector integration method selected'
#endif
       case(4)
!-----------------------
! SIB
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'SIB: Heun semi-implicit+projector+NC integration method selected'
#else
        write(6,'(a)') 'SIB: Heun semi-implicit+projector+NC integration method selected'
#endif
       case(1)
!-----------------------
! Euler
#ifdef CPP_MPI1
        if (irank.eq.0) write(6,'(a)') 'Euler integration method selected'
#else
        write(6,'(a)') 'Euler integration method selected'

#endif
       case(5)
!-----------------------
! SIB with error correction
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'SIB with error correction method selected'
#else
        write(6,'(a)') 'SIB with error correction method selected'
#endif
       case(6)
!-----------------------
! SIB without temperature
#ifdef CPP_MPI
        if (irank.eq.0) write(6,'(a)') 'SIB: Heun semi-implicit+projector+NC integration method selected (T=0)'
#else
        write(6,'(a)') 'SIB: Heun semi-implicit+projector+NC integration method selected (T=0)'
#endif
       case default
end select

end subroutine choice_solver

end module m_info_dynamics
