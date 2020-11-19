module m_input_types
use m_constants, only : pi
implicit none

private
public :: MC_input

type MC_input
    !parameters what is run how often
    integer     :: n_Tsteps=1
    integer     :: n_sizerelax=1
    integer     :: n_thousand=1000
    integer     :: restart_MC_steps=0
    integer     :: T_relax=1
    integer     :: T_auto=1
    integer     :: Total_MC_Steps=1000

    logical     :: i_restart=.false.
    logical     :: print_relax=.false.
    logical     :: Cor_log=.false.

    real(8)     :: cone=pi

    logical     :: ising=.false.
    logical     :: underrelax=.false.
    logical     :: overrelax=.false.
    logical     :: equi=.false.
    logical     :: sphere=.false.
end type


end module
