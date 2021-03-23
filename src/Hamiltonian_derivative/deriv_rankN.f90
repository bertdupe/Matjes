module m_deriv_rankN
use m_derived_types,only : lattice
use m_H_type

private
public :: t_deriv_l_N, t_deriv_r_N

type,extends(t_deriv_l)   :: t_deriv_l_N
    integer ::  order=0
contains
    procedure   :: get => get_lN
end type

type,extends(t_deriv_r)   :: t_deriv_r_N
    integer ::  order=0
contains
    procedure   :: get => get_rN
end type

contains
    subroutine get_lN(this,H,lat,vec)
        class(t_deriv_l_N),intent(in)   :: this
        class(t_H),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_r_red(lat,vec,this%order) 
    end subroutine

    subroutine get_rN(this,H,lat,vec)
        class(t_deriv_r_N),intent(in)   :: this
        class(t_H),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_l_red(lat,vec,this%order) 
    end subroutine
end module
