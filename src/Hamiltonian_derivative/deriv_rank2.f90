module m_deriv_rank2
!use m_deriv
use m_derived_types,only : lattice
use m_H_type

private
public t_deriv_l_1, t_deriv_l_1_sym, t_deriv_r_1 

type,extends(t_deriv_l)   :: t_deriv_l_1
contains
    procedure   :: get => get_l1
end type

type,extends(t_deriv_l)   :: t_deriv_l_1_sym
contains
    procedure   :: get => get_l1_sym
end type

type,extends(t_deriv_r)   :: t_deriv_r_1
contains
    procedure   :: get => get_r1
end type

contains
    subroutine get_l1(this,H,lat,vec)
        class(t_deriv_l_1),intent(in)   :: this
        class(t_H),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_r(lat,vec)
    end subroutine

    subroutine get_r1(this,H,lat,vec)
        class(t_deriv_r_1),intent(in)     :: this
        class(t_H),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_l(lat,vec)
    end subroutine

    subroutine get_l1_sym(this,H,lat,vec)
        class(t_deriv_l_1_sym),intent(in)     :: this
        class(t_H),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_r(lat,vec)
        vec=vec*2.0d0
    end subroutine
end module
