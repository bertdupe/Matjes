module m_deriv_null
use m_derived_types,only : lattice
use m_H_type
private
public :: t_deriv_l_null, t_deriv_r_null

type,extends(t_deriv_l)   :: t_deriv_l_null
contains
    procedure   :: get          => get_l_null
    procedure   :: get_single   => get_l_null_single
end type

type,extends(t_deriv_r)   :: t_deriv_r_null
contains
    procedure   :: get => get_r_null
    procedure   :: get_single   => get_r_null_single
end type

contains
    subroutine get_l_null(this,H,lat,vec)
        class(t_deriv_l_null),intent(in)   :: this
        class(t_H),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)
        
        vec=0.0d0
    end subroutine

    subroutine get_r_null(this,H,lat,vec)
        class(t_deriv_r_null),intent(in)   :: this
        class(t_H),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        vec=0.0d0
    end subroutine

    subroutine get_l_null_single(this,H,lat,site,vec)
        class(t_deriv_l_null),intent(in)    :: this
        class(t_H),intent(in)               :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        real(8),intent(inout)               :: vec(:)
        
        vec=0.0d0
    end subroutine

    subroutine get_r_null_single(this,H,lat,site,vec)
        class(t_deriv_r_null),intent(in)    :: this
        class(t_H),intent(in)               :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        real(8),intent(inout)               :: vec(:)
        
        vec=0.0d0
    end subroutine
end module
