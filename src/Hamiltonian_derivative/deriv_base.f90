module m_deriv_base
use m_H_type, only: t_H_base
use m_derived_types, only : lattice

private
public :: t_deriv_l, t_deriv_r, t_deriv_base 

type,abstract    ::  t_deriv_l
contains
    procedure(int_deriv_l_get),deferred :: get
    procedure(int_deriv_l_get_single),deferred :: get_single
end type

type,abstract    ::  t_deriv_r
contains
    procedure(int_deriv_r_get),deferred :: get
    procedure(int_deriv_r_get_single),deferred :: get_single
end type


type    :: t_deriv_base
    class(t_deriv_l),allocatable  :: l
    class(t_deriv_r),allocatable  :: r
contains
    procedure :: copy
    procedure :: get => get_deriv
    procedure :: get_single => get_deriv_single
end type

abstract interface
    subroutine int_deriv_l_get(this,H,lat,vec)
        import t_deriv_l, t_H_base, lattice
        class(t_deriv_l),intent(in) :: this
        class(t_H_base),intent(in)  :: H
        type(lattice),intent(in)    :: lat
        real(8),intent(inout)       :: vec(:)
    end subroutine
    subroutine int_deriv_r_get(this,H,lat,vec)
        import t_deriv_r, t_H_base, lattice
        class(t_deriv_r),intent(in) :: this
        class(t_H_base),intent(in)  :: H
        type(lattice),intent(in)    :: lat
        real(8),intent(inout)       :: vec(:)
    end subroutine
    subroutine int_deriv_l_get_single(this,H,lat,site,vec)
        import t_deriv_l, t_H_base, lattice
        class(t_deriv_l),intent(in) :: this
        class(t_H_base),intent(in)  :: H
        type(lattice),intent(in)    :: lat
        integer,intent(in)          :: site
        real(8),intent(inout)       :: vec(:)
    end subroutine
    subroutine int_deriv_r_get_single(this,H,lat,site,vec)
        import t_deriv_r, t_H_base, lattice
        class(t_deriv_r),intent(in) :: this
        class(t_H_base),intent(in)  :: H
        type(lattice),intent(in)    :: lat
        integer,intent(in)          :: site
        real(8),intent(inout)       :: vec(:)
    end subroutine
end interface
contains

    subroutine copy(this,deriv_out)
        class(t_deriv_base),intent(in)     :: this
        class(t_deriv_base),intent(inout)  :: deriv_out

        if(.not.allocated(this%l)) ERROR STOP "CANNOT COPY derivative type as left source is not allocated"
        if(.not.allocated(this%r)) ERROR STOP "CANNOT COPY derivative type as right source is not allocated"
        if(allocated(deriv_out%l)) deallocate(deriv_out%l)
        if(allocated(deriv_out%r)) deallocate(deriv_out%r)
        allocate(deriv_out%l,source=this%l)
        allocate(deriv_out%r,source=this%r)
    end subroutine

    subroutine get_deriv(this,H,lat,vec,tmp)
        class(t_deriv_base),intent(in)  :: this             !derive type with set procedure and order to derive with respect to
        class(t_H_base),intent(in)      :: H                !Hamiltonian that is derivated
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)           !add derivative to this vector (
        real(8),intent(inout)           :: tmp(size(vec))   !to prevent constant allocation/deallocation

        Call this%l%get(H,lat,tmp)
        vec=vec+tmp
        Call this%r%get(H,lat,tmp)
        vec=vec+tmp
    end subroutine

    subroutine get_deriv_single(this,H,lat,site,vec,tmp)
        class(t_deriv_base),intent(in)  :: this             !derive type with set procedure and order to derive with respect to
        class(t_H_base),intent(in)      :: H                !Hamiltonian that is derivated
        type(lattice),intent(in)        :: lat
        integer,intent(in)              :: site
        real(8),intent(inout)           :: vec(:)           !add derivative to this vector (
        real(8),intent(inout)           :: tmp(size(vec))   !to prevent constant allocation/deallocation

        Call this%l%get_single(H,lat,site,tmp)
        vec=vec+tmp
        Call this%r%get_single(H,lat,site,tmp)
        vec=vec+tmp
    end subroutine

end module
