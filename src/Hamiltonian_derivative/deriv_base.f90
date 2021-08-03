module m_deriv_base
use m_H_type, only: t_H_base
use m_derived_types, only : lattice
use m_work_ham_single, only:  work_ham_single, work_mode

private
public t_deriv

type    :: t_deriv
    integer :: order=0
    procedure(int_deriv_get),pointer        :: l=>uninitialized
    procedure(int_deriv_get),pointer        :: r=>uninitialized
    procedure(int_deriv_get_single),pointer :: l_single=>uninitialized_single
    procedure(int_deriv_get_single),pointer :: r_single=>uninitialized_single
contains
    procedure :: copy
    procedure :: get => get_deriv
    procedure :: get_single => get_deriv_single
    procedure :: mv
end type

abstract interface
    subroutine int_deriv_get(this,H,lat,vec,work)
        import t_deriv, t_H_base, lattice, work_mode
        class(t_deriv),intent(in)       :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)
        type(work_mode),intent(inout)   :: work
    end subroutine
    subroutine int_deriv_get_single(this,H,lat,site,work,vec)
        import  t_deriv,t_H_base, lattice, work_ham_single
        class(t_deriv),intent(in)           :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) ::  work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)
    end subroutine
end interface
contains

subroutine copy(this,deriv_out)
    class(t_deriv),intent(in)       :: this
    class(t_deriv),intent(inout)    :: deriv_out

    deriv_out%l=>this%l
    deriv_out%r=>this%r
    deriv_out%l_single=>this%l_single
    deriv_out%r_single=>this%r_single
end subroutine

subroutine mv(this,deriv_out)
    class(t_deriv),intent(in)       :: this
    class(t_deriv),intent(inout)    :: deriv_out

    deriv_out%l=>this%l
    deriv_out%r=>this%r
    deriv_out%l_single=>this%l_single
    deriv_out%r_single=>this%r_single
end subroutine

subroutine get_deriv(this,H,lat,vec,work)
    class(t_deriv),intent(in)       :: this             !derive type with set procedure and order to derive with respect to
    class(t_H_base),intent(in)      :: H                !Hamiltonian that is derivated
    type(lattice),intent(in)        :: lat
    real(8),intent(inout)           :: vec(:)           !add derivative to this vector (
    type(work_mode),intent(inout)   :: work

    Call this%l(H,lat,vec,work)
    Call this%r(H,lat,vec,work)
end subroutine

subroutine get_deriv_single(this,H,lat,site,work,vec,tmp)
    class(t_deriv),intent(in)           :: this             !derive type with set procedure and order to derive with respect to
    class(t_H_base),intent(in)          :: H                !Hamiltonian that is derivated
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: site
    type(work_ham_single),intent(inout) :: work             !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
    real(8),intent(inout)               :: vec(:)           !add derivative to this vector (
    real(8),intent(inout)               :: tmp(size(vec))   !to prevent constant allocation/deallocation !also remove with beta=1 as in get_deriv-case?

    Call this%l_single(H,lat,site,work,tmp)
    vec=vec+tmp
    Call this%r_single(H,lat,site,work,tmp)
    vec=vec+tmp
end subroutine

subroutine uninitialized(this,H,lat,vec,work)
    class(t_deriv),intent(in)       :: this
    class(t_H_base),intent(in)      :: H
    type(lattice),intent(in)        :: lat
    real(8),intent(inout)           :: vec(:)
    type(work_mode),intent(inout)   :: work
   
    ERROR STOP "DERIVATIVE POINTER HAS NOT BEEN SET"
end subroutine

subroutine uninitialized_single(this,H,lat,site,work,vec)
    class(t_deriv),intent(in)           :: this
    class(t_H_base),intent(in)          :: H
    type(lattice),intent(in)            :: lat
    integer,intent(in)                  :: site
    type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
    real(8),intent(inout)               :: vec(:)
    
    ERROR STOP "SINGLE DERIVATIVE POINTER HAS NOT BEEN SET"
end subroutine


end module
