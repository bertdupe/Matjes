module m_deriv_null
use m_deriv_base
use m_derived_types,only : lattice
use m_H_type,only:  t_H_base
use m_work_ham_single, only:  work_ham_single

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
        class(t_deriv_l_null),intent(in)    :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        real(8),intent(inout)               :: vec(:)
        
        vec=0.0d0
    end subroutine

    subroutine get_r_null(this,H,lat,vec)
        class(t_deriv_r_null),intent(in)    :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        real(8),intent(inout)               :: vec(:)

        vec=0.0d0
    end subroutine

    subroutine get_l_null_single(this,H,lat,site,work,vec)
        class(t_deriv_l_null),intent(in)    :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)
        
        vec=0.0d0
    end subroutine

    subroutine get_r_null_single(this,H,lat,site,work,vec)
        class(t_deriv_r_null),intent(in)    :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)
        
        vec=0.0d0
    end subroutine
end module
