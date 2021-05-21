module m_deriv_rankN
use m_derived_types,only : lattice, dim_modes_inner
use m_deriv_base
use m_H_type, only: t_H_base
use m_work_ham_single, only:  work_ham_single

private
public :: t_deriv_l_N, t_deriv_r_N

type,extends(t_deriv_l)   :: t_deriv_l_N
    integer ::  order=0
contains
    procedure   :: get          => get_lN
    procedure   :: get_single   => get_lN_single
end type

type,extends(t_deriv_r)   :: t_deriv_r_N
    integer ::  order=0
contains
    procedure   :: get          => get_rN
    procedure   :: get_single   => get_rN_single
end type

contains
    subroutine get_lN(this,H,lat,vec)
        class(t_deriv_l_N),intent(in)   :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_r_red(lat,vec,this%order) 
    end subroutine

    subroutine get_rN(this,H,lat,vec)
        class(t_deriv_r_N),intent(in)   :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_l_red(lat,vec,this%order) 
    end subroutine

    subroutine get_lN_single(this,H,lat,site,work,vec)
        class(t_deriv_l_N),intent(in)       :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)

        real(8),allocatable :: vec_full(:)
        integer     :: dim_l_single
        integer     :: comp_l   !should this be saved in the type

        comp_l=H%mode_l%get_comp(this%order)
        Call H%mode_l%get_mode_single_size(this%order,dim_l_single)
        allocate(vec_full(dim_l_single))
        Call H%mult_r_single(site,comp_l,lat,work,vec_full)     !always has to be first and only component
        STOP "MULTIPLY VEC_FULL WITH EXCLUDED LEFT MODE"
        Call H%mode_l%reduce_site_vec(comp_l,vec_full,vec)
    end subroutine

    subroutine get_rN_single(this,H,lat,site,work,vec)
        class(t_deriv_r_N),intent(in)       :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)

        real(8),allocatable :: vec_full(:)
        integer     :: dim_r_single
        integer     :: comp_r   !should this be saved in the type

        comp_r=H%mode_r%get_comp(this%order)
        Call H%mode_r%get_mode_single_size(this%order,dim_r_single)
        allocate(vec_full(dim_r_single))
        Call H%mult_l_single(site,comp_r,lat,work,vec_full)     !always has to be first and only component
        STOP "MULTIPLY VEC_FULL WITH EXCLUDED RIGHT MODE"
        Call H%mode_r%reduce_site_vec(comp_r,vec_full,vec)
    end subroutine
end module
