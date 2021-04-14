module m_deriv_rankN
use m_derived_types,only : lattice, dim_modes_inner
use m_deriv_base
use m_H_type, only: t_H_base

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
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_r_red(lat,vec,this%order) 
    end subroutine

    subroutine get_rN(this,H,lat,vec)
        class(t_deriv_r_N),intent(in)   :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_l_red(lat,vec,this%order) 
    end subroutine

    subroutine get_lN_single(this,H,lat,site,vec)
        class(t_deriv_l_N),intent(in)   :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        integer,intent(in)              :: site
        real(8),intent(inout)           :: vec(:)

        integer                         :: ind_vec(dim_modes_inner(this%order))

        integer,allocatable             :: ind_l(:)   !indices of considered mode that have to be reduced
        integer,allocatable             :: ind_r(:)   !indices of considered mode on which H has to be applied
        integer                         :: N_r

        real(8),allocatable             :: vec_l(:)
        real(8),allocatable             :: vec_r(:)
        real(8),allocatable             :: vec_applied(:)
        integer     :: comp_l   !should this be saved in the type
        integer     :: i

        comp_l=H%mode_l%get_comp(this%order)
        ind_vec=[((site-1)*dim_modes_inner(this%order)+i,i=1,dim_modes_inner(this%order))]
        !get all vector indices on the left side that contribute in the reduction to ind_vec
        Call H%mode_l%get_ind_site(comp_l,site,ind_l)   
        !get all vector indices on the right side that have contributions at ind_l after applying H
        N_r=size(ind_l)*10 !arbitrary, and think of something better with external work arrays
        allocate(ind_r(N_r),source=0)
        Call H%get_ind_mult_r(ind_l,N_r,ind_r)
        !get the mode components for ind_r
        Call H%mode_r%get_mode_disc(lat,ind_r(:N_r),vec_r)
        !apply H on the right side only considering ind_l and ind_r
        allocate(vec_applied(size(ind_l)))
        Call H%mult_r_disc_disc(ind_r(:N_r),vec_r,ind_l,vec_applied)
        !get left mode for ind_l
        allocate(vec_l(size(ind_l)))
        Call H%mode_l%get_mode_exc_disc(lat,comp_l,ind_l,vec_l)
        !multiply left mode and applied right mode and reduce to wanted mode
        vec_applied=vec_l*vec_applied
        Call H%mode_l%mode_reduce_comp_disc(ind_l,vec_applied,comp_l,ind_vec,vec)
    end subroutine

    subroutine get_rN_single(this,H,lat,site,vec)
        class(t_deriv_r_N),intent(in)   :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        integer,intent(in)              :: site
        real(8),intent(inout)           :: vec(:)

        integer                         :: ind_vec(dim_modes_inner(this%order))

        integer,allocatable             :: ind_l(:)  !indices of considered mode on which H has to be applied 
        integer,allocatable             :: ind_r(:)  !indices of considered mode that have to be reduced
        integer                         :: N_l

        real(8),allocatable             :: vec_l(:)
        real(8),allocatable             :: vec_r(:)
        real(8),allocatable             :: vec_applied(:)
        integer     :: comp_r   !should this be saved in the type
        integer     :: i

        comp_r=H%mode_r%get_comp(this%order)
        ind_vec=[((site-1)*dim_modes_inner(this%order)+i,i=1,dim_modes_inner(this%order))]
        !get all vector indices on the right side that contribute in the reduction to ind_vec
        Call H%mode_r%get_ind_site(comp_r,site,ind_r)   
        !get all vector indices on the left  side that have contributions at ind_r after applying H
        N_l=size(ind_r)*10 !arbitrary, and think of something better with external work arrays
        allocate(ind_l(N_l),source=0)
        Call H%get_ind_mult_l(ind_r,N_l,ind_l)
        !get the mode components for ind_l
        Call H%mode_l%get_mode_disc(lat,ind_l(:N_l),vec_l)
        !apply H on the left side only considering ind_r and ind_l
        allocate(vec_applied(size(ind_r)))
        Call H%mult_l_disc_disc(ind_l(:N_l),vec_l,ind_r,vec_applied)
        !get left mode for ind_r
        allocate(vec_r(size(ind_r)))
        Call H%mode_r%get_mode_exc_disc(lat,comp_r,ind_r,vec_r)
        !multiply right mode and applied left mode and reduce to wanted mode
        vec_applied=vec_r*vec_applied
        Call H%mode_r%mode_reduce_comp_disc(ind_r,vec_applied,comp_r,ind_vec,vec)
    end subroutine
end module
