module m_deriv_rank2
use m_deriv_base
use m_derived_types,only : lattice
use m_H_type, only: t_H_base

private
public t_deriv_l_1, t_deriv_l_1_sym, t_deriv_r_1 

type,extends(t_deriv_l)   :: t_deriv_l_1
contains
    procedure   :: get => get_l1
    procedure   :: get_single => get_l1_single
end type

type,extends(t_deriv_l)   :: t_deriv_l_1_sym
contains
    procedure   :: get => get_l1_sym
    procedure   :: get_single => get_l1_sym_single
end type

type,extends(t_deriv_r)   :: t_deriv_r_1
contains
    procedure   :: get => get_r1
    procedure   :: get_single => get_r1_single
end type

contains
    subroutine get_l1(this,H,lat,vec)
        class(t_deriv_l_1),intent(in)   :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_r(lat,vec)
    end subroutine

    subroutine get_r1(this,H,lat,vec)
        class(t_deriv_r_1),intent(in)     :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_l(lat,vec)
    end subroutine

    subroutine get_l1_sym(this,H,lat,vec)
        class(t_deriv_l_1_sym),intent(in)     :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        Call H%mult_r(lat,vec)
        vec=vec*2.0d0
    end subroutine

    subroutine get_l1_single(this,H,lat,site,vec)
        class(t_deriv_l_1),intent(in)   :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        integer,intent(in)              :: site
        real(8),intent(inout)           :: vec(:)
        integer     :: ind(size(vec))
        integer     :: i

        ind=[((site-1)*size(vec)+i,i=1,size(vec))]
        Call H%mult_r_ind(lat,size(vec),ind,vec)
    end subroutine

    subroutine get_r1_single(this,H,lat,site,vec)
        class(t_deriv_r_1),intent(in)     :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        integer,intent(in)              :: site
        real(8),intent(inout)           :: vec(:)
        integer     :: ind(size(vec))
        integer     :: i

        ind=[((site-1)*size(vec)+i,i=1,size(vec))]
        Call H%mult_l_ind(lat,size(vec),ind,vec)
    end subroutine

    subroutine get_l1_sym_single(this,H,lat,site,vec)
        use m_type_lattice, only: dim_modes_inner
        class(t_deriv_l_1_sym),intent(in)     :: this
        class(t_H_base),intent(in)           :: H
        type(lattice),intent(in)        :: lat
        integer,intent(in)              :: site
        real(8),intent(inout)           :: vec(:)
        integer     :: ind(size(vec))
        integer     :: i

        ind=[((site-1)*size(vec)+i,i=1,size(vec))]
        Call H%mult_r_ind(lat,size(vec),ind,vec)
        vec=vec*2.0d0
    end subroutine

end module
