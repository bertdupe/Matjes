module m_deriv_rank2
use m_deriv_base
use m_derived_types,only : lattice
use m_H_type, only: t_H_base
use m_work_ham_single, only:  work_ham_single, work_mode

integer,parameter       ::  comp=1  !mode component always 1 for rank2 derivative
private
public get_l1,get_r1,get_l1_sym, get_r1_sym
public get_l1_single,get_r1_single, get_l1_sym_single

contains
    subroutine get_l1(this,H,lat,vec,work)
        class(t_deriv),intent(in)       :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)
        type(work_mode),intent(inout)   :: work
        real(8),parameter               :: alpha=1.0d0
        real(8),parameter               :: beta=1.0d0

        Call H%mult_r(lat,vec,work,alpha,beta)
    end subroutine

    subroutine get_r1(this,H,lat,vec,work)
        class(t_deriv),intent(in)       :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        type(work_mode),intent(inout)   :: work
        real(8),intent(inout)           :: vec(:)
        real(8),parameter               :: alpha=1.0d0
        real(8),parameter               :: beta=1.0d0

        Call H%mult_l(lat,vec,work,alpha,beta)
    end subroutine

    subroutine get_l1_sym(this,H,lat,vec,work)
        class(t_deriv),intent(in)       :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        type(work_mode),intent(inout)   :: work
        real(8),intent(inout)           :: vec(:)
        real(8),parameter               :: alpha=2.0d0
        real(8),parameter               :: beta=1.0d0

        Call H%mult_r(lat,vec,work,alpha,beta)
    end subroutine

    subroutine get_r1_sym(this,H,lat,vec,work)
        class(t_deriv),intent(in)       :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        type(work_mode),intent(inout)   :: work
        real(8),intent(inout)           :: vec(:)
        real(8),parameter               :: alpha=2.0d0
        real(8),parameter               :: beta=1.0d0

        Call H%mult_l(lat,vec,work,alpha,beta)
    end subroutine


    subroutine get_r1_single(this,H,lat,site,work,vec)
        class(t_deriv),intent(in)           :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)

#define _dim_  H%dim_r_single(H%op_r(1))
#define _max_  H%col_max
#ifdef CPP_USE_WORK
        !temporary data slices
        integer,pointer,contiguous          :: ind_out(:)    !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer,pointer,contiguous          :: ind_sum(:)    !ind_mult index where the a vec-entry start and end
        integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
        real(8),pointer,contiguous          :: vec_mult(:)      !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
        !associate temporary arrays
        !!int vector slices
        ind_out (1:_dim_      )=>work%int_arr (1          :_dim_            )
        ind_sum (1:_dim_+1    )=>work%int_arr (1+_dim_    :_dim_* 2       +1)
        ind_mult(1:_dim_*_max_)=>work%int_arr (1+_dim_*2+1:_dim_*(2+_max_)+1)
        !!real vector slices
        mat_mult(1:_dim_*_max_)=>work%real_arr(1            :_dim_*_max_  )
        vec_mult(1:_dim_*_max_)=>work%real_arr(1+_dim_*_max_:_dim_*_max_*2)
#else
        !temporary arrays
        integer                     :: ind_out (_dim_      )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer                     :: ind_sum (_dim_+1    )     !ind_mult index where the a vec-entry start and end
        integer                     :: ind_mult(_dim_*_max_)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8)                     :: mat_mult(_dim_*_max_)     !matrix entries corresponding to ind_mult
        real(8)                     :: vec_mult(_dim_*_max_)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
#endif

        !get indices of the output vector 
        Call H%mode_r%get_ind_site(comp,site,_dim_,ind_out)  
        Call H%mult_l_disc(lat,_dim_,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
        write(*,*) 'Bertrand', vec
#undef _max_
#undef _dim_
    end subroutine

    subroutine get_l1_single(this,H,lat,site,work,vec)
        class(t_deriv),intent(in)           :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)
        integer     :: i

#define _dim_  H%dim_l_single(H%op_l(1))
#define _max_  H%row_max
#ifdef CPP_USE_WORK
        !temporary data slices
        integer,pointer,contiguous          :: ind_out(:)    !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer,pointer,contiguous          :: ind_sum(:)    !ind_mult index where the a vec-entry start and end
        integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
        real(8),pointer,contiguous          :: vec_mult(:)      !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
        !associate temporary arrays
        !!int vector slices
        ind_out (1:_dim_      )=>work%int_arr (1          :_dim_            )
        ind_sum (1:_dim_+1    )=>work%int_arr (1+_dim_    :_dim_* 2       +1)
        ind_mult(1:_dim_*_max_)=>work%int_arr (1+_dim_*2+1:_dim_*(2+_max_)+1)
        !!real vector slices
        mat_mult(1:_dim_*_max_)=>work%real_arr(1            :_dim_*_max_  )
        vec_mult(1:_dim_*_max_)=>work%real_arr(1+_dim_*_max_:_dim_*_max_*2)
#else
        !temporary arrays
        integer                     :: ind_out (_dim_      )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer                     :: ind_sum (_dim_+1    )     !ind_mult index where the a vec-entry start and end
        integer                     :: ind_mult(_dim_*_max_)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8)                     :: mat_mult(_dim_*_max_)     !matrix entries corresponding to ind_mult
        real(8)                     :: vec_mult(_dim_*_max_)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
#endif

        !get indices of the output vector 
        Call H%mode_l%get_ind_site(comp,site,_dim_,ind_out)  
        Call H%mult_r_disc(lat,_dim_,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
#undef _max_
#undef _dim_
    end subroutine

    subroutine get_l1_sym_single(this,H,lat,site,work,vec)
        use m_type_lattice, only: dim_modes_inner
        class(t_deriv),intent(in)           :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)
#define _dim_  H%dim_l_single(H%op_l(1))
#define _max_  H%row_max
#ifdef CPP_USE_WORK
        !temporary data slices
        integer,pointer,contiguous          :: ind_out(:)    !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer,pointer,contiguous          :: ind_sum(:)    !ind_mult index where the a vec-entry start and end
        integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
        real(8),pointer,contiguous          :: vec_mult(:)      !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
        !associate temporary arrays
        !!int vector slices
        ind_out (1:_dim_      )=>work%int_arr (1          :_dim_            )
        ind_sum (1:_dim_+1    )=>work%int_arr (1+_dim_    :_dim_* 2       +1)
        ind_mult(1:_dim_*_max_)=>work%int_arr (1+_dim_*2+1:_dim_*(2+_max_)+1)
        !!real vector slices
        mat_mult(1:_dim_*_max_)=>work%real_arr(1            :_dim_*_max_  )
        vec_mult(1:_dim_*_max_)=>work%real_arr(1+_dim_*_max_:_dim_*_max_*2)
#else
        !temporary arrays
        integer                     :: ind_out (_dim_      )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer                     :: ind_sum (_dim_+1    )     !ind_mult index where the a vec-entry start and end
        integer                     :: ind_mult(_dim_*_max_)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8)                     :: mat_mult(_dim_*_max_)     !matrix entries corresponding to ind_mult
        real(8)                     :: vec_mult(_dim_*_max_)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
#endif

        Call H%mode_l%get_ind_site(comp,site,_dim_,ind_out)
        Call H%mult_r_disc(lat,_dim_,ind_out,vec,ind_sum,ind_Mult,mat_mult,vec_mult)
        vec=vec*2.0d0
#undef _max_
#undef _dim_
    end subroutine
end module
