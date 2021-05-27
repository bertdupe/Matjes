module m_deriv_rankN
use m_derived_types,only : lattice, dim_modes_inner
use m_deriv_base
use m_H_type, only: t_H_base
use m_work_ham_single, only:  work_ham_single

private
public get_lN, get_rN, get_lN_single, get_rN_single

contains
    subroutine get_lN(this,H,lat,vec)
        class(t_deriv),intent(in)       :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        real(8)                         :: tmp(H%dimH(1))   !multipied, but not reduced
        real(8),parameter               :: beta=1.0d0
    
        Call H%mult_r(lat,tmp)
        Call H%mode_l%reduce_other_exc(lat,this%order,tmp,beta,vec)
    end subroutine

    subroutine get_rN(this,H,lat,vec)
        class(t_deriv),intent(in)       :: this
        class(t_H_base),intent(in)      :: H
        type(lattice),intent(in)        :: lat
        real(8),intent(inout)           :: vec(:)

        real(8)                         :: tmp(H%dimH(2))   !multipied, but not reduced
        real(8),parameter               :: beta=1.0d0
    
        Call H%mult_l(lat,tmp)
        Call H%mode_r%reduce_other_exc(lat,this%order,tmp,beta,vec)
    end subroutine

    subroutine get_lN_single(this,H,lat,site,work,vec)
        class(t_deriv),intent(in)           :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)

        integer     :: comp   !should this be saved in the type
#define _dim_  H%dim_l_single(this%order)
#define _max_  H%col_max
#ifdef CPP_USE_WORK
        !temporary data slices
        integer,pointer,contiguous          :: ind_out(:)    !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer,pointer,contiguous          :: ind_sum(:)    !ind_mult index where the a vec-entry start and end
        integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
        real(8),pointer,contiguous          :: vec_mult(:)   !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
        real(8),pointer,contiguous          :: mv_vec  (:)   !result vector matrix product without reduction
        real(8),pointer,contiguous          :: vec_red (:)   !right mode without reduction
        !associate temporary arrays
        !!int vector slices
        ind_out (1:_dim_      )=>work%int_arr (1          :_dim_            )
        ind_sum (1:_dim_+1    )=>work%int_arr (1+_dim_    :_dim_* 2       +1)
        ind_mult(1:_dim_*_max_)=>work%int_arr (1+_dim_*2+1:_dim_*(2+_max_)+1)
        !!real vector slices
        mat_mult(1:_dim_*_max_)=>work%real_arr(1                  :_dim_*_max_  )
        vec_mult(1:_dim_*_max_)=>work%real_arr(1+_dim_*(  _max_  ):_dim_*_max_*2)
        mv_vec  (1:_dim_      )=>work%real_arr(1+_dim_*(2*_max_  ):_dim_*(_max_*2+1))
        vec_red (1:_dim_      )=>work%real_arr(1+_dim_*(2*_max_+1):_dim_*(_max_*2+2))
#else
        !temporary arrays
        integer                     :: ind_out (_dim_      )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer                     :: ind_sum (_dim_+1    )     !ind_mult index where the a vec-entry start and end
        integer                     :: ind_mult(_dim_*_max_)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8)                     :: mat_mult(_dim_*_max_)     !matrix entries corresponding to ind_mult
        real(8)                     :: vec_mult(_dim_*_max_)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)

        real(8)                     :: mv_vec  (_dim_)           !result vector matrix product without reduction
        real(8)                     :: vec_red (_dim_)           !right mode without reduction
#endif

        comp=H%mode_l%get_comp(this%order)
        !get indices of left mode which has components of site site at component comp
        Call H%mode_l%get_ind_site(comp,site,_dim_,ind_out)  
        !get Hamiltonian.right_mode product for all left mode indices ind_out
        Call H%mult_r_disc(site,lat,_dim_,ind_out,mv_vec,ind_sum,ind_Mult,mat_mult,vec_mult)
        !get left mode indices ind_out excluding mode number comp
        Call H%mode_l%get_mode_exc_disc(lat,comp,_dim_,ind_out,vec_red)
        !mutliply Hamiltonian.right_mode values with left mode exluding comp values
        mv_vec=mv_vec*vec_red
        !reduce resulting vector according to rules of left mode
        Call H%mode_l%reduce_site_vec(comp,mv_vec,vec)
#undef _max_
#undef _dim_
    end subroutine

    subroutine get_rN_single(this,H,lat,site,work,vec)
        class(t_deriv),intent(in)           :: this
        class(t_H_base),intent(in)          :: H
        type(lattice),intent(in)            :: lat
        integer,intent(in)                  :: site
        type(work_ham_single),intent(inout) :: work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
        real(8),intent(inout)               :: vec(:)

        integer     :: comp   !should this be saved in the type
#define _dim_  H%dim_r_single(this%order)
#define _max_  H%row_max
#ifdef CPP_USE_WORK
        !temporary data slices
        integer,pointer,contiguous          :: ind_out(:)    !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer,pointer,contiguous          :: ind_sum(:)    !ind_mult index where the a vec-entry start and end
        integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
        real(8),pointer,contiguous          :: vec_mult(:)   !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
        real(8),pointer,contiguous          :: mv_vec  (:)   !result vector matrix product without reduction
        real(8),pointer,contiguous          :: vec_red (:)   !right mode without reduction
        !associate temporary arrays
        !!int vector slices
        ind_out (1:_dim_      )=>work%int_arr (1          :_dim_            )
        ind_sum (1:_dim_+1    )=>work%int_arr (1+_dim_    :_dim_* 2       +1)
        ind_mult(1:_dim_*_max_)=>work%int_arr (1+_dim_*2+1:_dim_*(2+_max_)+1)
        !!real vector slices
        mat_mult(1:_dim_*_max_)=>work%real_arr(1                  :_dim_*_max_  )
        vec_mult(1:_dim_*_max_)=>work%real_arr(1+_dim_*(  _max_  ):_dim_*_max_*2)
        mv_vec  (1:_dim_      )=>work%real_arr(1+_dim_*(2*_max_  ):_dim_*(_max_*2+1))
        vec_red (1:_dim_      )=>work%real_arr(1+_dim_*(2*_max_+1):_dim_*(_max_*2+2))
#else
        !temporary arrays
        integer                     :: ind_out (_dim_      )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
        integer                     :: ind_sum (_dim_+1    )     !ind_mult index where the a vec-entry start and end
        integer                     :: ind_mult(_dim_*_max_)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
        real(8)                     :: mat_mult(_dim_*_max_)     !matrix entries corresponding to ind_mult
        real(8)                     :: vec_mult(_dim_*_max_)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)

        real(8)                     :: mv_vec  (_dim_)           !result vector matrix product without reduction
        real(8)                     :: vec_red (_dim_)           !right mode without reduction
#endif
        comp=H%mode_r%get_comp(this%order)
        !get indices of right mode which has components of site site at component comp
        Call H%mode_r%get_ind_site(comp,site,_dim_,ind_out)
        !get left_mode.Hamiltonian product for all right mode indices ind_out
        Call H%mult_l_disc(site,lat,_dim_,ind_out,mv_vec,ind_sum,ind_Mult,mat_mult,vec_mult)
        !get right mode indices ind_out excluding mode number comp
        Call H%mode_r%get_mode_exc_disc(lat,comp,_dim_,ind_out,vec_red)
        !mutliply left_mode.Hamiltonian values with right mode exluding comp values
        mv_vec=mv_vec*vec_red
        !reduce resulting vector according to rules of right mode
        Call H%mode_r%reduce_site_vec(comp,mv_vec,vec)
#undef _max_
#undef _dim_
    end subroutine
end module
