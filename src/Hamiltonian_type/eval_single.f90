module eval_single
use m_derived_types, only : lattice, number_different_order_parameters
use m_H_type, only :t_H_base
use m_work_ham_single, only: work_ham_single
implicit none


type t_eval_single
    integer :: comp=0
    integer :: order=0
    procedure(int_eval_single),pointer      :: calc => eval_single_unset
contains
    procedure   :: set
    procedure   :: copy
end type

abstract interface
    subroutine int_eval_single(this,H,E,i_m,lat,work)
        import t_H_base, lattice ,work_ham_single, t_eval_single
        class(t_eval_single),intent(in)     :: this
        class(t_H_base), intent(in)         :: H
        type(lattice), intent(in)           :: lat
        integer, intent(in)                 :: i_m
        ! output
        real(8), intent(out)                :: E
        !temporary data
        type(work_ham_single),intent(inout) ::  work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
    end subroutine 
end interface


contains

subroutine set(this,H,order)
    !set the component and calc pointer for a given order
    class(t_eval_single),intent(inout)  :: this
    class(t_H_base), intent(in)         :: H
    integer,intent(in)                  :: order

    integer         :: Nl_order,Nr_order !numbern of order parameters occurances at left/right side of Hamiltonian

    this%order=order

    Nl_order=count(H%op_l==order); Nr_order=count(H%op_r==order)
    if(Nl_order>0.and.Nr_order>0)then
        !in case of Hamiltonian implementation with only one side efficient for mult_l/r_disc-evaluation
        !it makes sense to check row_ind (or col_ind) as simple check to see if one mult_r_disc is implemented
        !if neither is implemented it should still set the eval_single (here by default eval_single_right) 
        !for reasonable error messages asking for implementation
        !direct crash does not make sense as it is not clear if eval_single is really going to be evaluationed
        if(H%row_max>0)then
            this%comp=H%mode_l%get_comp(order)
            this%calc=>eval_single_left
        else
            write(*,*) allocated(H%mode_r)
            this%comp=H%mode_r%get_comp(order)
            this%calc=>eval_single_right
        endif
    elseif(Nl_order>0)then
        this%comp=H%mode_l%get_comp(order)
        this%calc=>eval_single_left
    elseif(Nr_order>0)then
        this%comp=H%mode_r%get_comp(order)
        this%calc=>eval_single_right
    else
        this%comp=0
        this%calc=>eval_single_null
    endif
end subroutine

subroutine copy(this,eval_out)
    class(t_eval_single),intent(in   )  :: this
    class(t_eval_single),intent(inout)  :: eval_out

    eval_out%order =  this%order
    eval_out%comp  =  this%comp
    eval_out%calc  => this%calc
end subroutine

subroutine eval_single_left(this,H,E,site,lat,work)
    !evaluates the single energy for Hamiltonians were the left mode contains the sought order-parameter
    ! input
    class(t_eval_single),intent(in)     :: this
    class(t_H_base), intent(in)         :: H
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: site
    ! output
    real(8), intent(out)                :: E
    !temporary data
    type(work_ham_single),intent(inout) ::  work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
#define _dim_  H%dim_l_single(this%order)
#define _max_  H%col_max
#ifdef CPP_USE_WORK
    !temporary data slices
    integer,pointer,contiguous          :: ind_out (:)   !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
    integer,pointer,contiguous          :: ind_sum (:)   !ind_mult index where the a vec-entry start and end
    integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
    real(8),pointer,contiguous          :: vec_mult(:)   !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
    real(8),pointer,contiguous          :: mv_vec  (:)   !result vector matrix product without reduction
    real(8),pointer,contiguous          :: vec_oth (:)   !other mode
    !associate temporary arrays
    !!int vector slices
    ind_out (1:_dim_      )=>work%int_arr (1          :_dim_            )
    ind_sum (1:_dim_+1    )=>work%int_arr (1+_dim_    :_dim_* 2       +1)
    ind_mult(1:_dim_*_max_)=>work%int_arr (1+_dim_*2+1:_dim_*(2+_max_)+1)
    !!real vector slices
    mat_mult(1:_dim_*_max_)=>work%real_arr(1                  :_dim_*(  _max_  ))
    vec_mult(1:_dim_*_max_)=>work%real_arr(1+_dim_*(  _max_  ):_dim_*(2*_max_  ))
    mv_vec  (1:_dim_      )=>work%real_arr(1+_dim_*(2*_max_  ):_dim_*(2*_max_+1))
    vec_oth (1:_dim_      )=>work%real_arr(1+_dim_*(2*_max_+1):_dim_*(2*_max_+2))
#else
    !temporary arrays
    integer                     :: ind_out (_dim_      )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
    integer                     :: ind_sum (_dim_+1    )     !ind_mult index where the a vec-entry start and end
    integer                     :: ind_mult(_dim_*_max_)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8)                     :: mat_mult(_dim_*_max_)     !matrix entries corresponding to ind_mult
    real(8)                     :: vec_mult(_dim_*_max_)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)

    real(8)                     :: mv_vec  (_dim_)           !result vector matrix product without reduction
    real(8)                     :: vec_oth (_dim_)           !other mode 
#endif
    !get indices of left mode which has components of site site at component comp
    Call H%mode_l%get_ind_site(this%comp, site, _dim_, ind_out)  
    !get Hamiltonian.right_mode product for all left mode indices ind_out
    Call H%mult_r_disc(site, lat, _dim_, ind_out, mv_vec, ind_sum, ind_Mult, mat_mult, vec_mult)
    !get left mode indices ind_out excluding mode number comp
    Call H%mode_l%get_mode_disc(lat, _dim_, ind_out, vec_oth)
    !energy is sum of left mode times Hamiltonian.right_mode values
    E=sum(mv_vec*vec_oth)
#undef _max_
#undef _dim_
end subroutine

subroutine eval_single_right(this,H,E,site,lat,work)
    !evaluates the single energy for Hamiltonians were the right mode contains the sought order-parameter
    ! input
    class(t_eval_single),intent(in)     :: this
    class(t_H_base), intent(in)         :: H
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: site
    ! output
    real(8), intent(out)                :: E
    !temporary data
    type(work_ham_single),intent(inout) ::  work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations
#define _dim_  H%dim_r_single(this%order)
#define _max_  H%row_max
#ifdef CPP_USE_WORK
    !temporary data slices
    integer,pointer,contiguous          :: ind_out (:)   !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
    integer,pointer,contiguous          :: ind_sum (:)   !ind_mult index where the a vec-entry start and end
    integer,pointer,contiguous          :: ind_mult(:)   !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8),pointer,contiguous          :: mat_mult(:)   !matrix entries corresponding to ind_mult
    real(8),pointer,contiguous          :: vec_mult(:)   !values of discontiguous left mode which has to be evaluated (indices of ind_mult)
    real(8),pointer,contiguous          :: mv_vec  (:)   !result vector matrix product without reduction
    real(8),pointer,contiguous          :: vec_oth (:)   !other mode
    !associate temporary arrays
    !!int vector slices
    ind_out (1:_dim_      )=>work%int_arr (1          :_dim_            )
    ind_sum (1:_dim_+1    )=>work%int_arr (1+_dim_    :_dim_* 2       +1)
    ind_mult(1:_dim_*_max_)=>work%int_arr (1+_dim_*2+1:_dim_*(2+_max_)+1)
    !!real vector slices
    mat_mult(1:_dim_*_max_)=>work%real_arr(1                  :_dim_*(  _max_  ))
    vec_mult(1:_dim_*_max_)=>work%real_arr(1+_dim_*(  _max_  ):_dim_*(2*_max_  ))
    mv_vec  (1:_dim_      )=>work%real_arr(1+_dim_*(2*_max_  ):_dim_*(2*_max_+1))
    vec_oth (1:_dim_      )=>work%real_arr(1+_dim_*(2*_max_+1):_dim_*(2*_max_+2))
#else
    !temporary arrays
    integer         :: ind_out (_dim_      )     !indices of all right mode entries which contain the order paramtere order of the site corresponding to i_m
    integer         :: ind_sum (_dim_+1    )     !ind_mult index where the a vec-entry start and end
    integer         :: ind_mult(_dim_*_max_)     !indices of the left array which have non-vanishing contributions to get the ind entries of the vec/mat product
    real(8)         :: mat_mult(_dim_*_max_)     !matrix entries corresponding to ind_mult
    real(8)         :: vec_mult(_dim_*_max_)     !values of discontiguous left mode which has to be evaluated (indices of ind_mult)

    real(8)         :: mv_vec  (_dim_)           !result vector matrix product without reduction
    real(8)         :: vec_oth (_dim_)           !other mode 
#endif
    !get indices of left mode which has components of site site at component comp
    Call H%mode_r%get_ind_site(this%comp, site, _dim_, ind_out)  
    !get Hamiltonian.right_mode product for all left mode indices ind_out
    Call H%mult_l_disc(site, lat, _dim_, ind_out, mv_vec, ind_sum, ind_Mult, mat_mult, vec_mult)
    !get left mode indices ind_out excluding mode number comp
    Call H%mode_r%get_mode_disc(lat, _dim_, ind_out, vec_oth)
    !energy is sum of left mode times Hamiltonian.right_mode values
    E=sum(mv_vec*vec_oth)
#undef _max_
#undef _dim_
    end subroutine

subroutine eval_single_null(this,H,E,site,lat,work)
    !In case the Hamiltonian has no entry in the considered space eval_single shall return zero energy
    class(t_eval_single),intent(in)     :: this
    class(t_H_base), intent(in)         :: H
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: site
    ! output
    real(8), intent(out)                :: E
    !temporary data
    type(work_ham_single),intent(inout) ::  work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations

    E=0.0d0
end subroutine

subroutine eval_single_unset(this,H,E,site,lat,work)
    class(t_eval_single),intent(in)     :: this
    class(t_H_base), intent(in)         :: H
    type(lattice), intent(in)           :: lat
    integer, intent(in)                 :: site
    ! output
    real(8), intent(out)                :: E
    !temporary data
    type(work_ham_single),intent(inout) ::  work    !data type containing the temporary data for this calculation to prevent constant allocations/deallocations

    ERROR STOP "The eval_set pointer of t_eval_single has not been set, programming mistake?"
end subroutine


end module
