module m_cuda_H_vec_interface
    implicit none
    public
    interface
        subroutine cuda_fvec_init(vec,size) bind( c, name="cuda_fvec_init" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)                :: vec
           integer,value                            :: size
        end subroutine

        subroutine cuda_fvec_destroy(vec) bind( c, name="cuda_fvec_destroy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)                :: vec
        end subroutine

        subroutine cuda_fvec_set(vec,arr_in) bind( c, name="cuda_fvec_set" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                   :: vec  !actually the value gets changed, but not the pointer (VOLATILE?)
           real( kind = c_double ),intent(in)       :: arr_in(*)
        end subroutine

        subroutine cuda_fvec_get(vec,arr_out) bind( c, name="cuda_fvec_get" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                   :: vec
           real( kind = c_double ),intent(inout)    :: arr_out(*)
        end subroutine

        subroutine cuda_fvec_alloccopy(vec_in,vec_out) bind( c, name="cuda_fvec_alloccopy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                   :: vec_in
           type(C_PTR),intent(inout)                :: vec_out
        end subroutine
    end interface
end module 
