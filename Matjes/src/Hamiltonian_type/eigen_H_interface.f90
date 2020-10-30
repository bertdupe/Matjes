module m_eigen_H_interface
    use,intrinsic :: iso_c_binding
    implicit none
    public
    interface
        subroutine eigen_H_init(nnz,dimH,row,col,val,H) bind( c, name="eigen_H_init" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),value           :: nnz
           integer( kind = c_int ),intent(in)      :: dimH(*)
           integer( kind = c_int ),intent(in)      :: row(*),col(*)
           real( kind = c_double ),intent(in)      :: val(*)
           type(C_PTR),intent(inout)               :: H
        end subroutine
        
        subroutine eigen_H_mult_mat_vec(mat,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine
        
        subroutine eigen_H_mult_vec_mat(mat,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_eval_single(ind,dim_mode,vec_l,vec_r,H,E) bind( c, name="eigen_H_eval_single" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),value           :: ind,dim_mode
           real( kind = c_double ),intent(in)      :: vec_l(*),vec_r(*)
           type(C_PTR),intent(in)                  :: H
           real( kind = c_double ),intent(out)     :: E
        end subroutine
        
        subroutine eigen_H_copy(H_in,H_out) bind( c, name="eigen_H_copy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)         :: H_in
           type(C_PTR),intent(inout)      :: H_out
        end subroutine
        
        subroutine eigen_H_add(H_sum,H_add) bind( c, name="eigen_H_add" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)      :: H_sum
           type(C_PTR),intent(in)         :: H_add
        end subroutine
        
        subroutine eigen_H_destroy(H) bind( c, name="eigen_H_destroy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)         :: H
        end subroutine
    
    end interface
    
end module 
