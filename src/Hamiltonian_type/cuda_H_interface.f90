module m_cuda_H_interface
    implicit none
    public
    interface
        subroutine cuda_create_handle(handle) bind( c, name="cuda_create_handle" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)                :: handle
        end subroutine

        subroutine cuda_H_init(nnz,dimH,row,col,val,mat,handle) bind( c, name="cuda_H_init" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),value           :: nnz
           integer( kind = c_int ),intent(in)      :: dimH(*)
           integer( kind = c_int ),intent(in)      :: row(*),col(*)
           real( kind = c_double ),intent(in)      :: val(*)
           type(C_PTR),intent(inout)               :: mat
           type(C_PTR),intent(in)                  :: handle
        end subroutine

        subroutine cuda_H_mult_mat_vec(mat,in_vec,out_vec,buffer,handle) bind( c, name="cuda_H_mult_mat_vec" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           type(C_PTR),intent(in)                  :: in_vec
           type(C_PTR),intent(in)                  :: out_vec
           type(C_PTR),intent(in)                  :: buffer
           type(C_PTR),intent(in)                  :: handle
        end subroutine

        subroutine cuda_H_copy(H_in,H_copy) bind( c, name="cuda_H_copy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)      		:: H_in
           type(C_PTR),intent(inout)   		:: H_copy
        end subroutine

        subroutine cuda_H_add(mat_1,mat_2,mat_s,handle) bind( c, name="cuda_H_add" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)      		:: mat_1
           type(C_PTR),intent(in)         	:: mat_2
           type(C_PTR),intent(inout)      	:: mat_s
           type(C_PTR),intent(in)			:: handle
        end subroutine

        subroutine cuda_H_destroy(mat) bind( c, name="cuda_H_destroy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)		:: mat
        end subroutine

        subroutine cuda_set_buffer(buffer,mat_1,mat_2,mat_s,handle) bind( c, name="cuda_set_buffer" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)   		:: buffer
           type(C_PTR),intent(in)      		:: mat_1
           type(C_PTR),intent(in)         	:: mat_2
           type(C_PTR),intent(inout)      	:: mat_s
           type(C_PTR),intent(in)			:: handle
        end subroutine

        subroutine cuda_free_buffer(buffer) bind( c, name="cuda_free_buffer" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)   		:: buffer
        end subroutine
    end interface
end module 
