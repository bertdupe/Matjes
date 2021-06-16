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

        subroutine cuda_H_mult_r(mat,in_vec,out_vec,alpha,beta,buffer,handle) bind( c, name="cuda_H_mult_r" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           type(C_PTR),intent(in)                  :: in_vec
           type(C_PTR),intent(in)                  :: out_vec
           real( kind = c_double ),intent(in)      :: alpha
           real( kind = c_double ),intent(in)      :: beta
           type(C_PTR),intent(in)                  :: buffer
           type(C_PTR),intent(in)                  :: handle
        end subroutine

        subroutine cuda_H_mult_l(mat,in_vec,out_vec,alpha,beta,buffer,handle) bind( c, name="cuda_H_mult_l" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           type(C_PTR),intent(in)                  :: in_vec
           type(C_PTR),intent(in)                  :: out_vec
           real( kind = c_double ),intent(in)      :: alpha
           real( kind = c_double ),intent(in)      :: beta
           type(C_PTR),intent(in)                  :: buffer
           type(C_PTR),intent(in)                  :: handle
        end subroutine

        subroutine cuda_H_copy(H_in,H_copy) bind( c, name="cuda_H_copy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)           :: H_in
           type(C_PTR),intent(inout)        :: H_copy
        end subroutine

        subroutine cuda_H_add(mat_1,mat_2,mat_s,handle) bind( c, name="cuda_H_add" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)           :: mat_1
           type(C_PTR),intent(in)           :: mat_2
           type(C_PTR),intent(inout)        :: mat_s
           type(C_PTR),intent(in)           :: handle
        end subroutine

        subroutine cuda_H_destroy(mat) bind( c, name="cuda_H_destroy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)        :: mat
        end subroutine

        subroutine cuda_set_buffer(buffer,mat,transp,in_vec,out_vec,handle) bind( c, name="cuda_set_buffer" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)        :: buffer
           type(C_PTR),intent(in)           :: mat
           logical(C_BOOL),value            :: transp
           type(C_PTR),intent(in)           :: in_vec
           type(C_PTR),intent(in)           :: out_vec
           type(C_PTR),intent(in)           :: handle
        end subroutine

        subroutine cuda_free_buffer(buffer) bind( c, name="cuda_free_buffer" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)        :: buffer
        end subroutine

!        subroutine cuda_H_mult_mat_disc_disc(mat,vec_in_dev,vec_out_dev,N_in,ind_in,vec_in,N_out,ind_out,vec_out) bind( c, name="cuda_H_mult_mat_disc_disc" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           type(C_PTR),intent(in)                  :: vec_in_dev
!           type(C_PTR),intent(in)                  :: vec_out_dev
!           integer(c_int),value                    :: N_in
!           integer(c_int),intent(in)               :: ind_in(*)
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           integer(c_int),intent(inout)            :: N_out
!           integer(c_int),intent(inout)            :: ind_out(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
    end interface
end module 
