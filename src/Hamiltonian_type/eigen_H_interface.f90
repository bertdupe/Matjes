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

        subroutine eigen_get_transpose(H,H_T) bind( c, name="eigen_get_transpose")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: H
           type(C_PTR),intent(inout)               :: H_T
        end subroutine
         
        subroutine eigen_H_mult_r(mat,vec_in,vec_out) bind( c, name="eigen_H_mult_r" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_l(mat,vec_in,vec_out) bind( c, name="eigen_H_mult_l" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_r_ind(mat,vec,N,ind_out,vec_out) bind(c,name="eigen_H_mult_r_ind")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           real( kind = c_double ),intent(in)      :: vec(*)
           integer(c_int),intent(in)               :: N
           integer(c_int),intent(in)               :: ind_out(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_l_ind(mat,vec,N,ind_out,vec_out) bind(c,name="eigen_H_mult_l_ind")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           real( kind = c_double ),intent(in)      :: vec(*)
           integer(c_int),intent(in)               :: N
           integer(c_int),intent(in)               :: ind_out(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
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

        subroutine eigen_get_dat(h,nnz,dimH,col,row,val) bind(c, name='eigen_get_dat')
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)               :: H
           integer( kind = c_int ),intent(inout)   :: nnz
           integer( kind = c_int ),intent(inout)   :: dimH(2)
           type(C_PTR),intent(inout)               :: col
           type(C_PTR),intent(inout)               :: row
           type(C_PTR),intent(inout)               :: val
        end subroutine

#ifdef CPP_MPI
        subroutine eigen_H_send(id,tag,mat,comm) bind( c, name="eigen_H_send" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),intent(in)       :: id    !MPI-rank to target
           integer( kind = c_int ),intent(in)       :: tag   !MPI-tag
           type(C_PTR),intent(in)                   :: mat   !matrix to send
           integer( kind = c_int ),intent(in)       :: comm  !MPI-communicator
        end subroutine

        subroutine eigen_H_recv(id,tag,mat,comm) bind( c, name="eigen_H_recv" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),intent(in)       :: id    !MPI-rank to target
           integer( kind = c_int ),intent(in)       :: tag   !MPI-tag
           type(C_PTR),intent(inout)                :: mat   !matrix to send
           integer( kind = c_int ),intent(in)       :: comm  !MPI-communicator
        end subroutine

        subroutine eigen_H_bcast(id,mas,ismas,mat,comm) bind( c, name="eigen_H_bcast" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),intent(in)       :: id    !MPI-rank of the processor
           integer( kind = c_int ),intent(in)       :: mas   !MPI-rank of the master
           logical( kind = c_bool ),intent(in)      :: ismas !is master to bcast from
           type(C_PTR),intent(inout)                :: mat   !matrix to broadcast
           integer( kind = c_int ),intent(in)       :: comm  !MPI-communicator  the processor
        end subroutine

        subroutine eigen_H_distribute(id,mas,ismas,mat,comm) bind( c, name="eigen_H_distribute" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),intent(in)       :: id    !MPI-rank of the processor
           integer( kind = c_int ),intent(in)       :: mas   !MPI-rank of the master
           logical( kind = c_bool ),intent(in)      :: ismas !is master to bcast from
           type(C_PTR),intent(inout)                :: mat   !matrix to broadcast
           integer( kind = c_int ),intent(in)       :: comm  !MPI-communicator  the processor
        end subroutine
#endif
    end interface
    
end module 
