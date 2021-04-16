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


        subroutine eigen_H_bcast(id,mas,ismas,mat,comm) bind( c, name="eigen_H_bcast" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),value            :: id    !MPI-rank of the processor
           integer( kind = c_int ),value            :: mas   !MPI-rank of the master
           logical( kind = c_bool ),value           :: ismas !is master to bcast from
           type(C_PTR),intent(inout)                :: mat   !matrix to broadcast
           integer( kind = c_int ),intent(in)       :: comm  !MPI-communicator  the processor
        end subroutine
        
        subroutine eigen_H_mult_mat_vec(mat,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_vec_mat_cont(mat,bnd_min,bnd_max,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat_cont" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           integer(c_int),value                    :: bnd_min,bnd_max               
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_vec_mat_disc(mat,N,ind,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat_disc" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           integer(c_int),value                    :: N 
           integer(c_int),intent(in)               :: ind(*)
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_mat_vec_cont(mat,bnd_min,bnd_max,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec_cont" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           integer(c_int),value                    :: bnd_min,bnd_max               
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_mat_vec_disc(mat,N,ind,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec_disc" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           integer(c_int),value                    :: N 
           integer(c_int),intent(in)               :: ind(*)
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_mat_disc_disc(mat,N_in,ind_in,vec_in,N_out,ind_out,vec_out) bind( c, name="eigen_H_mult_mat_disc_disc" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           integer(c_int),value                    :: N_in
           integer(c_int),intent(in)               :: ind_in(*)
           real( kind = c_double ),intent(in)      :: vec_in(*)
           integer(c_int),intent(out)              :: N_out
           integer(c_int),intent(inout)            :: ind_out(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine


        subroutine eigen_H_mult_mat_vec_single(mat,bnd_min,bnd_max,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec_single" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           integer(c_int),value                    :: bnd_min,bnd_max               
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine
        
        subroutine eigen_H_mult_vec_mat(mat,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_H_mult_vec_mat_single(mat,bnd_min,bnd_max,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat_single" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: mat
           integer(c_int),value                    :: bnd_min,bnd_max               
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

        subroutine eigen_H_get_ind_mult_r(H,N_in,ind_in,N_out,ind_out) bind(c, name="eigen_H_get_ind_mult_r")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: H
           integer( kind = c_int ),value           :: N_in
           integer( kind = c_int ),intent(in)      :: ind_in(*)
           integer( kind = c_int ),intent(inout)   :: N_out
           integer( kind = c_int ),intent(inout)   :: ind_out(*)
        end subroutine

        subroutine eigen_H_get_ind_mult_l(H,N_in,ind_in,N_out,ind_out) bind(c, name="eigen_H_get_ind_mult_l")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: H
           integer( kind = c_int ),value           :: N_in
           integer( kind = c_int ),intent(in)      :: ind_in(*)
           integer( kind = c_int ),intent(inout)   :: N_out
           integer( kind = c_int ),intent(inout)   :: ind_out(*)
        end subroutine

        subroutine eigen_mult_r_disc_disc(H,N_r,ind_r,vec_r,N_l,ind_l,vec_out) bind(c, name="eigen_mult_r_disc_disc")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: H
           integer( kind = c_int ),value           :: N_r
           integer( kind = c_int ),intent(in)      :: ind_r(*)
           real( kind = c_double ),intent(in)      :: vec_r(*)
           integer( kind = c_int ),value           :: N_l
           integer( kind = c_int ),intent(in)      :: ind_l(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_mult_l_disc_disc(H,N_l,ind_l,vec_l,N_r,ind_r,vec_out) bind(c, name="eigen_mult_l_disc_disc")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: H
           integer( kind = c_int ),value           :: N_l
           integer( kind = c_int ),intent(in)      :: ind_l(*)
           real( kind = c_double ),intent(in)      :: vec_l(*)
           integer( kind = c_int ),value           :: N_r
           integer( kind = c_int ),intent(in)      :: ind_r(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine
    
    end interface
    
end module 
