module m_eigen_mode_interface
    use,intrinsic :: iso_c_binding
    implicit none
    public
    interface

        subroutine eigen_modes_copy(modes_in,modes_out) bind( c, name="modes_copy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)         :: modes_in
           type(C_PTR),intent(inout)      :: modes_out
        end subroutine

        subroutine eigen_alloc_arr(N,modes) bind( c, name="modes_alloc" )
           use, intrinsic :: iso_c_binding
           integer( kind = c_int ),value           :: N
           type(C_PTR),intent(inout)               :: modes
        end subroutine

        subroutine eigen_init_mode(modes,i,nnz,dimH,row,col,val) bind( c, name="modes_init" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)               :: modes
           integer( kind = c_int ),value           :: i
           integer( kind = c_int ),value           :: nnz
           integer( kind = c_int ),intent(in)      :: dimH(*)
           integer( kind = c_int ),intent(in)      :: row(*),col(*)
           real( kind = c_double ),intent(in)      :: val(*)
        end subroutine

        subroutine eigen_get_mode_i(modes,i_mode,vec_in,vec_out) bind( c, name="get_mode_i")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: modes
           integer( kind = c_int ),value           :: i_mode
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_mode_reduce(modes,i_mode,vec_in,vec_out) bind( c, name="mode_reduce")
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: modes
           integer( kind = c_int ),value           :: i_mode
           real( kind = c_double ),intent(in)      :: vec_in(*)
           real( kind = c_double ),intent(inout)   :: vec_out(*)
        end subroutine

        subroutine eigen_mode_destroy(modes) bind( c, name="mode_destroy" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(inout)         :: modes
        end subroutine


!        subroutine eigen_mode_init(nnz,dimH,row,col,val,H) bind( c, name="mode_init" )
!           use, intrinsic :: iso_c_binding
!           integer( kind = c_int ),value           :: nnz
!           integer( kind = c_int ),intent(in)      :: dimH(*)
!           integer( kind = c_int ),intent(in)      :: row(*),col(*)
!           real( kind = c_double ),intent(in)      :: val(*)
!           type(C_PTR),intent(inout)               :: modes
!        end subroutine

!        subroutine eigen_H_bcast(id,mas,ismas,mat,comm) bind( c, name="eigen_H_bcast" )
!           use, intrinsic :: iso_c_binding
!           integer( kind = c_int ),value            :: id    !MPI-rank of the processor
!           integer( kind = c_int ),value            :: mas   !MPI-rank of the master
!           logical( kind = c_bool ),value           :: ismas !is master to bcast from
!           type(C_PTR),intent(inout)                :: mat   !matrix to broadcast
!           integer( kind = c_int ),intent(in)       :: comm  !MPI-communicator  the processor
!        end subroutine
!        
!        subroutine eigen_H_mult_mat_vec(mat,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
!
!        subroutine eigen_H_mult_vec_mat_cont(mat,bnd_min,bnd_max,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat_cont" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           integer(c_int),value                    :: bnd_min,bnd_max               
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
!
!        subroutine eigen_H_mult_vec_mat_disc(mat,N,ind,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat_disc" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           integer(c_int),value                    :: N 
!           integer(c_int),intent(in)               :: ind(*)
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
!
!        subroutine eigen_H_mult_mat_vec_cont(mat,bnd_min,bnd_max,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec_cont" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           integer(c_int),value                    :: bnd_min,bnd_max               
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
!
!        subroutine eigen_H_mult_mat_vec_disc(mat,N,ind,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec_disc" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           integer(c_int),value                    :: N 
!           integer(c_int),intent(in)               :: ind(*)
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
!
!        subroutine eigen_H_mult_mat_vec_single(mat,bnd_min,bnd_max,vec_in,vec_out) bind( c, name="eigen_H_mult_mat_vec_single" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           integer(c_int),value                    :: bnd_min,bnd_max               
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
!        
!        subroutine eigen_H_mult_vec_mat(mat,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
!
!        subroutine eigen_H_mult_vec_mat_single(mat,bnd_min,bnd_max,vec_in,vec_out) bind( c, name="eigen_H_mult_vec_mat_single" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(in)                  :: mat
!           integer(c_int),value                    :: bnd_min,bnd_max               
!           real( kind = c_double ),intent(in)      :: vec_in(*)
!           real( kind = c_double ),intent(inout)   :: vec_out(*)
!        end subroutine
!
!        subroutine eigen_H_eval_single(ind,dim_mode,vec_l,vec_r,H,E) bind( c, name="eigen_H_eval_single" )
!           use, intrinsic :: iso_c_binding
!           integer( kind = c_int ),value           :: ind,dim_mode
!           real( kind = c_double ),intent(in)      :: vec_l(*),vec_r(*)
!           type(C_PTR),intent(in)                  :: H
!           real( kind = c_double ),intent(out)     :: E
!        end subroutine
!        
!        
!        subroutine eigen_H_add(H_sum,H_add) bind( c, name="eigen_H_add" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(inout)      :: H_sum
!           type(C_PTR),intent(in)         :: H_add
!        end subroutine
!        
!        subroutine eigen_H_destroy(H) bind( c, name="eigen_H_destroy" )
!           use, intrinsic :: iso_c_binding
!           type(C_PTR),intent(inout)         :: H
!        end subroutine
    
    end interface
    
end module 
