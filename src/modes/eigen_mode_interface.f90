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

        subroutine eigen_get_mode_single_disc(modes,cptr,i_mode,N_in,ind_in,N_out,ind_out,val_out) bind( c, name="get_ind_disc" )
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: modes
           type(C_PTR),intent(in)                  :: cptr(*)
           integer( kind = c_int   ),value         :: i_mode
           integer( kind = c_int   ),value         :: N_in
           integer( kind = c_int   ),intent(in)    :: ind_in(*)
           integer( kind = c_int   ),value         :: N_out
           integer( kind = c_int   ),intent(inout) :: ind_out(*)
           real   ( kind = c_double),intent(inout) :: val_out(*)
        end subroutine

        subroutine eigen_get_disc(modes,cptr,N_ind,ind,vec) bind(c, name='get_disc')
           use, intrinsic :: iso_c_binding
           type(C_PTR),intent(in)                  :: modes
           type(C_PTR),intent(in)                  :: cptr(*)
           integer( kind = c_int   ),value         :: N_ind
           integer( kind = c_int   ),intent(in)    :: ind(*)
           real   ( kind = c_double),intent(inout) :: vec(*)
        end subroutine
    end interface
    
end module 
