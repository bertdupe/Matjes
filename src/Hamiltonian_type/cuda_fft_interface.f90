module m_cuda_fft
    use,intrinsic :: iso_c_binding
    implicit none
    public
    interface
        !subroutine eigen_H_init(nnz,dimH,row,col,val,H) bind( c, name="cuda_fft_init" )
        subroutine cuda_fft_init(dim_mode,N_rep,m_n,m_f,h_n,h_f) bind( c, name="cuda_fft_init" )
            use, intrinsic :: iso_c_binding
            integer,intent(in)                      :: dim_mode
            integer,intent(in)                      :: N_rep(3)
            type(C_PTR),intent(inout)               :: m_n
            type(C_PTR),intent(inout)               :: m_f
            type(C_PTR),intent(inout)               :: h_n
            type(C_PTR),intent(inout)               :: h_f
        end subroutine
    
        subroutine cuda_fft_set_real(length,arr_in,arr_device) bind( c, name="cuda_fft_set_real" )
            use, intrinsic :: iso_c_binding
            integer,intent(in)                      :: length
            real( kind = c_double ),intent(in)      :: arr_in(*)
            type(C_PTR),intent(inout)               :: arr_device
        end subroutine

        subroutine cuda_fft_get_real(length,arr_device,arr_out) bind( c, name="cuda_fft_get_real" )
            use, intrinsic :: iso_c_binding
            integer,intent(in),value                :: length
            real( kind = c_double ),intent(inout)   :: arr_out(*)
            type(C_PTR),intent(in)                  :: arr_device
        end subroutine

        subroutine set_cuda_fft_plan(dim_mode,N_rep,plan_fwd,plan_bwd) bind( c, name="set_cuda_fft_plan")
            use, intrinsic :: iso_c_binding
            integer,intent(in)                      :: dim_mode
            integer,intent(in)                      :: N_rep(3)
            type(C_PTR),intent(inout)               :: plan_fwd
            type(C_PTR),intent(inout)               :: plan_bwd
        end subroutine

        subroutine  cuda_fft_set_operator(dim_mode,N_rep,arr_n,op_f) bind( c, name='cuda_fft_set_operator')
            use, intrinsic :: iso_c_binding
            integer,intent(in)                      :: dim_mode
            integer,intent(in)                      :: N_rep(3)
            real( kind = c_double ),intent(in)      :: arr_n(*)
            type(C_PTR),intent(inout)               :: op_f
        end subroutine

        subroutine cuda_fft_calc_H(dim_mode,N_rep,m_n,m_f,k_f,h_n,h_f,plan_fwd,plan_bwd) bind( c, name='cuda_fft_calc_h')
            use, intrinsic :: iso_c_binding
            integer,intent(in)                      :: dim_mode
            integer,intent(in)                      :: N_rep(3)
            type(C_PTR),intent(inout)               :: m_n
            type(C_PTR),intent(inout)               :: m_f
            type(C_PTR),intent(inout)               :: k_f
            type(C_PTR),intent(inout)               :: h_n
            type(C_PTR),intent(inout)               :: h_f
            type(C_PTR),intent(inout)               :: plan_fwd
            type(C_PTR),intent(inout)               :: plan_bwd
        end subroutine

        subroutine cuda_fft_copy_cmplx(N,arr_in,arr_out) bind( c, name='cuda_fft_copy_cmplx')
            use, intrinsic :: iso_c_binding
            integer(C_INT),intent(in),value         :: N
            type(C_PTR),intent(in)                  :: arr_in
            type(C_PTR),intent(inout)               :: arr_out
        end subroutine
    end interface
end module 
