module m_laguerre_interface
    use,intrinsic :: iso_c_binding
    implicit none
    contains
    subroutine laguerre_polynomial(l, p, x, res)
        use, intrinsic :: iso_c_binding
        external laguerre
        integer(kind=c_int), intent(in)    :: l, p
        real(kind=c_double), intent(in)    :: x(:)
        real(kind=c_double), intent(inout) :: res(:)
        Call laguerre(l, p, x, res, size(x))
    end subroutine

    subroutine laguerre_polynomial_scalar(l, p, x, res)
        use, intrinsic :: iso_c_binding
        external laguerre
        integer(kind=c_int), intent(in)    :: l, p
        real(kind=c_double), intent(in)    :: x
        real(kind=c_double), intent(inout) :: res
        Call laguerre_scalar(l, p, x, res)
    end subroutine
end module
