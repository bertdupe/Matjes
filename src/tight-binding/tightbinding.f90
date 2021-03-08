
subroutine tightbinding(lat)
    use m_TB_types, only: parameters_TB
    use m_derived_types, only : lattice
    use m_tightbinding_r, only: tightbinding_r
    use m_tightbinding_k, only: tightbinding_k
    use m_rw_TB, only:  rw_TB

    implicit none
    ! internal parameter
    type(lattice), intent(in)   :: lat
    type(parameters_TB)         :: tb_par

    !read tight-binding io parameter from input and set TB_params(m_tb_params)
    call rw_TB(tb_par,'input')
    Call tb_par%init(lat)
    !do real-space tight-binding stuff
    if(tb_par%flow%do_r) Call tightbinding_r(lat,tb_par)   
    !do reciprocal-space tight-binding stuff
    if(tb_par%flow%do_k) Call tightbinding_k(lat,tb_par)

end subroutine tightbinding
