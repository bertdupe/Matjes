
subroutine tightbinding(lat)
    use m_tb_params, only: TB_params, set_TB_params
    use m_derived_types, only : lattice
    use m_tightbinding_r, only: tightbinding_r
    use m_tightbinding_k, only: tightbinding_k

    implicit none
    ! internal parameter
    type(lattice), intent(in) :: lat

    !read tight-binding io parameter from input and set them in TB_params(m_tb_params)
    Call set_TB_params(lat)
    !do real-space tight-binding stuff
    if(TB_params%flow%do_r) Call tightbinding_r(lat,TB_params%H)   
    !do reciprocal-space tight-binding stuff
    if(TB_params%flow%do_k) Call tightbinding_k(lat)

end subroutine tightbinding
