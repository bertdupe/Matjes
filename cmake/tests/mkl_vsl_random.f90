include 'mkl_vsl.f90' 

program randomdist
    use MKL_VSL_TYPE
    use MKL_VSL

    real(8) :: kt,sigma,r(1000)
    integer(4) :: stat
    type(VSL_STREAM_STATE)    :: stream
    integer :: brng,method,seed,n

    kt=5.0d0
    sigma=1.0d0
    n=1000
    r=0.0
    brng=VSL_BRNG_MT19937
    method=VSL_RNG_METHOD_GAUSSIAN_BOXMULLER
    seed=7

    ! initialization
    stat=vslnewstream( stream, brng,  seed )

    stat=vdrnggaussian(method,stream,n,r,kt,sigma)

    ! destroy
    stat=vsldeletestream( stream )

    ! write(*,*) r
end program 
