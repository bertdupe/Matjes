module m_energy_solve_dense
use m_tb_types
implicit none
public
contains




    subroutine Hr_eigval_zheev(h_par,Hr,eigval)
        type(parameters_TB_Hsolve),intent(in)     :: h_par
        complex(8),intent(in)                     ::  Hr(h_par%dimH,h_par%dimH)
        real(8),intent(out),allocatable           ::  eigval(:)

        integer                     :: dimH
        complex(8)                  :: H_loc(h_par%dimH,h_par%dimH)
        real(kind=8)                :: RWORK(3*h_par%dimH-2)
        integer                     :: info,l_work
        complex(kind=8),allocatable :: WORK(:)
       
        dimH=h_par%dimH
        allocate(eigval(h_Par%dimH),source=0.0d0)
        H_loc=Hr
        allocate(WORK(1))
        call ZHEEV( 'V', 'U', dimH, H_loc, dimH, eigval, WORK, -1, RWORK, INFO )
        l_work=int(work(1))
        deallocate(work)
        allocate(work(l_work),source=cmplx(0.0d0,0.0d0,8))
        call ZHEEV( 'N', 'U', dimH, H_loc, dimH, eigval, WORK, l_work, RWORK, INFO )
    end subroutine


    subroutine Hr_eigval_zheevd(h_par,Hr,eigval)
        type(parameters_TB_Hsolve),intent(in)     :: h_par
        complex(8),intent(in)       ::  Hr(h_par%dimH,h_par%dimH)
        real(8),intent(out),allocatable :: eigval(:)

        integer                     :: dimH
        complex(8)                  :: H_loc(h_par%dimH,h_par%dimH)
        integer                     :: info,lwork
        integer                     :: lrwork,liwork
        integer,allocatable         :: iwork(:)
        complex(kind=8),allocatable :: WORK(:)
        real(8),allocatable         :: RWORK(:)

        integer                     :: tmp_iwork(1)
        real(8)                     :: tmp_rwork(1)
        complex(8)                  :: tmp_work(1)

        dimH=h_par%dimH
        H_loc=Hr
        allocate(eigval(h_Par%dimH),source=0.0d0)
        call ZHEEVD( 'V', 'U', dimH, H_loc, dimH, eigval, tmp_WORK, -1, tmp_RWORK, -1,tmp_IWORK,-1,INFO )
        lwork=int(tmp_work(1))
        LIWORK=tmp_IWORK(1)
        LRWORK=int(tmp_rwork(1))
        allocate(work(lwork),source=cmplx(0.0d0,0.0d0,8))
        allocate(iwork(liwork),source=0)
        allocate(rwork(lrwork),source=0.0d0)
        call ZHEEVD( 'V', 'U', dimH, H_loc, dimH, eigval, WORK, lwork, RWORK, LRWORK,IWORK,LIWORK,INFO )
    end subroutine 

    subroutine Hr_eigval_feast(h_par,Hr,eigval)
        type(parameters_TB_Hsolve),intent(in)   :: h_par
        complex(8),intent(in)                   :: Hr(h_par%dimH,h_par%dimH)
        real(8),intent(out),allocatable         :: eigval(:)

        complex(8),allocatable      :: eigvec(:,:)

        Call Hr_eigvec_feast(h_par,Hr,eigvec,eigval)
    end subroutine 

    subroutine Hr_eigval_zheevr(h_par,Hr,eigval)
        type(parameters_TB_Hsolve),intent(in)     :: h_par
        complex(8),intent(in)                     ::  Hr(h_par%dimH,h_par%dimH)
        real(8),intent(out),allocatable           ::  eigval(:)

        integer                     :: dimH,lda
        complex(8)                  :: H_loc(h_par%dimH,h_par%dimH)
        real(8)                     :: abstol
        integer                     :: ldz,il,iu
        integer                     :: m
        integer,allocatable         :: isuppz(:)
        real(8)                     :: w(h_par%dimH)
        complex(8),allocatable      :: z(:,:)
        real(8)                     :: vl,vu
        !work data
        integer                     :: lwork,lrwork,liwork
        complex(8),allocatable      :: work(:)
        real(8),allocatable         :: rwork(:)
        integer,allocatable         :: iwork(:)
        
        integer                     :: info
        real(8),external            :: DLAMCH

        dimH=h_par%dimH
        lda=dimH
        vl=h_par%extE(1); vu=h_par%extE(2)
        H_loc=Hr
        abstol=DLAMCH('S')
        il=0;iu=0
        allocate(isuppz(2*dimH))
        ldz=1
        allocate(z(ldz,1))
#if 0
        lwork=2*dimH; lrwork=24*dimH; liwork=10*dimH
#else
        lwork=-1; lrwork=-1; liwork=-1
        allocate(work(1),rwork(1),iwork(1))
        Call ZHEEVR('N' , 'V' , 'U' , dimH , H_loc , lda , vl , vu , il , iu , abstol , m , w , z , ldz , isuppz , work , lwork , rwork , lrwork , iwork , liwork , info )
        if(info/=0) STOP "ZHEERV failed at setup step"
        lwork=int(work(1)); lrwork=int(rwork(1)); liwork=iwork(1)
        deallocate(work,rwork,iwork)
#endif
        allocate(work(lwork),rwork(lrwork),iwork(liwork))
        Call ZHEEVR('N' , 'V' , 'U' , dimH , H_loc , lda , vl , vu , il , iu , abstol , m , w , z , ldz , isuppz , work , lwork , rwork , lrwork , iwork , liwork , info )
        if(info/=0) STOP "ZHEERV failed at setup step"
        allocate(eigval,source=w(1:m))
    end subroutine


    subroutine Hr_eigvec_zheev(h_par,Hr,eigvec,eigval)
        type(parameters_TB_Hsolve),intent(in)   ::  h_par
        complex(8),intent(in)                   ::  Hr(h_par%dimH,h_par%dimH)
        complex(8),intent(out),allocatable      ::  eigvec(:,:)
        real(8),intent(out),allocatable         ::  eigval(:)

        real(kind=8)                :: RWORK(3*h_par%dimH-2)
        integer                     :: info,l_work
        complex(kind=8),allocatable :: WORK(:)

        allocate(eigvec,source=Hr)
        allocate(eigval(h_Par%dimH),source=0.0d0)
        allocate(work(1),source=cmplx(0.0d0,0.0d0,8))
        call ZHEEV( 'V', 'U', h_Par%dimH, eigvec, h_Par%dimH, eigval, WORK, -1, RWORK, INFO )
        l_work=int(work(1))
        deallocate(work)
        allocate(work(l_work),source=cmplx(0.0d0,0.0d0,8))
        call ZHEEV( 'V', 'U', h_Par%dimH, eigvec, h_Par%dimH, eigval, WORK, l_work, RWORK, INFO )
    end subroutine 

    subroutine Hr_eigvec_zheevd(h_par,Hr,eigvec,eigval)
        type(parameters_TB_Hsolve),intent(in)     :: h_par
        complex(8),intent(in)       ::  Hr(h_par%dimH,h_par%dimH)
        complex(8),intent(out),allocatable      ::  eigvec(:,:)
        real(8),intent(out),allocatable         ::  eigval(:)

        integer                     :: info,lwork
        integer                     :: lrwork,liwork
        integer,allocatable         :: iwork(:)
        complex(kind=8),allocatable :: WORK(:)
        real(8),allocatable         :: RWORK(:)

        integer                     :: tmp_iwork(1)
        real(8)                     :: tmp_rwork(1)
        complex(8)                  :: tmp_work(1)

        allocate(eigvec,source=Hr)
        allocate(eigval(h_Par%dimH),source=0.0d0)
        call ZHEEVD( 'V', 'U', h_par%dimH, eigvec, h_par%dimH, eigval, tmp_WORK, -1, tmp_RWORK, -1,tmp_IWORK,-1,INFO )
        lwork=int(tmp_work(1))
        LIWORK=tmp_IWORK(1)
        LRWORK=int(tmp_rwork(1))
        allocate(work(lwork),source=cmplx(0.0d0,0.0d0,8))
        allocate(iwork(liwork),source=0)
        allocate(rwork(lrwork),source=0.0d0)
        call ZHEEVD( 'V', 'U', h_par%dimH, eigvec, h_par%dimH, eigval, WORK, lwork, RWORK, LRWORK,IWORK,LIWORK,INFO )
    end subroutine 

    subroutine Hr_eigvec_feast(h_par,Hr,eigvec,eigval)
        type(parameters_TB_Hsolve),intent(in)   :: h_par
        complex(8),intent(in)                   :: Hr(h_par%dimH,h_par%dimH)
        complex(8),intent(out),allocatable      :: eigvec(:,:)
        real(8),intent(out),allocatable         :: eigval(:)

        integer                     :: dimH
        integer                     :: fpm(128)
        complex(8)                  :: Hin(h_par%dimH,h_par%dimH)
        real(8)                     :: emin,emax
        real(8)                     :: epsout
        integer                     :: loop
        integer                     :: m0,m
        real(8),allocatable         :: res(:)
        integer                     :: info
        complex(8),allocatable      :: x(:,:)
        real(8),allocatable         :: e(:)

        dimH=h_par%dimH
        Call feastinit(fpm) 
        fpm(1)=1
        fpm(2)=4
        emin=h_par%extE(1)
        emax=h_par%extE(2)
        m0=h_par%estNe
        if(m0==0.or.m0>dimH) m0=h_par%dimH
        allocate(res(m0),source=0.0d0)
        allocate(e(m0),source=0.0d0)
        allocate(x(h_par%dimh,m0),source=cmplx(0.0d0,0.0d0,8))
        Hin=Hr
        call zfeast_heev ( 'F' , dimH , Hin , dimH , fpm , epsout , loop , emin , emax , m0 , e , x , m , res , info )
        write(*,*) 'done'
        allocate(eigval,source=e(1:m))
        allocate(eigvec,source=x(1:dimH,1:m))
        deallocate(x,e)
        if(info/=0) STOP 'info of zfest_heev not zero'
    end subroutine 

    subroutine Hr_eigvec_zheevr(h_par,Hr,eigvec,eigval)
        type(parameters_TB_Hsolve),intent(in)     ::  h_par
        complex(8),intent(in)                     ::  Hr(h_par%dimH,h_par%dimH)
        real(8),intent(out),allocatable           ::  eigval(:)
        complex(8),intent(out),allocatable        ::  eigvec(:,:)

        integer                     :: dimH,lda
        complex(8)                  :: H_loc(h_par%dimH,h_par%dimH)
        real(8)                     :: abstol
        integer                     :: ldz,il,iu
        integer                     :: m
        integer,allocatable         :: isuppz(:)
        real(8)                     :: w(h_par%dimH)
        complex(8),allocatable      :: z(:,:)
        real(8)                     :: vl,vu
        !work data
        integer                     :: lwork,lrwork,liwork
        complex(8),allocatable      :: work(:)
        real(8),allocatable         :: rwork(:)
        integer,allocatable         :: iwork(:)
        
        integer                     :: info
        real(8),external            :: DLAMCH

        dimH=h_par%dimH
        lda=dimH
        vl=h_par%extE(1); vu=h_par%extE(2)
        H_loc=Hr
        abstol=DLAMCH('S')
        il=0;iu=0
        allocate(isuppz(2*dimH))
        ldz=dimH
        allocate(z(ldz,dimH)) !could also use estNe, if memory is scarse 
#if 0
        lwork=2*dimH; lrwork=24*dimH; liwork=10*dimH
#else
        lwork=-1; lrwork=-1; liwork=-1
        allocate(work(1),rwork(1),iwork(1))
        Call ZHEEVR('V' , 'V' , 'U' , dimH , H_loc , lda , vl , vu , il , iu , abstol , m , w , z , ldz , isuppz , work , lwork , rwork , lrwork , iwork , liwork , info )
        if(info/=0) STOP "ZHEERV failed at setup step"
        lwork=int(work(1)); lrwork=int(rwork(1)); liwork=iwork(1)
        deallocate(work,rwork,iwork)
#endif
        allocate(work(lwork),rwork(lrwork),iwork(liwork))
        Call ZHEEVR('V' , 'V' , 'U' , dimH , H_loc , lda , vl , vu , il , iu , abstol , m , w , z , ldz , isuppz , work , lwork , rwork , lrwork , iwork , liwork , info )
        if(info/=0) STOP "ZHEERV failed at setup step"
        allocate(eigval,source=w(1:m))
        allocate(eigvec,source=z(1:dimH,1:m))
    end subroutine


end module
