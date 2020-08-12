module m_energy_solve_dense
implicit none
public
contains

    subroutine Hr_eigval_zheev(dimH,Hr,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(in)       ::  Hr(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(8)                  :: H_loc(dimH,dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info,l_work
        complex(kind=8),allocatable :: WORK(:)
        complex(8)                  :: init_WORK(1)
        
        H_loc=Hr
        call ZHEEV( 'V', 'U', dimH, H_loc, dimH, eigval, init_WORK, -1, RWORK, INFO )
        l_work=int(init_work(1))
        allocate(work(l_work),source=cmplx(0.0d0,0.0d0,8))
        call ZHEEV( 'N', 'U', dimH, H_loc, dimH, eigval, WORK, l_work, RWORK, INFO )
    end subroutine


    subroutine Hr_eigval_zheevd(dimH,Hr,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(in)       ::  Hr(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(8)                  :: H_loc(dimH,dimH)
        integer                     :: info,lwork
        integer                     :: lrwork,liwork
        integer,allocatable         :: iwork(:)
        complex(kind=8),allocatable :: WORK(:)
        real(8),allocatable         :: RWORK(:)

        integer                     :: tmp_iwork(1)
        real(8)                     :: tmp_rwork(1)
        complex(8)                  :: tmp_work(1)

        H_loc=Hr
        eigval=0.0d0
        call ZHEEVD( 'V', 'U', dimH, H_loc, dimH, eigval, tmp_WORK, -1, tmp_RWORK, -1,tmp_IWORK,-1,INFO )
        lwork=int(tmp_work(1))
        LIWORK=tmp_IWORK(1)
        LRWORK=int(tmp_rwork(1))
        allocate(work(lwork),source=cmplx(0.0d0,0.0d0,8))
        allocate(iwork(liwork),source=0)
        allocate(rwork(lrwork),source=0.0d0)
        call ZHEEVD( 'V', 'U', dimH, H_loc, dimH, eigval, WORK, lwork, RWORK, LRWORK,IWORK,LIWORK,INFO )
    end subroutine 


    subroutine Hr_eigvec_zheev(dimH,Hr,eigvec,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(out)      ::  eigvec(dimH,dimH)
        complex(8),intent(in)       ::  Hr(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        complex(kind=8)             :: init_WORK(2*dimH)
        real(kind=8)                :: RWORK(3*dimH-2)
        integer                     :: info,l_work
        complex(kind=8),allocatable :: WORK(:)

        eigvec=Hr
        eigval=0.0d0
        call ZHEEV( 'V', 'U', dimH, eigvec, dimH, eigval, init_WORK, -1, RWORK, INFO )
        l_work=int(init_work(1))
        allocate(work(l_work),source=cmplx(0.0d0,0.0d0,8))
        call ZHEEV( 'V', 'U', dimH, eigvec, dimH, eigval, WORK, l_work, RWORK, INFO )
    end subroutine 

    subroutine Hr_eigvec_zheevd(dimH,Hr,eigvec,eigval)
        integer,intent(in)          ::  dimH
        complex(8),intent(in)       ::  Hr(dimH,dimH)
        complex(8),intent(out)      ::  eigvec(dimH,dimH)
        real(8),intent(out)         ::  eigval(dimH)

        integer                     :: info,lwork
        integer                     :: lrwork,liwork
        integer,allocatable         :: iwork(:)
        complex(kind=8),allocatable :: WORK(:)
        real(8),allocatable         :: RWORK(:)

        integer                     :: tmp_iwork(1)
        real(8)                     :: tmp_rwork(1)
        complex(8)                  :: tmp_work(1)

        eigvec=Hr
        eigval=0.0d0
        call ZHEEVD( 'V', 'U', dimH, eigvec, dimH, eigval, tmp_WORK, -1, tmp_RWORK, -1,tmp_IWORK,-1,INFO )
        lwork=int(tmp_work(1))
        LIWORK=tmp_IWORK(1)
        LRWORK=int(tmp_rwork(1))
        allocate(work(lwork),source=cmplx(0.0d0,0.0d0,8))
        allocate(iwork(liwork),source=0)
        allocate(rwork(lrwork),source=0.0d0)
        call ZHEEVD( 'V', 'U', dimH, eigvec, dimH, eigval, WORK, lwork, RWORK, LRWORK,IWORK,LIWORK,INFO )
    end subroutine 

    subroutine Hr_eigvec_feast(dimH,Hr,eigvec,eigval)
        integer,intent(in)          :: dimH
        complex(8),intent(in)       :: Hr(dimH,dimH)
        complex(8),intent(out)      :: eigvec(dimH,dimH)
        real(8),intent(out)         :: eigval(dimH)

        integer                     :: fpm(128)
        real(8)                     :: emin,emax
        complex(8)                  :: Hin(dimH,dimH)
        real(8)                     :: epsout
        integer                     :: loop
        integer                     :: m0,m
        real(8),allocatable         :: res(:)
        integer                     :: info

        Call feastinit(fpm) 
        fpm(1)=1
        emin=-1.0d1
        emax=1.0d1
        m0=dimH
        allocate(res(dimH),source=0.0d0)
        Hin=Hr
        eigval=0.0d0
        call zfeast_heev ( 'F' , dimH , Hin , dimH , fpm , epsout , loop , emin , emax , m0 , eigval , eigvec , m , res , info )
        write(*,*) maxval(res)
        write(*,*) eigval
        write(*,*) dimH,m
        if(m<dimH) STOP "not all eigenvalues found"
        if(info/=0) STOP 'info of zfest_heev not zero'
    end subroutine 

end module
