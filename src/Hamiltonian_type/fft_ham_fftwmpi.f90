module m_fft_H_fftwmpi
!module which contains the general discrete fourier transform Hamiltonian based on FFTW3-MPI

!the type is used by first calling init_shape, followed by init_op with the operator as described in more detail for the dipolar_fft interaction

use,intrinsic :: ISO_FORTRAN_ENV, only: error_unit
use m_fftw3
use m_fftwmpi    ! not supposed to be used in production. for tests only
use mpi_util
use m_fft_H_base, only: fft_H
use m_type_lattice,only: lattice
use m_H_type, only: len_desc
use mpi_basic

private
public fft_H_fftwmpi

type,extends(fft_H) ::  fft_H_fftwmpi
    private
    type(c_ptr)     :: plan_mag_F=c_null_ptr   !FFTW plan M_n -> M_F
    type(c_ptr)     :: plan_H_I=c_null_ptr     !FFTW plan H_F -> H_n

    !fourier transformated arrays
    complex(C_DOUBLE_COMPLEX),allocatable   ::  M_F(:,:)    !magnetization in fourier-space
    complex(C_DOUBLE_COMPLEX),allocatable   ::  K_F(:,:,:)  !demagnetization tensor in fourier-space
    complex(C_DOUBLE_COMPLEX),allocatable   ::  H_F(:,:)    !effective field fourier-space
    type(mpi_type)                          ::  com_fftw
    !data stored in base class
!    real(C_DOUBLE),allocatable              ::  M_n(:,:)    !magnetization in normal-space
!    real(C_DOUBLE),allocatable              ::  H_n(:,:)    !effective field normal-space
contains
    procedure,public    :: get_H            !get effective field
    procedure,public    :: get_H_single     !get effective field for a single site
    procedure,public    :: init_shape       !initializes the shape and the H and M operators, arrays
    procedure,public    :: init_op          !initializes the K_F array

    !utility functions
    procedure,public    :: destroy          !destroys all data of this type
    procedure,public    :: copy             !copy this to new instance
    procedure,public    :: mv               !mv this to new instance
    procedure,public    :: add              !adds two fft_H in same space by adding their K_F
    procedure,public    :: bcast=>bcast_fft !mv this to new instance


    !internal procedures
    procedure           :: set_M                !set internal magnetization in normal-space from lattice
end type

contains

subroutine add(this,H_in)
    class(fft_H_fftwmpi),intent(inout)  :: this
    class(fft_H),intent(in)             :: H_in

    logical :: is_set

    if(.not.H_in%is_set()) ERROR STOP "CANNOT ADD fft_H as H_in is not set"
    is_set=this%is_set()
    if(is_set)then
        if(.not.this%same_space(H_in)) ERROR STOP "CANNOT ADD fft_H as this and H_in do not act on same space"
        Call this%fft_H%add(H_in)
        select type(H_in)
        type is(fft_H_fftwmpi)
            this%K_F=this%K_F+H_in%K_F
        class default
            ERROR STOP "HAS TO HAVE SAME TYPE"
        end select
    else
        Call H_in%copy(this)
    endif
end subroutine

subroutine bcast_fft(this,comm)
    class(fft_H_fftwmpi),intent(inout)  ::  this
    type(mpi_type),intent(in)           ::  comm
    integer ::  shp(2)

    Call this%fft_H%bcast(comm)
    if(comm%ismas) shp=shape(this%M_F)
    Call bcast(shp,comm)
    if(.not.comm%ismas)then
        allocate(this%M_F(shp(1),shp(2)))
        allocate(this%H_F(shp(1),shp(2)))
    endif
    Call bcast_alloc(this%K_F,comm)
!    if(.not.comm%ismas) Call set_plans(this)
    call this%com_fftw%copy_base(comm)
    Call set_plans(this)
end subroutine

subroutine mv(this,H_out)
    class(fft_H_fftwmpi),intent(inout)  :: this
    class(fft_H),intent(inout)          :: H_out

    if(.not.this%set)then
        ERROR STOP "CANNOT MV UNINITIALIZED FFT_H"
    endif
    Call this%fft_H%mv(H_out)
    select type(H_out)
    type is(fft_H_fftwmpi)
        call move_alloc(this%M_F,H_out%M_F)
        call move_alloc(this%H_F,H_out%H_F)
        call move_alloc(this%K_F,H_out%K_F)
    class default
        ERROR STOP "HAS TO HAVE SAME TYPE"
    end select
    Call this%destroy()
end subroutine

subroutine copy(this,H_out)
    class(fft_H_fftwmpi),intent(in):: this
    class(fft_H),intent(inout)      :: H_out

    Call this%fft_H%copy(H_out)
    select type(H_out)
    type is(fft_H_fftwmpi)
        allocate(H_out%M_F,mold=this%M_F)
        allocate(H_out%H_F,mold=this%H_F)
        allocate(H_out%K_F,source=this%K_F)
        Call set_plans(H_out)
    class default
        ERROR STOP "HAS TO HAVE SAME TYPE"
    end select
end subroutine

subroutine destroy(this)
    class(fft_H_fftwmpi),intent(inout)     :: this

    Call this%fft_H%destroy()
#ifdef CPP_FFTWMPI
    if(c_associated(this%plan_mag_f)) Call fftw_destroy_plan(this%plan_mag_f)
    if(c_associated(this%plan_H_I))   Call fftw_destroy_plan(this%plan_H_I)
#else
    ERROR STOP "Requires FFTW3"
#endif
    if(allocated(this%M_F)) deallocate(this%M_F)
    if(allocated(this%K_F)) deallocate(this%K_F)
    if(allocated(this%H_F)) deallocate(this%H_F)
end subroutine


subroutine init_op(this,dim_mode,K_n,desc_in)
    !subroutine which initializes the fourier-transformed operator of K, while deallocating K_N
    class(fft_H_fftwmpi),intent(inout)         :: this
    integer,intent(in)                         :: dim_mode
    real(8),intent(inout),allocatable          :: K_n(:,:,:)
    character(len=*),intent(in),optional       :: desc_in

#ifdef CPP_FFTWMPI
    integer              :: Nk_tot           !number of state considered in FT (product of N_rep)
    integer(C_INT)       :: N_rep_rev(3)     !reversed N_rep necessary for fftw3 (col-major -> row-major)
    integer(C_INT)       :: howmany          !dimension of quantitiy which is fourier-transformed (see FFTW3)
    type(c_ptr)          :: plan_K_F         !plan for fourier transformation of K
    integer :: i,j,k

    Call this%fft_H%init_op(dim_mode,K_n,desc_in)

    !check shapes etc.
    if(.not.allocated(this%M_n))then
        write(error_unit,'(///A)') "Magnetic order parameter not set when tying to initialize operator in fft_H-type"
        write(error_unit,'(A)') "This is most probably the case because init_shape has not been called previously"
        ERROR STOP
    endif
    if(allocated(this%K_F))then
        write(error_unit,'(///A)') "Trying to initialize operator of fft_H, but the fourier-transform has already been set"
        write(error_unit,'(A)') "init_op shall only be called once"
        ERROR STOP
    endif
    if(size(K_N,3)/=product(this%N_rep))then
        write(error_unit,'(///A)') "Failed to set the fft-operator since the number of repetition shape does not fit"
        write(error_unit,'(A,I12,A,3I8)') "fft_H Nk_total",product(this%N_rep)," coming from this%N_rep:",this%N_rep
        write(error_unit,'(A,I12)')       "input Nk_total",size(K_n,2)
        ERROR STOP
    endif
    if(size(K_N,1)/=size(this%M_n,1).or.size(K_N,2)/=size(this%M_n,1))then
        write(error_unit,'(///A)') "Failed to set the fft-operator since the mode dimension does not fit"
        write(error_unit,'(A,I6)')       "fft_H dim_mode:     ",size(this%M_n,1)
        write(error_unit,'(A,2I6)')      "K_n input dim_modes:",size(K_n,1),size(K_n,2)
        ERROR STOP
    endif


    !calculate fourier transform of K and save it in dipole-type
    howmany=int(dim_mode**2,C_int)
    Nk_tot=product(this%N_rep)
    N_rep_rev=this%N_rep(size(this%N_rep):1:-1)
    allocate(this%K_F(dim_mode,dim_mode,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
    plan_K_F= fftw_plan_many_dft_r2c(int(3,C_INT), N_rep_rev,  howmany,&
                                    &K_n,          N_rep_rev,&
                                    &howmany,      int(1,C_int),&
                                    &this%K_F,     N_rep_rev,&
                                    &howmany,      int(1,C_int),&
                                    &FFTW_FORWARD+FFTW_ESTIMATE)
    Call fftw_execute_dft_r2c(plan_K_F, K_n, this%K_F)
    this%K_F=this%K_F/real(product(N_rep_rev),8)
    Call fftw_destroy_plan(plan_K_F)
    this%set=.true.
    !deallocate K_n, since at some point one might want to keep it in here
    deallocate(K_n)

#else
        ERROR STOP "CANNOT USE FFTW-MPI Hamiltonian without FFTW-MPI (CPP_FFTWMPI)"
#endif
end subroutine


subroutine init_shape(this,dim_mode,periodic,dim_lat,Kbd,N_rep)
    !initializes the arrays with whose the work will be done, sets the shapes, and returns shape data necessary for construction of the operator tensor
!$  use omp_lib
    class(fft_H_fftwmpi),intent(inout)  :: this
    integer,intent(in)                  :: dim_mode
    logical,intent(in)                  :: periodic(3)  !T: periodic boundary, F: open boundary
    integer,intent(in)                  :: dim_lat(3)
    integer,intent(out)                 :: Kbd(2,3)
    integer,intent(out)                 :: N_rep(3)

#ifdef CPP_FFTW3
    integer         :: i
    integer         :: Nk_tot           !number of state considered in FT (product of N_rep)


    if(allocated(this%M_n))then
        write(error_unit,'(///A)') "Trying to initialize shape of fft_H, but the fourier-transform has already been set"
        write(error_unit,'(A)') "init_shape shall only be called once"
        ERROR STOP
    endif

    Call this%fft_H%init_shape(dim_mode,periodic,dim_lat,Kbd,N_rep)
#ifdef CPP_FFTWMPI_THREAD
!$  Call fftw_mpi_plan_with_nthreads(omp_get_max_threads())
#endif
    !set order work arrays and fourier transform
    Nk_tot=product(N_rep)
    allocate(this%M_F(dim_mode,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
    allocate(this%H_F(dim_mode,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
    Call set_plans(this)

#else
        ERROR STOP "CANNOT USE FFTW DIPOL (mag_dip_fft) without FFTW (CPP_FFTW3)"
#endif
end subroutine

subroutine set_M(this,lat)
    !set this%M_n from lat%M%modes according to this%M_internal
    class(fft_H_fftwmpi),intent(inout)    ::  this
    type(lattice),intent(in)              ::  lat

    Call this%M_internal(this%M_n,lat%M%modes,lat%dim_lat,this%N_rep,size(this%M_n,1))
end subroutine

subroutine get_H(this,lat,Hout)
    !get effective field H
    class(fft_H_fftwmpi),intent(inout)    ::  this
    type(lattice),intent(in)              ::  lat
    real(8),intent(inout)                 ::  Hout(:,:)
#ifdef CPP_FFTWMPI
    integer ::  i,j,l

    Call this%set_M(lat)
    Call fftw_mpi_execute_dft_r2c(this%plan_mag_F, this%M_n, this%M_F)

    this%H_F=cmplx(0.0d0,0.0d0,8)
    do j=1,size(this%M_F,2)
        do i=1,lat%Nmag*3
            this%H_F(i,j)=sum(this%K_F(:,i,j)*this%M_F(:,j))
        enddo
    enddo

    Call fftw_mpi_execute_dft_c2r(this%plan_H_I, this%H_F, this%H_n)
    Call this%H_internal(this%H_n,Hout,lat%dim_lat,this%N_rep,size(Hout,1))
#else
    ERROR STOP "fft_H%get_H requires CPP_FFTW3"
#endif
end subroutine

subroutine get_H_single(this,lat,site,Hout)
    !get effective field H
    class(fft_H_fftwmpi),intent(inout)    ::  this
    type(lattice),intent(in)              ::  lat
    integer,intent(in)                    ::  site
    real(8),intent(inout)                 ::  Hout(3)
#ifdef CPP_FFTWMPI
    integer ::  i,j,l

    Call this%set_M(lat)
    stop 'toto'
    Call fftw_mpi_execute_dft_r2c(this%plan_mag_F, this%M_n, this%M_F)

    this%H_F=cmplx(0.0d0,0.0d0,8)
    do j=1,size(this%M_F,2)
        do i=1,lat%Nmag*3
            this%H_F(i,j)=sum(this%K_F(:,i,j)*this%M_F(:,j))
        enddo
    enddo
    Call fftw_mpi_execute_dft_c2r(this%plan_H_I, this%H_F, this%H_n)
    Call H_internal_single(this%H_n,Hout,site,lat%dim_lat,this%N_rep,lat%nmag)
#else
    ERROR STOP "fft_H_fftw%get_H_single requires CPP_FFTWMPI"
#endif
end subroutine

subroutine set_plans(this)
    class(fft_H_fftwmpi),intent(inout)  :: this


    integer(C_INTPTR_T) :: N_rep_rev(3)     !reversed N_rep necessary for fftw3 (col-major -> row-major)
    integer(C_INTPTR_T) :: howmany          !dimension of quantitiy which is fourier-transformed (see FFTW3)
    integer(C_INTPTR_T) :: M_offset,alloc_local,local_M
    type(C_PTR)         :: cdata, rdata
    complex(C_DOUBLE_COMPLEX), pointer :: data(:,:,:)
    real(C_DOUBLE), pointer :: in(:,:,:)
    integer :: dim_fft,errcode,ierror

#ifdef CPP_FFTWMPI

    if(.not.allocated(this%M_n)) ERROR STOP "cannot set fftw_plans as M_n not allocated"
    if(.not.allocated(this%M_F)) ERROR STOP "cannot set fftw_plans as M_F not allocated"
    if(.not.allocated(this%H_n)) ERROR STOP "cannot set fftw_plans as H_n not allocated"
    if(.not.allocated(this%H_F)) ERROR STOP "cannot set fftw_plans as H_F not allocated"

    N_rep_rev=this%N_rep(size(this%N_rep):1:-1)
    howmany=int(size(this%M_n,1),C_int)

    ! first check the dimension
    if (N_rep_rev(1)*N_rep_rev(2).eq.1) then     ! the supercell is dimension, parallelization can not be performed
       write(6,'(a)') 'cannot perform FFT for a one dimensional supercell'
       call mpi_abort(this%com_fftw%com,errcode,ierror)
    elseif(N_rep_rev(1).eq.1)  then                   ! redimension for a 2D MPI FFT
       dim_fft=2
    else                         ! redimension for a 3D MPI FFT
       dim_fft=3
    endif

    !   get local data size and allocate (note dimension reversal)
    alloc_local = fftw_mpi_local_size_many(dim_fft, N_rep_rev(:dim_fft), Howmany, &
                                & FFTW_MPI_DEFAULT_BLOCK,this%com_fftw%com,local_M,M_offset)
! allocate the real space in rdata
    rdata = fftw_alloc_real(2*howmany*alloc_local)
    call c_f_pointer(rdata, in, N_rep_rev)

! allocate the complex space in cdata
    cdata = fftw_alloc_complex(howmany*alloc_local)
    call c_f_pointer(cdata, data, N_rep_rev)

    this%plan_mag_F= fftw_mpi_plan_many_dft_r2c(dim_fft, N_rep_rev(:dim_fft), howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                 & in, data,this%com_fftw%com,FFTW_FORWARD+FFTW_MEASURE+FFTW_PATIENT)

    this%plan_H_I= fftw_mpi_plan_many_dft_c2r(dim_fft, N_rep_rev(:dim_fft), howmany, FFTW_MPI_DEFAULT_BLOCK, FFTW_MPI_DEFAULT_BLOCK, &
                                 & data, in,this%com_fftw%com,FFTW_BACKWARD+FFTW_MEASURE+FFTW_PATIENT)

#else
    ERROR STOP "CANNOT USE FFT_H without FFTW (CPP_FFTWMPI)"
#endif
end subroutine

subroutine H_internal_single(H,H_out,isite,dim_lat,N_rep,nmag)
    integer,intent(in)          :: isite
    integer,intent(in)          :: nmag
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(in)          :: H(3,Nmag,N_rep(1),N_rep(2),N_rep(3))
    real(8),intent(inout)       :: H_out(3)

!    integer     :: div(4),modu(4)
    integer     :: i4(4),i

    ! does not work
!    modu=[nmag,nmag*dim_lat(1),nmag*product(dim_lat(:2)),nmag*product(dim_lat)]
!    div=[(product(modu(:i-1)),i=1,4)]
!    i4=isite-1
!    i4=i4/div
!    i4=modulo(i4,modu)+1

    i4(4)=(isite-1)/(nmag*product(dim_lat(:2)))+1

    i=modulo((isite-1),nmag*product(dim_lat(:2)))+1
    i4(3)=(i-1)/(nmag*dim_lat(1))+1

    i=modulo(i-1,nmag*dim_lat(1))+1
    i4(2)=(i-1)/nmag+1

    i4(1)=modulo((isite-1),nmag)+1

    H_out=H(:,i4(1),i4(2),i4(3),i4(4))
end subroutine

end module
