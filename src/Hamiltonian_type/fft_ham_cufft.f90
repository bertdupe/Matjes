module m_fft_H_cufft
#ifdef CPP_CUDA
!module which contains the general discrete fourier transform Hamiltonian based on cuFFT

!the type is used by first calling init_shape, followed by init_op with the operator as described in more detail for the dipolar_fft interaction

use,intrinsic :: ISO_FORTRAN_ENV, only: error_unit
use, intrinsic :: iso_c_binding
use m_fft_H_base, only: fft_H
use m_type_lattice,only: lattice
use m_fft_H_internal, only: int_set_M, int_get_H
use m_H_type, only: len_desc
use m_cuda_fft
private
public fft_H_cufft

type,extends(fft_H) ::  fft_H_cufft
    private
    type(c_ptr)     :: plan_fwd=c_null_ptr   !cufftHandle which contains the handle
    type(c_ptr)     :: plan_bwd=c_null_ptr   !cufftHandle which contains the handle

    !arrays on device
    type(c_ptr)     :: M=  c_null_ptr   !mode in normal-space
    type(c_ptr)     :: M_F=c_null_ptr   !mode in fourier-space
    type(c_ptr)     :: K_F=c_null_ptr   !demagnetization tensor in fourier-space
    type(c_ptr)     :: H=  c_null_ptr   !effective field fourier-space
    type(c_ptr)     :: H_F=c_null_ptr   !effective field normal-space
contains
    !evaluation routines
    procedure,public    :: get_H            !get effective field
    procedure,public    :: get_H_single     !get effective field for a single site

    !initialization routines
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


subroutine init_shape(this,abbrev,dim_mode,periodic,dim_lat,Kbd,N_rep)
    !initializes the arrays with whose the work will be done, sets the shapes, and returns shape data necessary for construction of the operator tensor
    class(fft_H_cufft),intent(inout)  :: this
    character(1),intent(in)     :: abbrev       !Abbreviation for order parameter (M,E,U,...) according to order_parameter_abbrev 
    integer,intent(in)          :: dim_mode
    logical,intent(in)          :: periodic(3)  !T: periodic boundary, F: open boundary
    integer,intent(in)          :: dim_lat(3)
    integer,intent(out)         :: Kbd(2,3) 
    integer,intent(out)         :: N_rep(3)

#ifdef CPP_CUDA
    if(c_associated(this%M))then
        write(error_unit,'(///A)') "Trying to initialize shape of fft_H, but the fourier-transform has already been set"
        write(error_unit,'(A)') "init_shape shall only be called once"
        ERROR STOP
    endif

    Call this%fft_H%init_shape(abbrev,dim_mode,periodic,dim_lat,Kbd,N_rep)

    !set order work arrays and fourier transform
    Call cuda_fft_init(this%dim_mode,this%N_rep,this%M,this%M_F,this%H,this%H_F)
    Call set_cuda_fft_plan(this%dim_mode,this%N_rep,this%plan_fwd,this%plan_bwd)
#else
        ERROR STOP "CANNOT USE CUDA DIPOL (mag_dip_fft) without CUDA (CPP_CUDA)"
#endif
end subroutine

subroutine add(this,H_in)
    class(fft_H_cufft),intent(inout)  :: this
    class(fft_H),intent(in)     :: H_in

    logical :: is_set

    if(.not.H_in%is_set()) ERROR STOP "CANNOT ADD fft_H as H_in is not set"
    is_set=this%is_set()
    if(is_set)then
        if(.not.this%same_space(H_in)) ERROR STOP "CANNOT ADD fft_H as this and H_in do not act on same space"
        Call this%fft_H%add(H_in)
        select type(H_in)
        type is(fft_H_cufft)
            Call cuda_fft_add_cmplx(size(this%M_n),this%K_F,H_in%K_F) 
        class default
            ERROR STOP "HAS TO HAVE SAME TYPE"
        end select
    else
        Call H_in%copy(this)
    endif
end subroutine

subroutine bcast_fft(this,comm)
    use mpi_util
    class(fft_H_cufft),intent(inout)    ::  this
    type(mpi_type),intent(in)           ::  comm

    Call this%fft_H%bcast(comm)

    ERROR STOP "IMPLEMENT BCAST fft_H_cufft"    !parallelization with MPI and CUDA might me more tricky...
!    if(comm%ismas)then
!        shp=shape(this%M_n)
!    endif
!    Call bcast(shp,comm)
!    if(.not.comm%ismas)then
!        Call this%init_internal(this%periodic)
!        allocate(this%M_n(shp(1),shp(2))) 
!        allocate(this%M_F(shp(1),shp(2))) 
!        allocate(this%_n(shp(1),shp(2))) 
!        allocate(this%H_F(shp(1),shp(2))) 
!    endif
!    Call bcast_alloc(this%K_F,comm)
!    if(.not.comm%ismas) Call set_fftw_plans(this)
end subroutine

subroutine mv(this,H_out)
    class(fft_H_cufft),intent(inout)    :: this
    class(fft_H),intent(inout)          :: H_out

    select type(H_out)
    type is(fft_H_cufft)
        H_out%plan_fwd=this%plan_fwd
        H_out%plan_bwd=this%plan_bwd
        H_out%M       =this%M
        H_out%M_F     =this%M_F
        H_out%K_F     =this%K_F
        H_out%H       =this%H
        H_out%H_F     =this%H_F

        this%plan_fwd=c_null_ptr
        this%plan_bwd=c_null_ptr
        this%M       =c_null_ptr
        this%M_F     =c_null_ptr
        this%K_F     =c_null_ptr
        this%H       =c_null_ptr
        this%H_F     =c_null_ptr

        ERROR STOP "test if this works"
    class default
        ERROR STOP "HAS TO HAVE SAME TYPE"
    end select
    Call this%destroy()
end subroutine

subroutine copy(this,H_out)
    class(fft_H_cufft),intent(in):: this
    class(fft_H),intent(inout)   :: H_out

    integer(C_int)  :: size_k

    Call this%fft_H%copy(H_out)
    select type(H_out)
    type is(fft_H_cufft)
        Call cuda_fft_init(H_out%dim_mode,H_out%N_rep,H_out%M,H_out%M_F,H_out%H,H_out%H_F)
        size_k=int(H_out%dim_mode**2*product(H_out%N_rep),C_int)
        Call cuda_fft_copy_cmplx(size_K,this%K_F,H_out%K_F)
        Call set_cuda_fft_plan(H_out%dim_mode,H_out%N_rep,H_out%plan_fwd,H_out%plan_bwd)
    class default
        ERROR STOP "HAS TO HAVE SAME TYPE"
    end select
end subroutine

subroutine destroy(this)
    class(fft_H_cufft),intent(inout)     :: this

    Call this%fft_H%destroy()
#ifdef CPP_CUDA
    if(c_associated(this%plan_fwd)) Call cuda_fft_destroy_plan(this%plan_fwd)
    if(c_associated(this%plan_bwd)) Call cuda_fft_destroy_plan(this%plan_bwd)

    if(c_associated(this%M)) Call cuda_fft_destroy_real(this%M)
    if(c_associated(this%H)) Call cuda_fft_destroy_real(this%H)

    if(c_associated(this%M_F)) Call cuda_fft_destroy_real(this%M_F)
    if(c_associated(this%K_F)) Call cuda_fft_destroy_real(this%K_F)
    if(c_associated(this%H_F)) Call cuda_fft_destroy_real(this%H_F)
#else
    ERROR STOP "Requires CPP_CUDA"
#endif
end subroutine

subroutine init_op(this,dim_mode,K_n,desc_in)
    !subroutine which initializes the fourier-transformed operator of K, while deallocating K_N
    class(fft_H_cufft),intent(inout)         :: this
    integer,intent(in)                      :: dim_mode
    real(8),intent(inout),allocatable       :: K_n(:,:,:)
    character(len=*),intent(in),optional    :: desc_in
    integer(C_int)  ::  length

    integer(C_int)  :: shape_test(3)

    !some initial tests
    length=int(size(K_n),C_int)
    shape_test=[this%dim_mode, this%dim_mode, product(this%N_rep)]
    if(any(shape_test==0))then
        write(error_unit,'(///A)') "Nonsensical expected shape for operator input"
        write(error_unit,'(A,3I8)') "shape expectation: ", shape_test
        write(error_unit,'(A)') "Maybe init_shape has not been called?"
        ERROR STOP
    endif
    if(any(shape_test/=shape(K_N)))then
        write(error_unit,'(///A)') "Shape of input array for setting operator in fft_H_cufft seems to be wrong"
        write(error_unit,'(A,3I8)') "shape expectation: ", shape_test
        write(error_unit,'(A,3I8)') "shape input array: ", shape(K_N)
        write(error_unit,'(A)') "Maybe init_shape has not been called?"
        ERROR STOP
    endif
    if(.not.c_associated(this%plan_fwd))then
        write(error_unit,'(///A)') "Fourier transformation plan is not set when tying to initialize operator in fft_H_cufft-type"
        write(error_unit,'(A)') "This is most probably the case because init_shape has not been called previously"
        ERROR STOP
    endif
    if(.not.c_associated(this%M))then
        write(error_unit,'(///A)') "cuda magnetic moment field is not set when tying to initialize operator in fft_H_cufft-type"
        write(error_unit,'(A)') "This is most probably the case because init_shape has not been called previously"
        ERROR STOP
    endif
#ifdef CPP_CUDA

    !actual part which sets description and operator in fourier space (this%K_F)
    K_N=K_N/product(this%N_rep)
    Call this%fft_H%init_op(dim_mode,K_n,desc_in)
    Call cuda_fft_set_operator(this%dim_mode,this%N_rep,K_N,this%K_F)
    this%set=.true.

    !deallocate K_n, since at some point one might want to keep it in here
    deallocate(K_n)
#else
        ERROR STOP "CANNOT USE FFTW-Hamiltonian without FFTW (CPP_FFTW3)"
#endif
end subroutine

subroutine set_M(this,lat)
    !set this%M_n from lat%M%modes according to this%M_internal
    class(fft_H_cufft),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
#ifdef CPP_CUDA
    real(8),pointer,contiguous  :: M(:)
    

    integer(C_int)  ::  length
    integer         :: mode_id

    length=int(size(this%m_n),C_int)
    call this%get_mode_id(mode_id)
    select case (mode_id)
       case (1)
         Call this%M_internal(this%M_n,lat%M%modes,lat%dim_lat,this%N_rep,size(this%M_n,1))
       case (5)
         Call this%M_internal(this%M_n,lat%U%modes,lat%dim_lat,this%N_rep,size(this%M_n,1))
       case default
         stop 'node coded in set_M'
    end select
    Call cuda_fft_set_real(length,this%M_n,this%m)
#else
    ERROR STOP "fft_H_cufft%set_M requires CPP_CUDA"
#endif
end subroutine

subroutine get_H(this,lat,Hout)
    !get effective field H 
    class(fft_H_cufft),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    real(8),intent(inout)         ::  Hout(:,:)
#ifdef CPP_CUDA
    integer(C_int)  ::  length

    length=int(size(this%m_n),C_int)

    Call this%set_M(lat)

    Call cuda_fft_calc_H(this%dim_mode,this%N_rep,this%m,this%m_f,this%k_f,this%h,this%h_f,this%plan_fwd,this%plan_bwd)

    Call cuda_fft_get_real(length,this%h,this%H_n)
    Call this%H_internal(this%H_n,Hout,lat%dim_lat,this%N_rep,size(Hout,1))
#else
    ERROR STOP "fft_H_cufft%get_H requires CPP_CUDA"
#endif
end subroutine

subroutine get_H_single(this,lat,site,Hout)
    !get effective field H
    class(fft_H_cufft),intent(inout)    ::  this
    type(lattice),intent(in)            ::  lat
    integer,intent(in)                  ::  site
    real(8),intent(inout)               ::  Hout(3)

#ifdef CPP_CUDA
    integer(C_int)  ::  length

    length=int(size(this%m_n),C_int)
    Call this%set_M(lat)

    Call cuda_fft_calc_H(this%dim_mode,this%N_rep,this%m,this%m_f,this%k_f,this%h,this%h_f,this%plan_fwd,this%plan_bwd)

    Call cuda_fft_get_real(length,this%h,this%H_n)
    call this%get_mode_id(mode_id)
    select case (mode_id)
       case (1)
         Call H_internal_single(this%H_n,Hout,site,lat%dim_lat,this%N_rep,lat%nmag)
       case (5)
         Call H_internal_single(this%H_n,Hout,site,lat%dim_lat,This%N_rep,lat%nph)
       case default
         stop 'node coded in get_H_single'
    end select
#else
    ERROR STOP "fft_H_cufft%get_H_single requires CPP_CUDA"
#endif
end subroutine

subroutine H_internal_single(H,H_out,isite,dim_lat,N_rep,nmag)
    integer,intent(in)          :: isite
    integer,intent(in)          :: nmag
    integer,intent(in)          :: dim_lat(3)
    integer,intent(in)          :: N_rep(3)
    real(8),intent(in)          :: H(3,Nmag,N_rep(1),N_rep(2),N_rep(3))
    real(8),intent(inout)       :: H_out(3)

    integer     :: div(4),modu(4)
    integer     :: i4(4),i 

    modu=[nmag,nmag*dim_lat(1),nmag*product(dim_lat(:2)),nmag*product(dim_lat)]
    div=[(product(modu(:i-1)),i=1,4)]
    i4=isite-1
    i4=i4/div
    i4=modulo(i4,modu)+1

    H_out=H(:,i4(1),i4(2),i4(3),i4(4))
end subroutine
#endif
end module
