module m_fft_H_cufft
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


    !arrays on host
    real(C_DOUBLE),allocatable  ::  Mh_n(:,:)    !magnetization in normal-space
    real(C_DOUBLE),allocatable  ::  Hh_n(:,:)    !magnetization in normal-space

    !arrays on device
    type(c_ptr)     :: M_n=c_null_ptr   !magnetization in normal-space
    type(c_ptr)     :: M_F=c_null_ptr   !magnetization in fourier-space
    type(c_ptr)     :: K_F=c_null_ptr   !demagnetization tensor in fourier-space
    type(c_ptr)     :: H_N=c_null_ptr   !effective field fourier-space
    type(c_ptr)     :: H_F=c_null_ptr   !effective field normal-space

    procedure(int_set_M), pointer,nopass    ::  M_internal => null()    !function to set internal magnetization depending on periodic/open boundaries
    procedure(int_get_H), pointer,nopass    ::  H_internal => null()    !function to get effective field depending on periodic/open boundaries
contains
    procedure,public    :: get_H            !get effective field
!    procedure,public    :: get_H_single     !get effective field for a single site
    procedure,public    :: get_E            !get energy
!    procedure,public    :: get_E_single     !get energy for a single site
!    procedure,public    :: get_E_distrib    !get energy-distribution in Nmag*Ncell-space
    procedure,public    :: init_shape       !initializes the shape and the H and M operators, arrays
    procedure,public    :: init_op          !initializes the K_F array
    
!    !utility functions
!    procedure,public    :: destroy          !destroys all data of this type
    procedure,public    :: copy             !copy this to new instance
!    procedure,public    :: mv               !mv this to new instance
    procedure,public    :: add              !adds two fft_H in same space by adding their K_F
!    procedure,public    :: bcast=>bcast_fft !mv this to new instance
!
!
!    !internal procedures
    procedure           :: set_M                !set internal magnetization in normal-space from lattice
    procedure           :: init_internal        !initialize internal procedures
end type
contains 


subroutine init_shape(this,dim_mode,periodic,dim_lat,Kbd,N_rep)
    !initializes the arrays with whose the work will be done, sets the shapes, and returns shape data necessary for construction of the operator tensor
    class(fft_H_cufft),intent(inout)  :: this
    integer,intent(in)          :: dim_mode
    logical,intent(in)          :: periodic(3)  !T: periodic boundary, F: open boundary
    integer,intent(in)          :: dim_lat(3)
    integer,intent(out)         :: Kbd(2,3) 
    integer,intent(out)         :: N_rep(3)

#ifdef CPP_CUDA
    integer         :: i
    integer(C_int)  :: Nk_tot           !number of state considered in FT (product of N_rep)


    if(c_associated(this%M_n))then
        write(error_unit,'(///A)') "Trying to initialize shape of fft_H, but the fourier-transform has already been set"
        write(error_unit,'(A)') "init_shape shall only be called once"
        ERROR STOP
    endif

    Call this%fft_H%init_shape(dim_mode,periodic,dim_lat,Kbd,N_rep)

    Nk_tot=product(N_rep)
    !set order work arrays and fourier transform
    allocate(this%Mh_N(dim_mode,Nk_tot),source=0.0d0)
    allocate(this%Hh_N(dim_mode,Nk_tot),source=0.0d0)
    Call cuda_fft_init(this%dim_mode,this%N_rep,this%M_N,this%M_F,this%H_N,this%H_F)
    Call set_cuda_fft_plan(this%dim_mode,this%N_rep,this%plan_fwd,this%plan_bwd)
    Call this%init_internal(this%periodic)
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
        ERROR STOP "IMPLEMENT"
!        if(.not.this%same_space(H_in)) ERROR STOP "CANNOT ADD fft_H as this and H_in do not act on same space"
!        Call this%fft_H%add(H_in)
!        select type(H_in)
!        type is(fft_H_fftw)
!            this%K_F=this%K_F+H_in%K_F
!        class default
!            ERROR STOP "HAS TO HAVE SAME TYPE"
!        end select
    else
        Call H_in%copy(this)
    endif
end subroutine
!
!subroutine bcast_fft(this,comm)
!    use mpi_util
!    class(fft_H_fftw),intent(inout)  ::  this
!    type(mpi_type),intent(in)   ::  comm
!
!    integer ::  shp(2)
!
!    Call this%fft_H%bcast(comm)
!
!    if(comm%ismas)then
!        shp=shape(this%M_n)
!    endif
!    Call bcast(shp,comm)
!    if(.not.comm%ismas)then
!        Call this%init_internal(this%periodic)
!        allocate(this%M_n(shp(1),shp(2))) 
!        allocate(this%M_F(shp(1),shp(2))) 
!        allocate(this%H_n(shp(1),shp(2))) 
!        allocate(this%H_F(shp(1),shp(2))) 
!    endif
!    Call bcast_alloc(this%K_F,comm)
!    if(.not.comm%ismas) Call set_fftw_plans(this)
!end subroutine
!
!subroutine mv(this,H_out)
!    class(fft_H_fftw),intent(inout)  :: this
!    class(fft_H),intent(inout)   :: H_out
!
!    if(.not.this%set)then
!        ERROR STOP "CANNOT MV UNINITIALIZED FFT_H"
!    endif
!    Call this%fft_H%mv(H_out)
!    select type(H_out)
!    type is(fft_H_fftw)
!        call move_alloc(this%M_n,H_out%M_n)
!        call move_alloc(this%M_F,H_out%M_F)
!        call move_alloc(this%H_n,H_out%H_n)
!        call move_alloc(this%H_F,H_out%H_F)
!        call move_alloc(this%K_F,H_out%K_F)
!        H_out%M_internal=>this%M_internal
!        H_out%H_internal=>this%H_internal
!    class default
!        ERROR STOP "HAS TO HAVE SAME TYPE"
!    end select
!    Call this%destroy()
!end subroutine
!
subroutine copy(this,H_out)
    class(fft_H_cufft),intent(in):: this
    class(fft_H),intent(inout)   :: H_out

    integer(C_int)  :: Nk_tot           !number of state considered in FT (product of N_rep)
    integer(C_int)  :: size_k

    Call this%fft_H%copy(H_out)
    select type(H_out)
    type is(fft_H_cufft)
        Call cuda_fft_init(H_out%dim_mode,H_out%N_rep,H_out%M_N,H_out%M_F,H_out%H_N,H_out%H_F)
        allocate(H_out%Mh_n,mold=this%Mh_N)
        allocate(H_out%Hh_n,mold=this%Hh_N)
        size_k=int(H_out%dim_mode**2*product(H_out%N_rep),C_int)
        Call cuda_fft_copy_cmplx(size_K,this%K_F,H_out%K_F)
        Call set_cuda_fft_plan(H_out%dim_mode,H_out%N_rep,H_out%plan_fwd,H_out%plan_bwd)
        Call H_out%init_internal(H_out%periodic)
    class default
        ERROR STOP "HAS TO HAVE SAME TYPE"
    end select
end subroutine
!
!subroutine destroy(this)
!    class(fft_H_fftw),intent(inout)     :: this
!
!    Call this%fft_H%destroy()
!#ifdef CPP_FFTW3
!    if(c_associated(this%plan_mag_f)) Call fftw_destroy_plan(this%plan_mag_f)
!    if(c_associated(this%plan_H_I))   Call fftw_destroy_plan(this%plan_H_I)
!#else
!    ERROR STOP "Requires FFTW3"
!#endif
!    if(allocated(this%M_n)) deallocate(this%M_n)
!    if(allocated(this%M_F)) deallocate(this%M_F)
!    if(allocated(this%K_F)) deallocate(this%K_F)
!    if(allocated(this%H_F)) deallocate(this%H_F)
!    if(allocated(this%H_n)) deallocate(this%H_n)
!    nullify(this%M_internal, this%M_internal)
!end subroutine
!
subroutine init_op(this,dim_mode,K_n,desc_in)
    !subroutine which initializes the fourier-transformed operator of K, while deallocating K_N
    class(fft_H_cufft),intent(inout)         :: this
    integer,intent(in)                      :: dim_mode
    real(8),intent(inout),allocatable       :: K_n(:,:,:)
    character(len=*),intent(in),optional    :: desc_in
    integer(C_int)  ::  length

    integer(C_int)  :: shape_test(3)
    integer     :: i

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
    if(.not.c_associated(this%M_n))then
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



subroutine init_internal(this,periodic)
    !initialize internal procedures M_internal and H_internal
    use m_fft_H_internal
    class(fft_H_cufft),intent(inout)    :: this
    logical,intent(in)            :: periodic(3)  !T: periodic boundary, F: open boundary

    if(all(periodic))then
        this%M_internal=>set_M_period_TTT
        this%H_internal=>set_H_period_TTT
    elseif(all(periodic(1:2)))then
        this%M_internal=>set_M_period_TTF
        this%H_internal=>set_H_period_TTF
    elseif(all(periodic(1:1)))then
        this%M_internal=>set_M_period_TFF
        this%H_internal=>set_H_period_TFF
    else
        this%M_internal=>set_M_period_FFF
        this%H_internal=>set_H_period_FFF
    endif
end subroutine

subroutine set_M(this,lat)
    !set this%M_n from lat%M%modes according to this%M_internal
    class(fft_H_cufft),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat

    integer(C_int)  ::  length

    length=int(size(this%mh_n),C_int)
    Call this%M_internal(this%MH_n,lat%M%modes,lat%dim_lat,this%N_rep,size(this%Mh_n,1))
    Call cuda_fft_set_real(length,this%MH_n,this%m_n)
end subroutine
!
subroutine get_H(this,lat,Hout)
    !get effective field H 
    class(fft_H_cufft),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    real(8),intent(inout)         ::  Hout(:,:)
#ifdef CPP_CUDA
    integer ::  i,j,l
    integer(C_int)  ::  length

    length=int(size(this%mh_n),C_int)

    Call this%set_M(lat)

    Call cuda_fft_calc_H(this%dim_mode,this%N_rep,this%m_n,this%m_f,this%k_f,this%h_n,this%h_f,this%plan_fwd,this%plan_bwd)

    do i=1,size(this%Hh_n,2)
        this%Hh_n(:,i)=real(i,8)
    enddo
    Call cuda_fft_get_real(length,this%h_n,this%Hh_n)
    Call this%H_internal(this%Hh_n,Hout,lat%dim_lat,this%N_rep,size(Hout,1))
#else
    ERROR STOP "fft_H_cufft%get_H requires CPP_CUDA"
#endif
end subroutine
!
!subroutine get_H_single(this,lat,site,Hout)
!    !get effective field H
!    class(fft_H_fftw),intent(inout)    ::  this
!    type(lattice),intent(in)      ::  lat
!    integer,intent(in)            ::  site
!    real(8),intent(inout)         ::  Hout(3)
!#ifdef CPP_FFTW3
!    integer ::  i,j,l
!
!    Call this%set_M(lat)
!    Call fftw_execute_dft_r2c(this%plan_mag_F, this%M_n, this%M_F)
!
!    this%H_F=cmplx(0.0d0,0.0d0,8)
!    do j=1,size(this%M_F,2)
!        do i=1,lat%Nmag*3
!            do l=1,lat%Nmag*3
!                this%H_F(i,j)=this%H_F(i,j)+this%K_F(i,l,j)*this%M_F(l,j)   !change to matmul
!            enddo
!        enddo
!    enddo
!    Call fftw_execute_dft_c2r(this%plan_H_I, this%H_F, this%H_n)
!    Call H_internal_single(this%H_n,Hout,site,lat%dim_lat,this%N_rep,lat%nmag)
!#else
!    ERROR STOP "fft_H%get_H_single requires CPP_FFTW"
!#endif
!end subroutine
!
!subroutine get_E_distrib(this,lat,Htmp,E)
!    class(fft_H_fftw),intent(inout)    ::  this
!    type(lattice),intent(in)      ::  lat
!    real(8),intent(inout)         ::  Htmp(:,:)
!    real(8),intent(out)           ::  E(:)
!
!    Call this%get_H(lat,Htmp)
!    Htmp=Htmp*lat%M%modes_v
!    E=sum(reshape(Htmp,[3,lat%Nmag*lat%Ncell]),1)*2.0d0 !not sure about *2.0d0
!end subroutine
!
subroutine get_E(this,lat,Htmp,E)
    class(fft_H_cufft),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    real(8),intent(inout)         ::  Htmp(:,:)
    real(8),intent(out)           ::  E

    Call this%get_H(lat,Htmp)
    Htmp=Htmp*lat%M%modes_v
    E=sum(Htmp)
end subroutine
!
!subroutine get_E_single(this,lat,site,E)
!    !get energy by getting the effective field caused by all sites and multiply the considered site's H with the moment there
!    !it might be faster to do the explicit folding in real space to get the effective field caused by the considered site and multiply then with the entire magnetization (might be worth testing)
!    class(fft_H_fftw),intent(inout)  :: this
!    type(lattice),intent(in)    :: lat
!    integer,intent(in)          :: site
!    real(8),intent(out)         :: E
!
!    real(8)                     :: Htmp(3)
!
!    Call this%get_H_single(lat,site,Htmp)
!    Htmp=Htmp*lat%M%modes_3(:,site)
!    E=sum(Htmp)*2.0d0
!end subroutine
!
!subroutine set_fftw_plans(this)
!    class(fft_H_fftw),intent(inout)  :: this
!
!    integer(C_INT)  :: N_rep_rev(3)     !reversed N_rep necessary for fftw3 (col-major -> row-major)
!    integer(C_int)  :: howmany          !dimension of quantitiy which is fourier-transformed (see FFTW3)
!    
!#ifdef CPP_FFTW3
!    if(.not.allocated(this%M_n)) ERROR STOP "cannot set fftw_plans as M_n not allocated"
!    if(.not.allocated(this%M_F)) ERROR STOP "cannot set fftw_plans as M_F not allocated"
!    if(.not.allocated(this%H_n)) ERROR STOP "cannot set fftw_plans as H_n not allocated"
!    if(.not.allocated(this%H_F)) ERROR STOP "cannot set fftw_plans as H_F not allocated"
!
!    N_rep_rev=this%N_rep(size(this%N_rep):1:-1)
!    howmany=int(size(this%M_n,1),C_int)
!    this%plan_mag_F= fftw_plan_many_dft_r2c(int(3,C_INT), N_rep_rev, howmany,&
!                                           &this%M_n,     N_rep_rev,&
!                                           &howmany,      int(1,C_int), &
!                                           &this%M_F,     N_rep_rev,&
!                                           &howmany,      int(1,C_int), &
!                                           &FFTW_FORWARD+FFTW_MEASURE+FFTW_PATIENT)
!
!    this%plan_H_I= fftw_plan_many_dft_c2r(int(3,C_INT), N_rep_rev, howmany,&
!                                         &this%H_F,     N_rep_rev,&
!                                         &howmany,      int(1,C_int), &
!                                         &this%H_n,     N_rep_rev,&
!                                         &howmany,      int(1,C_int), &
!                                         &FFTW_BACKWARD+FFTW_MEASURE+FFTW_PATIENT)
!#else
!    ERROR STOP "CANNOT USE FFT_H without FFTW (CPP_FFTW3)"
!#endif
!end subroutine
!
!subroutine H_internal_single(H,H_out,isite,dim_lat,N_rep,nmag)
!    integer,intent(in)          :: isite
!    integer,intent(in)          :: nmag
!    integer,intent(in)          :: dim_lat(3)
!    integer,intent(in)          :: N_rep(3)
!    real(8),intent(in)          :: H(3,Nmag,N_rep(1),N_rep(2),N_rep(3))
!    real(8),intent(inout)       :: H_out(3)
!
!    integer     :: div(4),modu(4)
!    integer     :: i4(4),i 
!
!    modu=[nmag,nmag*dim_lat(1),nmag*product(dim_lat(:2)),nmag*product(dim_lat)]
!    div=[(product(modu(:i-1)),i=1,4)]
!    i4=isite-1
!    i4=i4/div
!    i4=modulo(i4,modu)+1
!
!    H_out=H(:,i4(1),i4(2),i4(3),i4(4))
!end subroutine

end module
