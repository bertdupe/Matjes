module m_fft_ham
!module which contains the general discrete fourier transform Hamiltonian based on FFTW3

!the type is used by first calling init_shape, followed by init_op with the operator as described in more detail for the dipolar_fft interaction

use,intrinsic :: ISO_FORTRAN_ENV, only: error_unit
use m_fftw3
use m_type_lattice,only: lattice
use m_fft_H_internal, only: int_set_M, int_get_H
private
public fft_H

type        ::  fft_H
    private
    logical         :: set=.false.  !all parameters have been initialized
    integer         :: N_rep(3)=0   !size in each dimension of each considered field
    logical         :: periodic(3)=.true.   !consider as periodic or open boundary condition along each direction (T:period, F:open)
                                            ! (dim_lat(i)=1->period(i)=T, since the calculation in the periodic case is easier, but choice of supercell_vec still does not consider periodicity)
    type(c_ptr)     :: plan_mag_F=c_null_ptr   !FFTW plan M_n -> M_F
    type(c_ptr)     :: plan_H_I=c_null_ptr     !FFTW plan H_F -> H_n
    real(C_DOUBLE),allocatable              ::  M_n(:,:)    !magnetization in normal-space
    complex(C_DOUBLE_COMPLEX),allocatable   ::  M_F(:,:)    !magnetization in fourier-space

    complex(C_DOUBLE_COMPLEX),allocatable   ::  K_F(:,:,:)  !demagnetization tensor in fourier-space

    complex(C_DOUBLE_COMPLEX),allocatable   ::  H_F(:,:)    !effective field fourier-space 
    real(C_DOUBLE),allocatable              ::  H_n(:,:)    !effective field normal-space
    procedure(int_set_M), pointer,nopass    ::  M_internal => null()    !function to set internal magnetization depending on periodic/open boundaries
    procedure(int_get_H), pointer,nopass    ::  H_internal => null()    !function to get effective field depending on periodic/open boundaries
contains
    procedure,public    :: get_H            !get effective field
    procedure,public    :: get_E            !get energy
    procedure,public    :: get_E_distrib    !get energy-distribution in Nmag*Ncell-space
    procedure,public    :: is_set           !returns set
    procedure,public    :: init_shape       !initializes the shape and the H and M operators, arrays
    procedure,public    :: init_op          !initializes the K_F array
    
    !utility functions
    procedure,public    :: destroy          !destroys all data of this type
    procedure,public    :: copy             !copy this to new instance
    procedure,public    :: mv               !mv this to new instance
    procedure,public    :: bcast=>bcast_fft !mv this to new instance


    !internal procedures
    procedure           :: set_M            !set internal magnetization in normal-space from lattice
    procedure           :: init_internal    !initialize internal procedures

end type
contains 

subroutine bcast_fft(this,comm)
    use mpi_util
    class(fft_H),intent(inout)  ::  this
    type(mpi_type),intent(in)   ::  comm

    integer ::  shp(2)

    Call bcast(this%N_rep,comm)
    Call bcast(this%periodic,comm)
    if(comm%ismas)then
        shp=shape(this%M_n)
    endif
    Call bcast(shp,comm)
    if(.not.comm%ismas)then
        Call this%init_internal(this%periodic)
        allocate(this%M_n(shp(1),shp(2))) 
        allocate(this%M_F(shp(1),shp(2))) 
        allocate(this%H_n(shp(1),shp(2))) 
        allocate(this%H_F(shp(1),shp(2))) 
    endif
    Call bcast_alloc(this%K_F,comm)
    if(.not.comm%ismas) Call set_fftw_plans(this)
    Call bcast(this%set,comm)
end subroutine

subroutine mv(this,H_out)
    class(fft_H),intent(inout)  :: this
    type(fft_H),intent(inout)   :: H_out

    if(.not.this%set)then
        ERROR STOP "CANNOT MV UNINITIALIZED FFT_H"
    endif

    H_out%N_rep=this%N_rep
    H_out%periodic=this%periodic
    call move_alloc(this%M_n,H_out%M_n)
    call move_alloc(this%M_F,H_out%M_F)
    call move_alloc(this%H_n,H_out%H_n)
    call move_alloc(this%H_F,H_out%H_F)
    call move_alloc(this%K_F,H_out%K_F)
    H_out%M_internal=>this%M_internal
    H_out%H_internal=>this%H_internal

    H_out%set=this%set
    Call this%destroy()
end subroutine

subroutine copy(this,H_out)
    class(fft_H),intent(in)     :: this
    type(fft_H),intent(inout)   :: H_out

    if(.not.this%set)then
        ERROR STOP "CANNOT COPY UNINITIALIZED FFT_H"
    endif
    H_out%N_rep=this%N_rep
    H_out%periodic=this%periodic
    allocate(H_out%M_n,mold=this%M_n)
    allocate(H_out%M_F,mold=this%M_F)
    allocate(H_out%H_n,mold=this%H_n)
    allocate(H_out%H_F,mold=this%H_F)
    allocate(H_out%K_F,source=this%K_F)
    H_out%M_internal=>this%M_internal
    H_out%H_internal=>this%H_internal
    Call set_fftw_plans(H_out)

    H_out%set=this%set
end subroutine

subroutine destroy(this)
    class(fft_H),intent(inout)     :: this

    this%N_rep=0
    this%set=.false.
    this%periodic=.true.
#ifdef CPP_FFTW3
    if(c_associated(this%plan_mag_f)) Call fftw_destroy_plan(this%plan_mag_f)
    if(c_associated(this%plan_H_I))   Call fftw_destroy_plan(this%plan_H_I)
#else
    ERROR STOP "Requires FFTW3"
#endif
    if(allocated(this%M_n)) deallocate(this%M_n)
    if(allocated(this%M_F)) deallocate(this%M_F)
    if(allocated(this%K_F)) deallocate(this%K_F)
    if(allocated(this%H_F)) deallocate(this%H_F)
    if(allocated(this%H_n)) deallocate(this%H_n)
    nullify(this%M_internal, this%M_internal)

end subroutine

pure function is_set(this)result(set)
    class(fft_H),intent(in)       :: this
    logical     ::  set
    set=this%set
end function

subroutine init_op(this,dim_mode,K_n)
    !subroutine which initializes the fourier-transformed operator of K, while deallocating K_N
    class(fft_H),intent(inout)          :: this
    integer,intent(in)                  :: dim_mode
    real(8),intent(inout),allocatable   :: K_n(:,:,:)

#ifdef CPP_FFTW3
    integer         :: Nk_tot           !number of state considered in FT (product of N_rep)
    integer(C_INT)  :: N_rep_rev(3)     !reversed N_rep necessary for fftw3 (col-major -> row-major)
    integer(C_int)  :: howmany          !dimension of quantitiy which is fourier-transformed (see FFTW3)
    type(c_ptr)     :: plan_K_F !plan for fourier transformation of K

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
        ERROR STOP "CANNOT USE FFTW-Hamiltonian without FFTW (CPP_FFTW3)"
#endif
end subroutine

subroutine init_shape(this,dim_mode,periodic,dim_lat,Kbd,N_rep)
    !initializes the arrays with whose the work will be done, sets the shapes, and returns shape data necessary for construction of the operator tensor
!$  use omp_lib    
    class(fft_H),intent(inout)  :: this
    integer,intent(in)          :: dim_mode
    logical,intent(in)          :: periodic(3)  !T: periodic boundary, F: open boundary
    integer,intent(in)          :: dim_lat(3)
    integer,intent(out)         :: Kbd(2,3) 
    integer,intent(out)         :: N_rep(3)

#ifdef CPP_FFTW3
    integer         :: i
    integer         :: Nk_tot           !number of state considered in FT (product of N_rep)


    if(allocated(this%M_n))then
        write(error_unit,'(///A)') "Trying to initialize shape of fft_H, but the fourier-transform has already been set"
        write(error_unit,'(A)') "init_shape shall only be called once"
        ERROR STOP
    endif

    N_rep=dim_lat
    !set K-boundaries for periodic boundaries
    Kbd(1,:)=[0,0,0]
    Kbd(2,:)=dim_lat-1
    !set K-boundaries for open boundaries if open
    do i=1,3
        if(.not.periodic(i))then
            N_rep(i)=2*N_rep(i)
            Kbd(:,i)=[-dim_lat(i)+1,dim_lat(i)-1]
        endif
    enddo
    this%N_rep=N_rep
    this%periodic=periodic
    Nk_tot=product(N_rep)

!$  Call fftw_plan_with_nthreads(omp_get_max_threads())
    !set order work arrays and fourier transform
    allocate(this%M_N(dim_mode,Nk_tot),source=0.0d0)
    allocate(this%M_F(dim_mode,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
    allocate(this%H_F(dim_mode,Nk_tot),source=cmplx(0.0d0,0.0d0,8))
    allocate(this%H_n(dim_mode,Nk_tot),source=0.0d0)
    Call set_fftw_plans(this)

    Call this%init_internal(periodic)
#else
        ERROR STOP "CANNOT USE FFTW DIPOL (mag_dip_fft) without FFTW (CPP_FFTW3)"
#endif
end subroutine


subroutine init_internal(this,periodic)
    !initialize internal procedures M_internal and H_internal
    use m_fft_H_internal
    class(fft_H),intent(inout)    :: this
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
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat

    Call this%M_internal(this%M_n,lat%M%modes,lat%dim_lat,this%N_rep,size(this%M_n,1))
end subroutine

subroutine get_H(this,lat,Hout)
    !get effective field H 
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    real(8),intent(inout)         ::  Hout(:,:)
#ifdef CPP_FFTW3
    integer ::  i,j,l

    Call this%set_M(lat)
    Call fftw_execute_dft_r2c(this%plan_mag_F, this%M_n, this%M_F)

    this%H_F=cmplx(0.0d0,0.0d0,8)
    do j=1,size(this%M_F,2)
        do i=1,lat%Nmag*3
            do l=1,lat%Nmag*3
                this%H_F(i,j)=this%H_F(i,j)+this%K_F(i,l,j)*this%M_F(l,j)
            enddo
        enddo
    enddo
    Call fftw_execute_dft_c2r(this%plan_H_I, this%H_F, this%H_n)
    Call this%H_internal(this%H_n,Hout,lat%dim_lat,this%N_rep,size(Hout,1))
#else
    ERROR STOP "fft_H%get_H requires CPP_FFTW"
#endif
end subroutine

subroutine get_E_distrib(this,lat,Htmp,E)
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    real(8),intent(inout)         ::  Htmp(:,:)
    real(8),intent(out)           ::  E(:)

    Call this%get_H(lat,Htmp)
    Htmp=Htmp*lat%M%modes_v
    E=sum(reshape(Htmp,[3,lat%Nmag*lat%Ncell]),1)*2.0d0 !not sure about *2.0d0
end subroutine


subroutine get_E(this,lat,Htmp,E)
    class(fft_H),intent(inout)    ::  this
    type(lattice),intent(in)      ::  lat
    real(8),intent(inout)         ::  Htmp(:,:)
    real(8),intent(out)           ::  E

    Call this%get_H(lat,Htmp)
    Htmp=Htmp*lat%M%modes_v
    E=sum(Htmp)
end subroutine


subroutine set_fftw_plans(this)
    class(fft_H),intent(inout)  :: this

    integer(C_INT)  :: N_rep_rev(3)     !reversed N_rep necessary for fftw3 (col-major -> row-major)
    integer(C_int)  :: howmany          !dimension of quantitiy which is fourier-transformed (see FFTW3)
    
#ifdef CPP_FFTW3
    if(.not.allocated(this%M_n)) ERROR STOP "cannot set fftw_plans as M_n not allocated"
    if(.not.allocated(this%M_F)) ERROR STOP "cannot set fftw_plans as M_F not allocated"
    if(.not.allocated(this%H_n)) ERROR STOP "cannot set fftw_plans as H_n not allocated"
    if(.not.allocated(this%H_F)) ERROR STOP "cannot set fftw_plans as H_F not allocated"

    N_rep_rev=this%N_rep(size(this%N_rep):1:-1)
    howmany=int(size(this%M_n,1),C_int)
    this%plan_mag_F= fftw_plan_many_dft_r2c(int(3,C_INT), N_rep_rev, howmany,&
                                           &this%M_n,     N_rep_rev,&
                                           &howmany,      int(1,C_int), &
                                           &this%M_F,     N_rep_rev,&
                                           &howmany,      int(1,C_int), &
                                           &FFTW_FORWARD+FFTW_MEASURE+FFTW_PATIENT)

    this%plan_H_I= fftw_plan_many_dft_c2r(int(3,C_INT), N_rep_rev, howmany,&
                                         &this%H_F,     N_rep_rev,&
                                         &howmany,      int(1,C_int), &
                                         &this%H_n,     N_rep_rev,&
                                         &howmany,      int(1,C_int), &
                                         &FFTW_BACKWARD+FFTW_MEASURE+FFTW_PATIENT)
#else
    ERROR STOP "CANNOT USE FFT_H without FFTW (CPP_FFTW3)"
#endif
end subroutine

end module
