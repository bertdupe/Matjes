module m_FT_Ham_coo
use m_FT_Ham_base
use m_H_type
use m_parameters_FT_Ham
implicit none

private :: mode
public :: H_inp_k_coo

#if defined CPP_CUDA
integer  ::        mode=6
#elif defined CPP_MKL
integer  ::        mode=1
#elif defined CPP_EIGEN
integer  ::        mode=2
#elif defined CPP_BLAS
integer  ::        mode=4
#else
integer  ::        mode=5
#endif

type,abstract,extends(t_H_base) :: H_inp_k_coo
    private
    type(parameters_FT_HAM_IO)  :: ham_io
    integer                     :: nnz=0  !number of entries in sparse matrix
    complex(8),allocatable      :: val(:)
    integer,allocatable         :: row(:),col(:)

    contains
!    procedure :: create
    procedure   :: destroy_coo
    procedure   :: init
    procedure   :: pop_par

end type

contains

subroutine init(this,val,row,col)
    !constructs FT Hamiltonian manually for a given sparse coo-input
    class(H_inp_k_coo),intent(inout)       :: this
    complex(8),intent(in)               :: val(:)  !all entries between 2 cell sites of considered orderparameter
    integer,intent(in)                  :: row(size(val)),col(size(val))


    if(this%is_set()) STOP "cannot set FT hamiltonian as it is already set"
    this%nnz=size(val)
    this%val=val
    this%row=row
    this%col=col

end subroutine

subroutine pop_par(this,nnz,val,row,col)
    class(H_inp_k_coo),intent(inout)     :: this
    integer,intent(out)                 :: nnz
    complex(8),allocatable,intent(out)  :: val(:)
    integer,allocatable,intent(out)     :: row(:),col(:)

    nnz=this%nnz
    this%dimH=0
    this%nnz=0
    Call MOVE_ALLOC(this%val,val)
    Call MOVE_ALLOC(this%row,row)
    Call MOVE_ALLOC(this%col,col)
    Call this%destroy
end subroutine

subroutine destroy_coo(this)
    class(H_inp_k_coo),intent(inout)    :: this

    if(this%is_set())then
        this%nnz=0
        if(allocated(this%val   )) deallocate(this%val   )
        if(allocated(this%col   )) deallocate(this%col   )
        if(allocated(this%row   )) deallocate(this%row   )
    endif
end subroutine

end module m_FT_Ham_coo
