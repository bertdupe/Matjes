module m_FT_Ham_coo
use m_FT_Ham_base
use m_H_type
implicit none

private
public :: H_inp_k_coo

type,abstract,extends(t_H_base) :: H_inp_k_coo
    private
    integer             :: nnz=0  !number of entries in sparse matrix
    complex(8),allocatable :: val(:)
    integer,allocatable :: row(:),col(:)

    contains
!    procedure :: create
!    procedure :: destroy
    procedure :: init

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

end module m_FT_Ham_coo
