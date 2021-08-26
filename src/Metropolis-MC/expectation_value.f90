module m_mc_exp_val
use m_constants, only : pi
use m_mc_track_val, only: track_val
use m_fluct,only: fluct_parameters, eval_fluct
use mpi_basic
implicit none
private
public :: exp_val, measure_add!,measure_eval,print_av
!public MPI_routines
public :: measure_scatterv, measure_gatherv, measure_bcast, measure_reduce
public :: measure_save, measure_read, print_spatial_fluct
character(*),parameter  :: save_name='expval_save.dat'
integer,protected       :: MPI_custom_type          ! mpi variable for direct mpi operations with this type 
logical,protected       :: MPI_custom_set=.false.   !check is MPI_custom_type is set

!THESE ENTRIES HAVE TO BE UPDATED EVERYTIME EXP_VAL IS MODIFIED, OTHERWISE MPI-STUFF WILL FAIL
integer,parameter       :: blocks(*)=[1, 1,1, 3,3, 1,1,1,1,1,3, 1,1,1,1, 1] !size of each element of exp_val
integer,parameter       :: bnd_real(2) =[ 1, 11] !initial and final entry of real numbers
integer,parameter       :: bnd_cmplx(2)=[12, 15] !initial and final entry of complex numbers
integer,parameter       :: bnd_int(2)  =[16, 16] !initial and final entry of integer numbers

type exp_val
    sequence !should be sequence for MPI-type
    !! contribution of the different energy parts
    real(8) :: kt=0.0d0 !temperature should always be at first position so it does not get added measure_reduce
    ! energy
    real(8) :: E=0.0d0
    real(8) :: E_sq=0.0d0
    ! magnetization
    real(8) :: M(3)=0.0d0
    real(8) :: M_sq(3)=0.0d0
    ! topo-stuff
    real(8) :: Q_sq=0.0d0
    real(8) :: Qp=0.0d0
    real(8) :: Qm=0.0d0
    real(8) :: Qp_sq=0.0d0
    real(8) :: Qm_sq=0.0d0
    real(8) :: vortex(3)=0.0d0
    !fluctuation stuff
    complex(8)  :: MipMip =cmplx(0.0d0,0.0d0,8) !<Mi+Mi+>
    complex(8)  :: MipMim =cmplx(0.0d0,0.0d0,8) !<Mi+Mi-> (is real)
    complex(8)  :: MipMjp =cmplx(0.0d0,0.0d0,8) !<Mi+Mj+>

    !<Mj+Mi->
    complex(8),allocatable:: MjpMim_ij(:,:) 

    integer :: N_add=0 !counts how often values have been added
end type

contains 

subroutine measure_save(this)
    use, intrinsic :: iso_fortran_env, only : output_unit
    type(exp_val),intent(in)    :: this(:)
    integer :: io_unit
    integer :: i

    write(output_unit,'(/2A/)') "Saving expectation values in file: ",save_name 
    open(newunit=io_unit,file=save_name)
    write(io_unit,*) size(this)
    ERROR STOP "MEASURE WRITE NO LONGER WORKS SINCE FLUCTUATIONS HAVE BEEN ADDED"
    !do i=1,size(this)
    !    write(io_unit,*) this(i)
    !end do
    close(io_unit)
end subroutine

subroutine measure_read(this)
    use, intrinsic :: iso_fortran_env, only : output_unit,error_unit
    type(exp_val),intent(out)   :: this(:)
    integer :: io_unit
    integer :: i,NT
    logical :: fexist

    inquire(file=save_name,exist=fexist)
    if(.not.fexist)then
        write(error_unit,'(3A)') "Trying to read expectation value from file ",save_name," but it does not exist"
        write(error_unit,'(3A)') "Remove 'expval_save' from input or find file?"
        STOP  
    endif
    write(output_unit,'(/2A/)') "Reading expectation values from file: ",save_name 
    open(newunit=io_unit,file=save_name)
    read(io_unit,*) NT
    if(NT/=size(this)) ERROR STOP "Trying to read expectation values, but the number of temperatures has changes"
    ERROR STOP "MEASURE READ NO LONGER WORKS SINCE FLUCTUATIONS HAVE BEEN ADDED"
    !do i=1,size(this)
    !    read(io_unit,*) this(i)
    !end do
    close(io_unit)
end subroutine

#ifdef CPP_MPI    
subroutine set_custom_type()
    if(.not.MPI_custom_set)then
        Call set_MPI_type(blocks,bnd_real,bnd_cmplx,bnd_int,MPI_custom_type)
        MPI_custom_set=.true.
    endif
end subroutine
#endif

subroutine measure_gatherv(this,com)
    type(exp_val),intent(inout)     :: this(:)
    class(mpi_distv),intent(in)     :: com
#ifdef CPP_MPI    
    integer                         :: ierr

    Call set_custom_type()
    Call MPI_Gatherv(this(1),com%cnt(com%id+1),MPI_custom_type,this(1),com%cnt,com%displ,MPI_custom_type,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine measure_scatterv(this,com)
    type(exp_val),intent(inout)     :: this(:)
    class(mpi_distv),intent(in)     :: com
#ifdef CPP_MPI    
    integer                         :: ierr

    Call set_custom_type()

    Call MPI_Scatterv(this(1),com%cnt,com%displ,MPI_custom_type,this(1),com%cnt(com%id+1),MPI_custom_type,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

subroutine measure_reduce(this,com)
    !sums all values except for the temperature over the MPI-communicator(kt)
    type(exp_val),intent(inout)     :: this
    class(mpi_type),intent(in)      :: com
#ifdef CPP_MPI    
    integer         :: ierr
    integer         :: nsize(3)

    if(com%NP/=1)then   !only necessary is there is something to reduce
        nsize(1)=sum(blocks(bnd_real (1)+1:bnd_real (2)))    !+1 to ignore kt
        nsize(2)=sum(blocks(bnd_cmplx(1)  :bnd_cmplx(2)))
        nsize(3)=sum(blocks(bnd_int  (1)  :bnd_int  (2)))
        if(com%ismas)then
            Call MPI_Reduce(MPI_IN_PLACE,this%E     ,nsize(1),MPI_DOUBLE_PRECISION,MPI_SUM,com%mas,com%com,ierr)
            Call MPI_Reduce(MPI_IN_PLACE,this%MipMim,nsize(2),MPI_DOUBLE_COMPLEX  ,MPI_SUM,com%mas,com%com,ierr)
            Call MPI_Reduce(MPI_IN_PLACE,this%N_add ,nsize(3),MPI_Int             ,MPI_SUM,com%mas,com%com,ierr)
        else
            Call MPI_Reduce(this%E      ,this%E     ,nsize(1),MPI_DOUBLE_PRECISION,MPI_SUM,com%mas,com%com,ierr)
            Call MPI_Reduce(this%MipMim ,this%MipMip,nsize(2),MPI_DOUBLE_COMPLEX  ,MPI_SUM,com%mas,com%com,ierr)
            Call MPI_Reduce(this%N_add  ,this%N_add ,nsize(3),MPI_Int             ,MPI_SUM,com%mas,com%com,ierr)
        endif
    endif
#else
    continue
#endif
end subroutine

subroutine measure_bcast(this,com)
    type(exp_val),intent(inout)     :: this
    class(mpi_type),intent(in)      :: com
#ifdef CPP_MPI    
    integer                         :: ierr

    Call set_custom_type()
    Call MPI_Bcast(this,1,MPI_custom_type,com%mas,com%com,ierr)
#else
    continue
#endif
end subroutine

!!!!! print spatial distribution for <Mj+Mi->
subroutine print_spatial_fluct(this,com,io_unit_in)
    use m_constants, only : k_b
    use m_convert
    use mpi_basic, only: mpi_type
    use m_io_files_utils, only: open_file_write,close_file
    type(exp_val),intent(inout)     :: this
    class(mpi_type),intent(in)      :: com
    integer,optional                :: io_unit_in
  ! internal
    integer             :: shape_MiMj(2),io_file(2),iomp
    character(len=50)   :: file_name(2),form
    real(8)             :: av_Nadd

    av_Nadd=1.0d0/real(this%N_add,8)
    shape_MiMj=shape(this%MjpMim_ij)
    form=convert('(',shape_MiMj(1),'(E20.12E3,2x))')

    file_name(1)=convert('fluct_Re_MjpMim_per_site_',this%kT/k_B,'.dat')
    io_file(1)=open_file_write(file_name(1))
    write(io_file(1),form) real(this%MjpMim_ij)*av_Nadd
    call close_file(file_name(1),io_file(1))

    file_name(2)=convert('fluct_Im_MjpMim_per_site_',this%kT/k_B,'.dat')
    io_file(2)=open_file_write(file_name(2))
    write(io_file(2),form) aimag(this%MjpMim_ij(:,iomp))*av_Nadd
    call close_file(file_name(2),io_file(2))
end subroutine


subroutine measure_add(this,lat,state_prop,Q_neigh,fluct_val)
    use m_topo_commons
    use m_derived_types,only: lattice
    type(exp_val),intent(inout)             :: this
    type(lattice),intent(in)                :: lat
    type(track_val),intent(in)              :: state_prop
    integer,intent(in),allocatable          :: Q_neigh(:,:)
    type(fluct_parameters),intent(in)       :: fluct_val

    !put that into state_prop as well?
    real(8)     :: dumy(5),qeulerp,qeulerm


    this%N_add=this%N_add+1

    this%E=this%E+state_prop%E_total
    this%E_sq=this%E_sq+state_prop%E_total**2

    this%M=this%M+state_prop%Magnetization
    this%M_sq=this%M_sq+state_prop%Magnetization**2

    ! calculate the topocharge
    if(allocated(Q_neigh))then
        dumy=get_charge(lat,Q_neigh)
        qeulerp=dumy(1)
        qeulerm=-dumy(2)

        this%Qp=this%Qp+qeulerp
        this%Qm=this%Qm+qeulerm
        this%Q_sq=this%Q_sq+(qeulerp-qeulerm)**2
        this%Qp_sq=this%Qp+qeulerp**2
        this%Qm_sq=this%Qm+qeulerm**2
        this%vortex=this%vortex+dumy(3:5)
    endif

    if(fluct_val%l_use) Call eval_fluct(this%MjpMim_ij, &
                                       &this%MipMip, &
                                       &this%MipMim, &
                                       &this%MipMjp, &
                                       &lat,fluct_val)
end subroutine
end module
