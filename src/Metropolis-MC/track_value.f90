module m_mc_track_val
use m_derived_types, only : lattice,number_different_order_parameters
use m_H_public, only: energy_all
use m_H_type, only:t_H
implicit none
private
public track_val
type track_val
    !values to track the changes during MC-step
    real(8) :: magnetization(3)=0.0d0
    real(8) :: E_total=0.0d0
    ! statistics on the MC
    real(8) :: acc=0.0d0
    real(8) :: rate=0.0d0
    real(8) :: nb=0.0d0
    real(8) :: cone=0.0d0
    !parameters for single evaluation
    integer ::  dim_bnd(2,number_different_order_parameters)
    procedure(int_sample),pointer   :: spin_sample => null()
    contains
    procedure :: init
end type


abstract interface
    function int_sample(this,i_spin,lat,hams)result(spin_new)
        import track_val,lattice,t_H
        class(track_val),intent(inout)  :: this
        integer,intent(in)              :: i_spin
        type(lattice),intent(in)        :: lat
        class(t_H), intent(in)          :: Hams(:)
        real(8)                         :: spin_new(3)
    end function
end interface


contains 
subroutine init(this,lat,hams,io_MC)
    use m_MC_io,only: MC_input
    class(track_val),intent(out)    :: this
    type(lattice),intent(in)        :: lat
    class(t_H), intent(in)          :: Hams(:)
    class(MC_input),intent(in)      :: io_MC
    integer     ::  i
    logical     ::  used(number_different_order_parameters)

    this%magnetization=sum(lat%M%modes_3,2) !sum over magnetization of all magnetic atoms without magnetic moment, probably not what is really wanted
    this%E_total=energy_all(Hams,lat)
    this%cone=io_MC%cone
    this%dim_bnd=0
    Call lat%used_order(used)
    do i=1,number_different_order_parameters
        if(.not.used(i)) cycle 
        this%dim_bnd(1,i)=1
        this%dim_bnd(2,i)=lat%get_order_dim(i)
    enddo

    if(io_MC%ising)then
        this%spin_sample=>sample_ising
    elseif(io_MC%underrelax)then
        this%spin_sample=>sample_underrelax
    elseif(io_MC%overrelax)then
        this%spin_sample=>sample_overrelax
    elseif(io_MC%equi)then
        this%spin_sample=>sample_equi
    elseif(io_MC%sphere)then
        this%spin_sample=>sample_sphere
    endif
END subroutine 

!wrappers for sampling
function sample_ising(this,i_spin,lat,hams)result(spin_new)
    class(track_val),intent(inout)  :: this
    integer,intent(in)              :: i_spin
    type(lattice),intent(in)        :: lat
    class(t_H), intent(in)          :: Hams(:)
    real(8)                         :: spin_new(3)

    spin_new=-lat%M%modes_3(:,i_spin)
end function


function sample_underrelax(this,i_spin,lat,hams)result(spin_new)
    use m_relaxtyp, only: underrelax 
    class(track_val),intent(inout)  :: this
    integer,intent(in)              :: i_spin
    type(lattice),intent(in)        :: lat
    class(t_H), intent(in)          :: Hams(:)
    real(8)                         :: spin_new(3)

    spin_new=underrelax(i_spin,lat,Hams)
end function


function sample_overrelax(this,i_spin,lat,hams)result(spin_new)
    use m_relaxtyp, only: overrelax 
    class(track_val),intent(inout)  :: this
    integer,intent(in)              :: i_spin
    type(lattice),intent(in)        :: lat
    class(t_H), intent(in)          :: Hams(:)
    real(8)                         :: spin_new(3)

    spin_new=overrelax(i_spin,lat,Hams)
end function


function sample_equi(this,i_spin,lat,hams)result(spin_new)
    use m_sampling, only: equirep
    class(track_val),intent(inout)  :: this
    integer,intent(in)              :: i_spin
    type(lattice),intent(in)        :: lat
    class(t_H), intent(in)          :: Hams(:)
    real(8)                         :: spin_new(3)

    Call cone_update(this%cone,this%rate)
    spin_new=equirep(lat%M%modes_3(:,i_spin),this%cone)
end function


function sample_sphere(this,i_spin,lat,hams)result(spin_new)
    use m_sampling, only: sphereft
    class(track_val),intent(inout)  :: this
    integer,intent(in)              :: i_spin
    type(lattice),intent(in)        :: lat
    class(t_H), intent(in)          :: Hams(:)
    real(8)                         :: spin_new(3)

    Call cone_update(this%cone,this%rate)
    spin_new=sphereft(lat%M%modes_3(:,i_spin),this%cone)
end function

pure subroutine  cone_update(cone,rate)
    use m_constants, only : pi
    real(8),intent(inout)   :: cone
    real(8),intent(in)      :: rate

    ! cone angle update and maximal magnetic moment change
    if ((rate.gt.0.5d0).and.(cone.lt.pi))then
         cone=cone+0.0001d0
    elseif ((rate.lt.0.50d0).and.(cone.gt.0.01d0)) then
         cone=cone-0.0001d0
    endif
end subroutine


end module
