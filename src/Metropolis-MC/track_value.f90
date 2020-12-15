module m_mc_track_val
use m_derived_types, only : lattice,number_different_order_parameters
use m_H_public, only: t_H,energy_all
implicit none
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
    contains
    procedure :: init
end type

contains 
subroutine init(this,lat,hams,cone_in)
    class(track_val),intent(out)    :: this
    type(lattice),intent(in)        :: lat
    class(t_H), intent(in)          :: Hams(:)
    real(8),intent(in)              :: cone_in
    integer     ::  i
    logical     ::  used(number_different_order_parameters)

    this%magnetization=sum(lat%M%modes_3,2) !sum over magnetization of all magnetic atoms without magnetic moment, probably not what is really wanted
    this%E_total=energy_all(Hams,lat)
    this%cone=cone_in
    this%dim_bnd=0
    Call lat%used_order(used)
    do i=1,number_different_order_parameters
        if(.not.used(i)) cycle 
        this%dim_bnd(1,i)=1
        this%dim_bnd(2,i)=lat%get_order_dim(i)
    enddo
END subroutine 
end module
