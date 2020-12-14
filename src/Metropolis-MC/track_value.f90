module m_mc_track_val
use m_derived_types, only : lattice
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
    contains
    procedure :: init
end type

contains 
subroutine init(this,lat,hams,cone_in)
    class(track_val),intent(out)    :: this
    type(lattice),intent(in)        :: lat
    class(t_H), intent(in)          :: Hams(:)
    real(8),intent(in)              :: cone_in

    this%magnetization=sum(lat%M%modes_3,2) !sum over magnetization of all magnetic atoms without magnetic moment, probably not what is really wanted
    this%E_total=energy_all(Hams,lat)
    this%cone=cone_in
END subroutine 
end module
