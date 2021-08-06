module m_init_heavyside
use m_derived_types
implicit none
private
public :: init_heavyside
contains


subroutine init_heavyside(lat,dim_mode,state)
    type (lattice),intent(in)       :: lat
    integer, intent(in)             :: dim_mode
    real(8),pointer,intent(inout)   :: state(:)
    ! internal variable
    integer         :: nmag
    real(8),pointer :: state_x(:,:,:,:)
   
    nmag=lat%nmag
    state_x(1:3,1:nmag,1:lat%dim_lat(1),1:lat%dim_lat(2)*lat%dim_lat(3))=>state

    state_x=0.0d0
    state_x(3,:,:,:)=-1.0d0
    state_x(3,:,1:lat%dim_lat(1)/2,:)=1.0d0
end subroutine 
end module m_init_heavyside
