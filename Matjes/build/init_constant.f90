module m_init_constant
use m_derived_types

private
public :: init_constant_config


contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! value that initialized a constant vector field (typically for the B_field, E_field...)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine init_constant_config(my_lattice,mode_name,start,end,ext_param)
implicit none
type (lattice), intent(inout) :: my_lattice
integer, intent(in) :: start,end
character(len=*), intent(in) :: mode_name
type (simulation_parameters), intent(in) :: ext_param
!internal
integer :: shape_lattice(4),i,j,k,l
real(kind=8), allocatable :: value(:)

allocate(value(end-start+1))
call get_value_field(value,mode_name,ext_param)

shape_lattice=shape(my_lattice%l_modes)

do l=1,shape_lattice(4)
   do k=1,shape_lattice(3)
      do j=1,shape_lattice(2)
         do i=1,shape_lattice(1)
           my_lattice%l_modes(i,j,k,l)%w(start:end)=value
         enddo
      enddo
   enddo
enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! subroutine that finds the field depending on the mode name
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_value_field(value,mode_name,ext_param)
implicit none
character(len=*), intent(in) :: mode_name
type (simulation_parameters), intent(in) :: ext_param
real(kind=8), intent(inout) :: value(:)
! internal

select case (adjustl(mode_name))
  case('Efield')
    value=ext_param%E_ext%value
  case('Bfield')
    value=ext_param%H_ext%value
  case('temperature')
    value(1)=ext_param%ktini%value
  case default
    stop 'field name could not be found in get_value_field'
end select

end subroutine

end module m_init_constant
