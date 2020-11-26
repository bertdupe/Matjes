#if 0
module m_modes_variables
use m_basic_types

! related to the number of atoms in the shell
type sites
     type (vec_point), allocatable :: site(:)
     integer, allocatable :: num(:)
end type sites

type shell
     type (sites), allocatable :: link(:)
     integer, allocatable :: num(:)
end type shell

type order
     type (shell), allocatable :: shell_num(:)
     integer, allocatable :: num(:)
end type order

! old variable
type point_shell_mode
     type(vec_point), allocatable, dimension(:) :: shell
end type point_shell_mode

end module m_modes_variables
#endif
