module m_input_H_types
implicit none
public


type :: io_H_base
    logical :: is_set=.false.
end type

type,extends(io_H_base) :: io_H_aniso
    real(8),allocatable :: val(:)
end type

type,extends(io_H_base) :: io_H_zeeman
    real(8)     :: c_zeeman=-1.0d0 !constant factor to furthermore rescale zeeman energy 
end type

type,extends(io_H_base) :: io_H_J
    real(8), allocatable   ::  val(:)   !REMEMBER TO CHANGE GET_SHELL_NUMBER IF MORE MAGNETIC ATOMS ARE ADDED
end type

type,extends(io_H_base) :: io_H_D
    real(8), allocatable   ::  val(:)   !REMEMBER TO CHANGE GET_SHELL_NUMBER IF MORE MAGNETIC ATOMS ARE ADDED
end type

type,extends(io_H_base) :: io_H_ME_J
    real(8), allocatable   ::  val(:)
end type

type,extends(io_H_base) :: io_H_ME_D
    real(8), allocatable   ::  val(:)
end type

type,extends(io_H_base) :: io_H_F
    real(8), allocatable   ::  val(:) !Force constant for next nearest neighbor shell
end type

!type,extends(io_H_base) :: io_H_ME
!    real(8), allocatable   ::  sym(:),asym(:)   !REMEMBER TO CHANGE GET_SHELL_NUMBER IF MORE MAGNETIC ATOMS ARE ADDED
!end type

type :: io_H
    type(io_H_aniso)    :: aniso
    type(io_H_zeeman)   :: zeeman
    type(io_H_ME_J)     :: ME_J
    type(io_H_ME_D)     :: ME_D
    type(io_H_J)        :: J
    type(io_H_D)        :: D
    type(io_H_F)        :: F
    contains
    procedure :: get_shell_number
end type

contains

pure function get_shell_number(this)result(N_shell)
    integer                 :: N_shell
    class(io_H),intent(in)  :: this

    N_shell=0
    if(this%J%is_set) N_shell=max(N_shell,size(this%J%val))
    if(this%D%is_set) N_shell=max(N_shell,size(this%D%val))
    if(this%ME_J%is_set) N_shell=max(N_shell,size(this%ME_J%val))
    if(this%ME_D%is_set) N_shell=max(N_shell,size(this%ME_D%val))
    if(this%F%is_set) N_shell=max(N_shell,size(this%F%val))
end function

end module
