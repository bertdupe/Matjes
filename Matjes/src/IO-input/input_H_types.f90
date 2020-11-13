module m_input_H_types
implicit none
public


type :: io_H_base
    logical :: is_set=.false.
end type

type,extends(io_H_base) :: io_H_aniso
    real(8),allocatable :: val(:)
end type

type :: io_H
    type(io_H_aniso)    :: aniso
end type

end module
