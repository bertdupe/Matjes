module m_H_type
!module containing the basic polymorphic Hamiltonian class t_H
use m_derived_types, only : lattice
use m_derived_types, only: operator_real_order_N
implicit none

type,abstract :: t_H
contains
procedure(int_eval_single),deferred :: eval_single
procedure(int_eval_all),deferred    :: eval_all
procedure(int_set_H),deferred       :: set_H
end type
private
public t_H

interface
subroutine int_set_H(this,energy_in,lat)
    import t_H,lattice, operator_real_order_N
	!need to get rid of dim_mode input
	class(t_H),intent(inout)  :: this
	type(operator_real_order_N)     :: energy_in
	type(lattice),intent(in)    	:: lat
end subroutine
subroutine int_eval_all(this,E, lat)
    import t_H,lattice
    class(t_H),intent(in)    ::  this
    type(lattice),intent(in)    ::  lat
    real(8),intent(out)			::	E
end subroutine
subroutine int_eval_single(this,E,i_m, lat)
    import t_H,lattice
    class(t_H),intent(in)    ::  this
    type(lattice),intent(in)    ::  lat
	integer,intent(in)			::	i_m
    real(8),intent(out)			::	E
end subroutine
end interface

end module
