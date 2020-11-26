!
!subroutine that carries out molecular dynamics

module m_molecular_dynamics
implicit none
contains

subroutine molecular_dynamics(my_lattice,io_simu,Hams)
   use m_derived_types, only : lattice,io_parameter,number_different_order_parameters
   use m_energy_output_contribution, only:Eout_contrib_init, Eout_contrib_write
   use m_constants, only : k_b
   use m_energyfield, only : write_energy_field
   use m_solver_order
   use m_H_public

   !!!!!!!!!!!!!!!!!!!!!!!
   ! arguments
   !!!!!!!!!!!!!!!!!!!!!!!

   type(lattice), intent(inout) :: my_lattice
   type(io_parameter), intent(in) :: io_simu
   class(t_H), intent(in) :: Hams(:)

   !!!!!!!!!!!!!!!!!!!!!!!
   ! internal variables
   !!!!!!!!!!!!!!!!!!!!!!!
   real(8) :: Edy
   logical :: used(number_different_order_parameters)

   Call my_lattice%used_order(used)

   Edy=energy_all(Hams,my_lattice)


end subroutine

end module
