module m_orders_initialize
implicit none
public 

contains
subroutine orders_initialize(lat)
    use m_type_lattice,only: number_different_order_parameters,order_parameter_name
    use m_type_lattice,only: lattice
    use m_init_config, only: init_config_new
    use m_punch_init, only: punch_lattice
    type(lattice), intent(inout) :: lat

    logical         ::  initialized(number_different_order_parameters)

    Call lat%read_order(fexist_out=initialized)
    Call init_config_new(lat,initialized)
    Call punch_lattice(lat)
end subroutine



end module
