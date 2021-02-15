module m_orders_initialize
implicit none
public 
contains
subroutine orders_initialize(lat,extpar_io)
    !subroutine which sets the initial values for all order parameters of the lattice, except for the very basic initialization in lat%init_order
    !checks if "*_init.dat" files are provided, if not tries to read from init.config
    !afterwards, tries to punch lattice
    use m_type_lattice,only: number_different_order_parameters,order_parameter_name
    use m_type_lattice,only: lattice
    use m_init_config, only: init_config_lattice
    use m_punch_init, only: punch_lattice
    use m_input_types, only: extpar_input
    type(lattice), intent(inout)    :: lat
    type(extpar_input),intent(in)   :: extpar_io
    !internal
    logical         :: initialized(number_different_order_parameters)
    logical         :: fexist

    !legacy checks for discontinues input files
    inquire(file='SpinSTMi.dat',exist=fexist)
    if (fexist) stop 'please adjust the file "SpinSTMi.dat" to "*_init.dat", eg."magnetic_init.dat"'
    inquire(file='init-modes.dat',exist=fexist)
    if (fexist) stop 'please rename the file "init-modes.dat" to "*_init.dat", eg."magnetic_init.dat"'

    initialized=.false.
    Call lat%read_order(isinit_opt=initialized)
    Call lat%read_order("_init.dat",isinit_opt=initialized) !backwards compatible name, remove at some point?
    Call init_config_lattice(lat,initialized,extpar_io)
    Call punch_lattice(lat)
end subroutine
end module
