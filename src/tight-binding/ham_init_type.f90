module m_ham_init_type
!parameter which contains all essential parameters for shape of Hamiltonian
use m_TB_types, only: parameters_TB_IO_H
private
public parameters_ham_init
type parameters_ham_init
    !parameter with all initialization types for the Hamiltonian creation
    integer  ::  nspin
    integer  ::  ncell
    integer  ::  norb
    integer  ::  nsc

    real(8)  ::  Ebnd(2)
    integer  ::  estNe
    real(8)  ::  diag_acc
contains
    procedure       :: init => init_ham_init
end type


contains

subroutine init_ham_init(this,io)
    class(parameters_ham_init),intent(out)  ::  this
    class(parameters_TB_IO_H),intent(in)    ::  io
    this%nspin   = io%nspin    
    this%ncell   = io%ncell   
    this%norb    = io%norb     
    this%nsc     = io%nsc     

    this%Ebnd    = io%Ebnd    
    this%estNe   = io%estNe    
    this%diag_acc= io%diag_acc
end subroutine
end module
