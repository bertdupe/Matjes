module m_tightbinding_r
use m_tb_types
use m_derived_types, only : lattice
use m_occupation, only: calc_occupation
use m_occupation_mult, only: occupation_mult
use m_fermi, only: calc_fermi 
use m_distribution, only: int_distrib,fermi_distrib,dE_fermi_distrib
use m_save_state_r,only: TB_write_states_r, TB_read_states_r
use m_init_Hr
use m_H_tb_public
use m_dos_util
use, intrinsic :: iso_fortran_env, only : output_unit
implicit none
private
public :: tightbinding_r
contains

subroutine tightbinding_r(lat,tb_par)
    type(lattice), intent(in)               :: lat
    type(parameters_TB),intent(in)          :: tb_par

    type(parameters_TB_Hsolve)              :: h_par

    real(8),allocatable         ::  eigval(:)
    complex(8),allocatable      ::  eigvec(:,:)

    real(8)                     :: E_f
    logical                     :: calc_eigval,calc_eigvec

    logical                     :: read_success
    class(H_tb),allocatable     :: H
    integer                     :: i

    procedure(int_distrib),pointer  :: dist_ptr => null()

    h_par=tb_par%H

    !check what to calculate (eigenvalue/eigenvector)
    calc_eigval=tb_par%flow%dos_r.or.tb_par%flow%occ_r.or.tb_par%flow%spec_r.or.tb_par%flow%fermi_r
    calc_eigvec=(tb_par%flow%dos_r.and.(tb_par%is_sc.or.allocated(tb_par%io_dos%bnd))).or.tb_par%flow%occ_r

    !possibly read previous solution
    if(tb_par%flow%read_solution_r)then
        Call TB_read_states_r(eigval,eigvec,read_success)
        calc_eigvec=calc_eigvec.and..not.allocated(eigvec)
        calc_eigval=calc_eigval.and..not.allocated(eigval)
    else
        read_success=.false.
    endif

    !initialize Hamiltonian
    if(calc_eigvec.or.calc_eigval) Call get_Hr(lat,tb_par%io_H,H)

    !calculate eigenvalue/eigenvector
    if(calc_eigvec)then
        write(output_unit,'(//A)') 'Start calculating eigenvectors'
        Call H%get_evec(eigval,eigvec)
        read_success=.false.
    elseif(calc_eigval)then
        write(output_unit,'(//A)') 'Start calculating eigenvalues'
        Call H%get_eval(eigval)
        read_success=.false.
    endif

    !write eigenvalue/eigenvector to file for later reuse
    if(.not.read_success.and.tb_par%flow%write_solution_r) Call TB_write_states_r(eigval,eigvec)

    !write spectrum
    if(tb_par%flow%spec_r)then
        write(output_unit,'(/A)') 'start write spectrum'
        open(newunit=i,file='eigval.dat'); write(i,'(E16.8)') eigval; close(i)
    endif

    !Calculate Fermi energy (only useful without SC)
    E_f=tb_par%io_ef%E_F_in
    if(tb_par%flow%fermi_r)then
        write(output_unit,'(/A)') 'start calculate Fermi energy'
        if(tb_par%is_sc)then
            STOP "calculation of Fermi energy doesn't work when using superconductivity"
        else
            Call calc_fermi(eigval, tb_par%io_EF%N_electrons*h_par%ncell, tb_par%io_ef%kt, E_f)
        endif
    endif

    if(tb_par%flow%dos_r)then
        write(output_unit,'(/A)') 'start calculate DOS'
        if(tb_par%is_sc)then
            Call write_dos_sc(eigval,eigvec,lat,tb_par%io_dos)
        else
            Call write_dos_nc(eigval,eigvec,lat,tb_par%io_dos)
        endif
    endif

    if(tb_par%flow%occ_r)then
        write(output_unit,'(/A)') 'start calculate occupation'
         !maybe use different smearing than EF input
        dist_ptr=>fermi_distrib
        Call calc_occupation(h_par,eigvec,eigval,E_f,tb_par%io_Ef%kt,'occ.dat',dist_ptr)
        dist_ptr=>dE_fermi_distrib
        Call calc_occupation(h_par,eigvec,eigval,E_f,tb_par%io_Ef%kt,'occ_dE.dat',dist_ptr)
    endif

    if(tb_par%flow%occ_mult_r)then
        write(output_unit,'(/A)') 'start calculate multiple occupations'
        Call occupation_mult(h_par,tb_par%io_occ_mult,eigval,eigvec)
    endif

end subroutine 

subroutine write_dos_nc(eval,evec,lat,io_dos) 
    real   (8),intent(in)       :: eval(:)
    complex(8),intent(in)       :: evec(:,:)
    type(lattice),intent(in)    :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    type(dos_nc)                            :: dos

    type(dos_bnd_nc),allocatable            :: dos_bnd(:)
    integer                                 :: idos,Ndos
    character(len=3)                        :: bnd_id
    Call dos%init(io_dos)

    Ndos=0
    if(allocated(io_dos%bnd))then
        Ndos=size(io_dos%bnd,2)
        allocate(dos_bnd(Ndos))
        do idos=1,Ndos
            Call dos_bnd(idos)%init_bnd(io_dos,io_dos%bnd(:,idos))
        enddo
    endif

    Call dos%add(eval)
    do idos=1,Ndos
        Call dos_bnd(idos)%add(eval,evec)
    enddo

    Call dos%print("dos_r_nc.dat")
    do idos=1,Ndos
        write(bnd_id,'(I0.3)') idos
        Call dos_bnd(idos)%print("dos_r_nc_bnd_"//bnd_id//".dat")
    enddo
end subroutine


subroutine write_dos_sc(eval,evec,lat,io_dos) 
    real   (8),intent(in)       :: eval(:)
    complex(8),intent(in)       :: evec(:,:)
    type(lattice),intent(in)    :: lat
    type(parameters_TB_IO_DOS),intent(in)   :: io_dos

    type(dos_sc)                            :: dos

    type(dos_bnd_sc),allocatable            :: dos_bnd(:)
    integer                                 :: idos,Ndos
    character(len=3)                        :: bnd_id
    Call dos%init(io_dos)

    Ndos=0
    if(allocated(io_dos%bnd))then
        Ndos=size(io_dos%bnd,2)
        allocate(dos_bnd(Ndos))
        do idos=1,Ndos
            Call dos_bnd(idos)%init_bnd(io_dos,io_dos%bnd(:,idos))
        enddo
    endif

    Call dos%add(eval,evec)
    do idos=1,Ndos
        Call dos_bnd(idos)%add(eval,evec)
    enddo

    Call dos%print("dos_r_sc.dat")
    do idos=1,Ndos
        write(bnd_id,'(I0.3)') idos
        Call dos_bnd(idos)%print("dos_r_sc_bnd_"//bnd_id//".dat")
    enddo
end subroutine

end module
