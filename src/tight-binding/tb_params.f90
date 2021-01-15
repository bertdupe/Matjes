module m_tb_params
use m_tb_types
use m_ham_arrange
implicit none
public

type(parameters_TB),public,protected :: TB_params

contains
subroutine set_TB_params(lat)
    use m_rw_TB, only: get_parameters_io_TB
    use m_derived_types, only: lattice
    type(lattice), intent(in) :: lat
    integer ::  i

    Call get_parameters_io_TB(TB_params)

    !use superconducting version if any delta is set
    TB_params%is_sc=allocated(TB_params%io_H%del)
    !use magnetic version if any atom has both a magnetic moment and an orbital or if superconducting is set
    TB_params%is_mag=any(lat%cell%atomic%moment/=0.0d0.and.lat%cell%atomic%orbitals>0).or.TB_params%is_sc
    !set dimension of Hamiltonian parameters
    associate( par=>TB_params%io_H)
        if(TB_params%is_sc) par%nsc=2
        if(TB_params%is_mag) par%nspin=2
        par%ncell=lat%Ncell
        par%norb_at=lat%cell%atomic%orbitals
        par%norb_at_off=[(sum(par%norb_at(1:i-1)),i=1,size(par%norb_at))]
        par%norb=sum(par%norb_at)
        par%dimH=par%nsc*par%nspin*par%norb*par%ncell
        do i=1,size(par%hop_io)
            Call par%hop_io(i)%check(lat)
        enddo
        Call par%hop%set(par%hop_io,lat,par%nspin)
    end associate

    !REMOVE-> move to H_io
    !set dimensions of the Hamiltonian
    if(TB_params%is_sc) TB_params%H%nsc=2
    if(TB_params%is_mag) TB_params%H%nspin=2
    TB_params%H%ncell=lat%Ncell
    TB_params%H%norb=sum(lat%cell%atomic%orbitals)
    Call TB_params%H%upd()

    TB_params%H%sparse=TB_params%io_H%sparse
    TB_params%H%i_diag=TB_params%io_H%i_diag
    TB_params%H%estNe=TB_params%io_H%estNe
    TB_params%H%Ebnd=TB_params%io_H%Ebnd
    TB_params%H%diag_acc=TB_params%io_H%diag_acc

    !TB_params%H%rearrange=TB_params%io_H%rearrange !this is not really implementend and most probably breakes dos or something else
end subroutine

end module
