module m_TB_types
!module which contains the main tight-binding data type parameters_TB

!io sub-types
use m_parameters_rw_high,           only: parameters_IO_highs
use m_parameters_TB_IO_DOS,         only: parameters_TB_IO_DOS
use m_parameters_TB_IO_flow,        only: parameters_TB_IO_flow
use m_parameters_TB_IO_ef,          only: parameters_TB_IO_ef
use m_parameters_TB_IO_OCC_mult,    only: parameters_TB_IO_OCC_mult
use m_parameters_TB_IO_kint,        only: parameters_TB_IO_kint
use m_parameters_TB_IO_H,           only: parameters_TB_IO_H

implicit none
public
private init_parameters_TB, bcast_parameters_TB

type parameters_TB
    !main tight-binding data type which contains all the io-input and some derived quantities
    type(parameters_TB_IO_H)            :: io_H
    type(parameters_TB_IO_EF)           :: io_ef
    type(parameters_TB_IO_DOS)          :: io_dos
    type(parameters_IO_HIGHS)           :: io_highs
    type(parameters_TB_IO_OCC_MULT)     :: io_occ_mult
    type(parameters_TB_IO_kint)         :: io_kmesh
    type(parameters_TB_IO_flow)         :: flow
    logical         ::  is_mag=.False. !Hamiltonian has spins
    logical         ::  is_sc=.False. !Hamiltonian is superconducting-> everything doubles to include creators and destructors
contains
    procedure :: read_file => parameters_TB_read 
    procedure :: init   => init_parameters_TB
    procedure :: bcast  => bcast_parameters_TB
end type

contains

subroutine parameters_TB_read(TB_params,fname)
    use m_io_utils
    use m_io_files_utils
    use m_io_read_util
    class(parameters_TB),intent(inout)  :: TB_params
    character(len=*), intent(in)        :: fname
    integer     ::  io
    
    io=open_file_read(fname)
    Call TB_params%io_H%read_file    (io,fname)
    Call TB_params%io_ef%read_file   (io,fname)
    Call TB_params%io_dos%read_file  (io,fname)
    Call TB_params%flow%read_file   (io,fname)
    Call TB_params%io_highs%read_file('k',io,fname)
    Call TB_params%io_occ_mult%read_file(io,fname)
    Call TB_params%io_kmesh%read_file(io,fname)
    call close_file(fname,io)
end subroutine

subroutine bcast_parameters_TB(this,comm)
    use mpi_basic
    use mpi_util
    class(parameters_TB),intent(inout)  :: this
    type(mpi_type),intent(in)           :: comm
#ifdef CPP_MPI
    Call bcast(this%is_mag,comm)
    Call bcast(this%is_sc, comm)

    Call this%io_H%bcast        (comm)
    Call this%io_ef%bcast       (comm)
    Call this%io_dos%bcast      (comm)
    Call this%io_highs%bcast    (comm)
    Call this%io_occ_mult%bcast (comm)
    Call this%io_kmesh%bcast    (comm)
    Call this%flow%bcast        (comm)
#else
    continue
#endif
end subroutine

subroutine init_parameters_TB(TB_params,lat)
    !subroutine which does operation to prepare the TB-hamiltonian 
    ! starting from the pure input using lattice properties etc.
    use m_derived_types, only: lattice
    use m_dos_io, only: dos_get_ind, dos_get_orb
    class(parameters_TB),intent(inout)  :: TB_params
    type(lattice), intent(in)           :: lat

    integer :: i
    logical :: fexist

    !use superconducting version if any delta is set
    TB_params%is_sc=allocated(TB_params%io_H%del_io)
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
        inquire(file='delta_onsite.inp',exist=fexist)
        if(fexist)then
            Call par%del%read_file('delta_onsite.inp',lat)
        else
            if(allocated(par%del_io)) Call par%del%set(par%del_io,lat,par%norb_at_off)
        endif
        par%use_scf=allocated(par%del_scf_io)
        if(par%use_scf) Call par%del%set_scf(par%del_scf_io,lat,par%norb_at_off)

        if(.not.par%wann_io%is_set)then
            Call par%wann_io%combine_updn(par%wann_io_up,par%wann_io_dn)
        endif
        if(par%wann_io%is_set.and..true.) Call par%wann_io%rearrange_spin() !add additional parameter to control spin-rearangement
    end associate
    if(allocated(TB_params%io_dos%bnd_io))then
        do i=1,size(TB_params%io_dos%bnd_io)
            Call TB_params%io_dos%bnd_io(i)%check(lat)
        enddo
        Call dos_get_ind(TB_params%io_dos%bnd_io,lat,TB_params%io_H%nspin,TB_params%io_H%norb_at,TB_params%io_H%norb_at_off,TB_params%io_dos%bnd)
    endif
    if(allocated(TB_params%io_dos%orb_io))then
        do i=1,size(TB_params%io_dos%orb_io)
            Call TB_params%io_dos%orb_io(i)%check(lat)
        enddo
        Call dos_get_orb(TB_params%io_dos%orb_io,lat,TB_params%io_H%nspin,TB_params%io_H%norb_at,TB_params%io_H%norb_at_off,TB_params%io_dos%orb)
    endif
end subroutine

end module
