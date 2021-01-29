module m_rw_TB
use m_TB_types
use m_io_utils
use m_io_files_utils
use m_io_read_util
implicit none
private
public :: rw_TB
contains

subroutine rw_TB(TB_params,fname)
    type(parameters_TB),intent(out) :: TB_params
    character(len=*), intent(in)    :: fname
    integer     ::  io
    
    io=open_file_read(fname)
    Call read_TB_H       (io,fname,TB_params%io_H)
    Call read_TB_EF      (io,fname,TB_params%io_ef)
    Call read_TB_dos     (io,fname,TB_params%io_dos)
    Call read_TB_flow    (io,fname,TB_params%flow)
    Call read_TB_highs   (io,fname,TB_params%io_highs)
    Call read_TB_occ_mult(io,fname,TB_params%io_occ_mult)
    call close_file(fname,io)
end subroutine

subroutine read_TB_highs(io,fname,highs)
    use m_convert
    integer,intent(in)                          :: io
    character(len=*), intent(in)                :: fname
    type(parameters_TB_IO_highs),intent(out)    :: highs

    real(8),allocatable     ::  tmp(:,:)

    call get_parameter(io,fname,'N_highsym',highs%N_highsym)
    if(highs%N_highsym < 1)then
        return
    endif
    allocate(highs%k_highs(3,highs%N_highsym),source=0.0d0)
    allocate(tmp(highs%N_highsym,3),source=0.0d0)  !necessary because the order of reading in get_parameter
    call get_parameter(io,fname,'k_highs_pts',highs%N_highsym,3,tmp)
    highs%k_highs=transpose(tmp)
    deallocate(tmp)
    allocate(highs%k_highs_frac(highs%N_highsym),source=1)
    call get_parameter(io,fname,'k_highs_frac',highs%N_highsym,highs%k_highs_frac)
    call get_parameter(io,fname,'k_highs_dist',highs%aim_dist)
end subroutine 



subroutine read_TB_flow(io,fname,flow)
    integer,intent(in)                      :: io
    character(len=*), intent(in)            :: fname
    type(parameters_TB_IO_flow),intent(out) :: flow

    call get_parameter(io,fname,'do_TB_r',flow%do_r)
    call get_parameter(io,fname,'do_dos_r',flow%dos_r)
    call get_parameter(io,fname,'do_occ_r',flow%occ_r)
    call get_parameter(io,fname,'do_spec_r',flow%spec_r)
    call get_parameter(io,fname,'do_fermi_r',flow%fermi_r)
    call get_parameter(io,fname,'do_occ_mult_r',flow%occ_mult_r)

    call get_parameter(io,fname,'TB_read_solution_r',flow%read_solution_r)
    call get_parameter(io,fname,'TB_write_solution_r',flow%write_solution_r)

    call get_parameter(io,fname,'do_TB_k',flow%do_k)
    call get_parameter(io,fname,'do_dos_k',flow%dos_k)
    call get_parameter(io,fname,'do_fermi_k',flow%fermi_k)
    call get_parameter(io,fname,'do_highs_k',flow%highs_k)
end subroutine

subroutine read_TB_dos(io,fname,io_dos)
    integer,intent(in)                      :: io
    character(len=*), intent(in)            :: fname
    type(parameters_TB_IO_dos),intent(out)  :: io_dos

    call get_parameter(io,fname,'dos_sigma',io_dos%sigma)
    call get_parameter(io,fname,'dos_E_ext',io_dos%E_ext)
    call get_parameter(io,fname,'dos_dE',io_dos%dE)
    call get_parameter(io,fname,'dos_kgrid',io_dos%kgrid)
end subroutine


subroutine read_TB_occ_mult(io,fname,io_occ_mult)
    integer,intent(in)                              :: io
    character(len=*), intent(in)                    :: fname
    type(parameters_TB_IO_OCC_MULT),intent(out)     :: io_occ_mult

    call get_parameter(io, fname, 'occ_mult_dE', io_occ_mult%dE)
    call get_parameter(io, fname, 'occ_mult_E_ext',2, io_occ_mult%E_ext)
    call get_parameter(io, fname, 'occ_mult_kt',  io_occ_mult%kt)
end subroutine


subroutine read_TB_EF(io,fname,io_ef)
    integer,intent(in)                      :: io
    character(len=*), intent(in)            :: fname
    type(parameters_TB_IO_EF),intent(out)   :: io_ef

    call get_parameter(io, fname, 'N_electrons', io_ef%N_electrons)
    call get_parameter(io, fname, 'fermi_kt', io_ef%kt)
    call get_parameter(io, fname, 'TB_EF', io_ef%E_F_in)
end subroutine
    
subroutine read_TB_H(io,fname,TB_params)
    use m_convert
    use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
    integer,intent(in)                       :: io
    character(len=*), intent(in)             :: fname
    type(parameters_TB_IO_H),intent(out)     :: TB_params
    ! Internal variables
    integer :: N

    Call number_Hpar(io,fname,'TB_hopping',N)
    if(N>0)then
        allocate(TB_params%hop_io(N))
        Call read_Hpar(io,fname,'TB_hopping',TB_params%hop_io)
    else
        STOP "FAILED TO READ tight-binding hopping Hamiltonian (TB_hopping), but orbitals are set"
    endif

    Call number_Hpar(io,fname,'TB_delta',N)
    if(N>0)then
        allocate(TB_params%del_io(N))
        Call read_Hpar(io,fname,'TB_delta',TB_params%del_io)
    else
        write(output_unit,'(/A/)') "No tight-binding superconducting delta found"
    endif

    Call number_Hpar(io,fname,'TB_scfdelta',N)
    if(N>0)then
        allocate(TB_params%del_scf_io(N))
        Call read_Hpar(io,fname,'TB_scfdelta',TB_params%del_scf_io)
    else
        write(output_unit,'(/A/)') "No tight-binding self-consistent superconducting delta found"
    endif

    Call number_Hpar(io,fname,'TB_Jsd',N)
    if(N>0)then
        allocate(TB_params%Jsd(N))
        Call read_Hpar(io,fname,'TB_Jsd',TB_params%Jsd)
    else
        write(output_unit,'(/A/)') "No tight-binding Jsd-coupling found"
    endif

    Call number_Hpar(io,fname,'TB_defect',N)
    if(N>0)then
        allocate(TB_params%defect(N))
        Call read_Hpar(io,fname,'TB_defect',TB_params%defect)
    else
        write(output_unit,'(/A/)') "No tight-binding defect found"
    endif


    call get_parameter(io,fname,'TB_scf_print',TB_params%scf_print)
    call get_parameter(io,fname,'TB_scf_loopmax',TB_params%scf_loopmax)
    call get_parameter(io,fname,'TB_scf_diffconv',TB_params%scf_diffconv)
    call get_parameter(io,fname,'TB_scf_Ecut',TB_params%scf_Ecut)
    call get_parameter(io,fname,'TB_scf_kgrid',TB_params%scf_kgrid)
    call get_parameter(io,fname,'TB_sparse',TB_params%sparse)
    call get_parameter(io,fname,'TB_diag',TB_params%i_diag)
    call get_parameter(io,fname,'TB_diag_acc',TB_params%diag_acc)
    call get_parameter(io,fname,'TB_diag_Ebnd',TB_params%Ebnd)
    call get_parameter(io,fname,'TB_diag_Emin',TB_params%Ebnd(1))
    call get_parameter(io,fname,'TB_diag_Emax',TB_params%Ebnd(2))
    if(TB_params%Ebnd(1)>=TB_params%Ebnd(2))then
        write(error_unit,'(2/A/2(E16.8/))') "WARNING, tight binding minimal energy bound is smaller than maximal energy bound:", TB_params%Ebnd
        STOP "Fix input"
    endif
    call get_parameter(io,fname,'TB_diag_estNe',TB_params%estNe)
end subroutine

subroutine number_Hpar(io,fname,var_name,Nnonzero)
    use, intrinsic :: iso_fortran_env, only : output_unit

    integer, intent(in)                         :: io
    character(len=*), intent(in)                :: fname,var_name
    integer,intent(out)                         :: Nnonzero

    ! internal variable
    integer :: Nentry
    integer :: i
    integer :: stat
    logical :: success
    class(TB_H_par),allocatable :: tmp
    character(len=100) :: str

    Nentry=0; Nnonzero=0
    Call set_pos_entry(io,fname,var_name,success)
    if(.not.success) return
    read(io,*) !read line with var_name, there is no additional information

    !We start to read the input
    !Find out how many entries there are 
    Call alloc_TB_H(tmp,var_name)
    do 
        read(io,'(a)',iostat=stat) str
        if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
        read(str,*,iostat=stat) tmp
        if (stat /= 0) exit
        Nentry=Nentry+1
        if(.not.tmp%is_zero()) Nnonzero=Nnonzero+1
    enddo
    Call write_info_number_found(Nentry,Nnonzero,var_name)
    success=.true.
    !rewind to read actual data
    do i=1,Nentry+1
        backspace(io)
    enddo
end subroutine

subroutine read_Hpar(io,fname,var_name,par)
    use, intrinsic :: iso_fortran_env, only : output_unit
    integer, intent(in)                         :: io
    character(len=*), intent(in)                :: fname,var_name
    class(TB_H_par), intent(inout)              :: par(:)

    ! internal variable
    integer :: ii
    integer :: stat
    character(len=100) :: str

    ii=1
    do while (ii<=size(par))
        read(io,'(a)',iostat=stat) str
        read(str,*,iostat=stat) par(ii)
        if(par(ii)%is_zero()) cycle
        write(output_unit,'(2A,I6,A)') var_name,' entry no.',ii,':'
        Call par(ii)%print_std()
        ii=ii+1
    enddo 
    Call check_further_entry(io,fname,var_name)
end subroutine
end module  m_rw_TB
