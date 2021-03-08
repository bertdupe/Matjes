module m_anisotropy_heisenberg
private
public :: get_anisotropy_H,read_anisotropy_input

contains

subroutine read_anisotropy_input(io_unit,fname,io)
    use m_io_utils
    use m_input_H_types, only: io_H_aniso
    use, intrinsic :: iso_fortran_env, only : output_unit
    integer,intent(in)              :: io_unit
    character(len=*), intent(in)    :: fname
    type(io_H_aniso),intent(out)    :: io

    character(len=*), parameter     :: var_name='magnetic_anisotropy'
    ! internal variable
    integer :: attype
    real(8) :: vec(3)

    integer :: Nentry,Nnonzero
    integer :: nread,check, length_string
    integer :: i,ii
    integer :: stat
    character(len=100) :: str
    
    nread=0
    length_string=len_trim(var_name)
    rewind(io_unit)
    do
        read (io_unit,'(a)',iostat=stat) str
        if (stat /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle
        if (str(1:1) == '#' ) cycle

        !We start to read the input
        if ( str(1:length_string) == var_name(1:length_string)) then
            Nentry=0; Nnonzero=0
            nread=nread+1
            do 
                read(io_unit,'(a)',iostat=stat) str
                if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
                read(str,*,iostat=stat) attype,vec
                if (stat /= 0) exit
                Nentry=Nentry+1
                if(norm2(vec)/=0.0d0) Nnonzero=Nnonzero+1
            enddo
            if(Nentry<1)then
                write(output_unit,'(/2A/A/)') "Found no entries for ",var_name,' although the keyword is specified'
                ERROR STOP "INPUT PROBABLY WRONG"
            endif
            if(Nnonzero<1)then
                write(output_unit,'(/2A/A/)') "Found no nonzero entries for ",var_name,' although the keyword is specified'
                ERROR STOP "INPUT PROBABLY WRONG"
            endif
            write(output_unit,'(/A,I6,2A)') "Found ",Nnonzero," nonzero entries for Hamiltonian ",var_name
            io%is_set=.true.
            allocate(io%attype(Nnonzero),source=0)
            allocate(io%val(3,Nnonzero) ,source=0.0d0)
            do i=1,Nentry+1
                backspace(io_unit)
            enddo
            ii=1
            do i=1,Nentry
                read(io_unit,*,iostat=stat) attype,vec
                if(norm2(vec)==0.0d0) cycle
                io%attype(ii)=attype
                io%val(:,ii)=vec
                write(output_unit,'(2A,I6,A)') var_name,' entry no.',ii,':'
                write(output_unit,'(A,I6)')    '  atom type:', io%attype(ii)
                write(output_unit,'(A,3F16.8)') '  vector   :', io%val(:,ii)
                ii=ii+1
            enddo 
            write(output_unit,'(/)')
        endif
    enddo

    check=check_read(nread,var_name,fname)
end subroutine

subroutine get_anisotropy_H(Ham,io,lat)
    !get anisotropy in t_H Hamiltonian format
    use m_H_public, only: t_H
    use m_derived_types, only: lattice
    use m_setH_util,only: get_coo
    use m_input_H_types, only: io_H_aniso

    class(t_H),intent(inout)    :: Ham
    type(io_H_aniso),intent(in) :: io
    type(lattice),intent(in)    :: lat
    !local 
    integer :: i,j
    real(8),allocatable :: Htmp(:,:)
    real(8),allocatable :: val_tmp(:)
    integer,allocatable :: ind_tmp(:,:)
    integer,allocatable :: ind_at(:)    !atom indices in space of all atoms
    integer             :: ind_mag      !atom index in space of magnetic atoms
    integer,allocatable :: connect(:,:)

    if(io%is_set)then
        !set local Hamiltonian 
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode),source=0.d0)
        do i=1,size(io%attype)
            Call lat%cell%ind_attype(io%attype(i),ind_at)
            do j=1,size(ind_at)
                ind_mag=lat%cell%ind_mag(ind_at(j))
                Htmp((ind_mag-1)*3+1,(ind_mag-1)*3+1)=io%val(1,i)
                Htmp((ind_mag-1)*3+2,(ind_mag-1)*3+2)=io%val(2,i)
                Htmp((ind_mag-1)*3+3,(ind_mag-1)*3+3)=io%val(3,i)
            enddo
        enddo
        !get local Hamiltonian in coo format
        Call get_coo(Htmp,val_tmp,ind_tmp)
        !Anisotropy only has simple onsite terms
        allocate(connect(2,lat%Ncell))
        do i=1,lat%Ncell
            connect(:,i)=i
        enddo
        Call Ham%init_connect(connect,val_tmp,ind_tmp,"MM",lat,1)
        Ham%desc="magnetic anisotropy"
    endif
end subroutine

end module m_anisotropy_heisenberg
