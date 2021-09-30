module m_anisotropy_heisenberg
use, intrinsic :: iso_fortran_env, only : output_unit, error_unit
use m_derived_types, only: lattice
use m_input_H_types, only: io_H_aniso
private
public :: get_anisotropy_H,read_anisotropy_input, get_anisotropy_fft

character(len=*),parameter  :: ham_desc="magnetic anisotropy"

contains

subroutine read_anisotropy_input(io_unit,fname,io)
    use m_input_H_types, only: io_H_aniso
    use m_io_utils, only: get_parameter
    integer,intent(in)              :: io_unit
    character(len=*), intent(in)    :: fname
    type(io_H_aniso),intent(out)    :: io
    real(8),allocatable     :: tmp(:,:)
    integer                 :: i,j
    logical                 :: success

    !read parameters for anisotropy in cartesian coordinates assuming 4 entries, where 1:3 is unit vector direction and 4 is magnitude
    Call read_int_realarr(io_unit,fname,'magnetic_anisotropy',4,io%attype,io%val,success)

    !read parameters for anisotropy in cartesian coordinates using 3-real format where is magnitude is given by the norm of the entry with the sign of the first nonzero entry
    if(.not.success)then
        Call read_int_realarr(io_unit,fname,'magnetic_anisotropy',3,io%attype,tmp)
        if(allocated(tmp))then
            allocate(io%val(4,size(tmp,2)))
            do i=1,size(tmp,2)
                io%val(4,i)=norm2(tmp(:,i))
                io%val(1:3,i)=tmp(:,i)
                do j=1,3
                    if(tmp(j,i)/=0.0d0)then
                        io%val(4,i)=sign(io%val(4,i),tmp(j,i))
                    endif
                enddo
            enddo
            deallocate(tmp)
        endif
    endif

    !read parameters for anisotropy in normalized lattice parameters
    Call read_int_realarr(io_unit,fname,'magnetic_anisotropy_lat',4,io%attype_lat,io%val_lat)
    if(allocated(io%val_lat))then
        if(any(norm2(io%val_lat(1:3,:),dim=1)==0))then
            write(error_unit,'(/A)') 'ERROR, entry of "magnetic_anisotropy_lat" has vanishing real array(1:3) components which encode the direction' 
            STOP "CHECK INPUT"
        endif
    endif

    if(allocated(io%val).or.allocated(io%val_lat)) io%is_set=.true.
    
    Call get_parameter(io_unit,fname,'magnetic_anisotropy_fft',io%fft) 
end subroutine

subroutine read_int_realarr(io_unit,fname,var_name,Nreal,ints,reals,success)
    !reads lines after variable name with 1 integer and Nreal reals and return all non-vanishing entries allocated to the correct size
    use m_io_utils
    integer,intent(in)                  :: io_unit
    character(len=*)                    :: var_name
    character(len=*)                    :: fname
    integer,intent(in)                  :: Nreal
    integer,allocatable,intent(inout)   :: ints(:)
    real(8),allocatable,intent(inout)   :: reals(:,:)
    logical,optional,intent(out)        :: success  !return if reading was successfull

    ! internal variable
    integer :: int_tmp
    real(8) :: vec_tmp(Nreal)

    integer :: Nentry,Nnonzero
    integer :: nread,check, length_string
    integer :: i,ii
    integer :: stat
    character(len=100) :: str

    if(allocated(ints)) deallocate(ints)
    if(allocated(reals)) deallocate(reals)
    if(present(success)) success=.false.

    nread=0
    length_string=len_trim(var_name)
    rewind(io_unit)
    do
        read (io_unit,'(a)',iostat=stat) str
        if (stat /= 0) exit
        str= trim(adjustl(str))
        if (len_trim(str)==0) cycle
        if (str(1:1) == '#' ) cycle
        if (len_trim(str)/=length_string) cycle

        !We start to read the input
        if ( str(1:length_string) == var_name(1:length_string)) then
            Nentry=0; Nnonzero=0
            nread=nread+1
            do 
                read(io_unit,'(a)',iostat=stat) str
                if (stat /= 0) STOP "UNEXPECTED END OF INPUT FILE"
                read(str,*,iostat=stat) int_tmp,vec_tmp
                if (stat /= 0) exit
                Nentry=Nentry+1
                if(norm2(vec_tmp)/=0.0d0) Nnonzero=Nnonzero+1
            enddo
            if(present(success).and.Nentry<1)then
                !failed to read, but success is supplied so no complaining
                return
            endif
            if(Nentry<1)then
                write(error_unit,'(/2A/A/)') "Found no entries for ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT            
                ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
                return
            endif
            if(Nnonzero<1)then
                write(error_unit,'(/2A/A/)') "Found no nonzero entries for ",var_name,' although the keyword is specified'
#ifndef CPP_SCRIPT            
                ERROR STOP "INPUT PROBABLY WRONG (disable with CPP_SCRIPT preprocessor flag)"
#endif
                return
            endif
            write(output_unit,'(/A,I6,2A)') "Found ",Nnonzero," nonzero entries for Hamiltonian ",var_name
            allocate(ints(Nnonzero),source=0)
            allocate(reals(Nreal,Nnonzero) ,source=0.0d0)
            do i=1,Nentry+1
                backspace(io_unit)
            enddo
            ii=0
            write(str,'(I16)') Nreal
            do i=1,Nentry
                read(io_unit,*,iostat=stat) int_tmp,vec_tmp
                if(norm2(vec_tmp)==0.0d0) cycle
                ii=ii+1
                ints(ii)=int_tmp
                reals(:,ii)=vec_tmp
                write(output_unit,'(2A,I6,A)') var_name,' entry no.',ii,':'
                write(output_unit,'(A,I6)')    '  atom type:', ints(ii)
                write(output_unit,'(A,'//trim(adjustl(str))//'F16.8)') '  vector   :', reals(:,ii)
            enddo 
            write(output_unit,'(/)')
        endif
    enddo
    check=check_read(nread,var_name,fname)
    if(present(success)) success=.true.
end subroutine

subroutine get_anisotropy_H(Ham,io,lat)
    !get anisotropy in t_H Hamiltonian format
    use m_H_public, only: t_H
    use m_setH_util,only: get_coo
    use m_mode_public

    class(t_H),intent(inout)    :: Ham
    type(io_H_aniso),intent(in) :: io
    type(lattice),intent(in)    :: lat
    !local 
    integer :: i
    real(8),allocatable :: Htmp(:,:)
    real(8),allocatable :: val_tmp(:)
    integer,allocatable :: ind_tmp(:,:)
    integer,allocatable :: connect(:,:)

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting Hamiltonian: ", ham_desc
        !set local Hamiltonian 
        allocate(Htmp(lat%M%dim_mode,lat%M%dim_mode),source=0.d0)
        !add anisotropy input in cartesian input coordinates to unit-cell Hamiltonian
        Call get_Haniso_unitcell(io,lat,Htmp)            
        !get local Hamiltonian in coo format
        Call get_coo(Htmp,val_tmp,ind_tmp)
        !Anisotropy only has simple onsite terms
        allocate(connect(2,lat%Ncell))
        do i=1,lat%Ncell
            connect(:,i)=i
        enddo
        Call Ham%init_connect(connect,val_tmp,ind_tmp,"MM",lat,1)
        Ham%desc=ham_desc
        
        !set modes
        Call mode_set_rank1(Ham%mode_l,lat,"M")
        Call mode_set_rank1(Ham%mode_r,lat,"M")
    endif
end subroutine

subroutine get_anisotropy_fft(H_fft,io,lat)
    !get anisotropy in the fft_H format. 
    !Since the anisotropy is super localized in normal space this makes absolutely no sense unless used in combination with a delocalized Hamiltonian (dipolar-interaction) so that the evaluation is for free 
    use m_fft_H_base, only: fft_H
    use m_mode_public

    class(fft_H),intent(inout)  :: H_fft
    type(io_H_aniso),intent(in) :: io
    type(lattice),intent(in)    :: lat


    integer         :: Nmag             !number of magnetic atoms per unit-cell
    logical         :: period(3)        !consider as periodic or open boundary condition along each direction (T:period, F:open)
                                        ! (dim_lat(i)=1->period(i)=T, since the calculation in the periodic case is easier, but choice of supercell_vec still does not consider periodicity)
    integer         :: N_rep(3)         !number of states in each direction in the fourier transformation
    integer         :: Nk_tot           !number of state considered in FT (product of N_rep)
    integer         :: Kbd(2,3)    !boundaries of the K-operator

    real(8),allocatable :: Karr(:,:,:)  !K-operator to be FT'd (1:3*Nmag,1:3*Nmag,1:Nk_tot)

    !local 

    if(io%is_set)then
        write(output_unit,'(/2A)') "Start setting fft-Hamiltonian: ", ham_desc
        !set some initial parameters locally for convencience
        Nmag=lat%nmag
        period=lat%periodic.or.lat%dim_lat==1

        !set shape-dependent quantities of fft_H and get Kdb,N_rep
        Call H_fft%init_shape(3*lat%nmag,period,lat%dim_lat,Kbd,N_rep)
        Nk_tot=product(N_rep)

        !set local Hamiltonian 
        allocate(Karr(3*Nmag,3*Nmag,Nk_tot),source=0.0d0)
        Call get_Haniso_unitcell(io,lat,Karr(:,:,1))            
        Call H_fft%init_op(3*Nmag,Karr,ham_desc)
    endif
end subroutine

subroutine get_Haniso_unitcell(io,lat,H)
    !translates the io_H_aniso input to the anisotropy contributions as a dense Hamiltonian (H) within the unit-cell
    type(io_H_aniso),intent(in) :: io
    type(lattice),intent(in)    :: lat
    real(8),intent(inout)       :: H(lat%M%dim_mode,lat%M%dim_mode)
    integer             ::  i,j
    integer,allocatable :: ind_at(:)    !atom indices in space of all atoms
    integer             :: ind_mag      !atom index in space of magnetic atoms
    real(8)             :: direction(3),magnitude
    real(8)             :: lat_norm(3,3)
    integer             :: i1,i2

    !add anisotropy input in cartesian input coordinates to unit-cell Hamiltonian
    if(allocated(io%attype))then
        do i=1,size(io%attype)
            Call lat%cell%ind_attype(io%attype(i),ind_at)
            magnitude=io%val(4,i)
            direction=io%val(1:3,i)/norm2(io%val(1:3,i))
            do j=1,size(ind_at)
                ind_mag=lat%cell%ind_mag(ind_at(j))
                do i1=1,3
                    do i2=1,3
                        H((ind_mag-1)*3+i1,(ind_mag-1)*3+i2)=H((ind_mag-1)*3+i1,(ind_mag-1)*3+i2)+magnitude*direction(i1)*direction(i2)
                    enddo
                enddo
            enddo
        enddo
    endif
    !add anisotropy input in real-space lattice coordinates to unit-cell Hamiltonian
    if(allocated(io%attype_lat))then
        lat_norm=transpose(lat%areal)
        do i=1,3
            lat_norm(:,i)=lat_norm(:,i)/norm2(lat_norm(:,i))
        enddo
        do i=1,size(io%attype_lat)
            Call lat%cell%ind_attype(io%attype_lat(i),ind_at)
            magnitude=io%val_lat(4,i)
            direction=io%val_lat(1:3,i)/norm2(io%val_lat(1:3,i))
            direction=matmul(lat_norm,direction)
            do j=1,size(ind_at)
                ind_mag=lat%cell%ind_mag(ind_at(j))
                do i1=1,3
                    do i2=1,3
                        H((ind_mag-1)*3+i1,(ind_mag-1)*3+i2)=H((ind_mag-1)*3+i1,(ind_mag-1)*3+i2)+magnitude*direction(i1)*direction(i2)
                    enddo
                enddo
            enddo
        enddo
    endif
end subroutine


end module m_anisotropy_heisenberg
