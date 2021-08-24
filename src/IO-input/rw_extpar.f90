module m_rw_extpar
use m_input_types, only: extpar_input
implicit none
public
contains
subroutine rw_extpar(extpar_io,areal,io_in,fname_in)
    use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
    use m_io_files_utils, only : open_file_read,close_file
    use m_io_utils,only: get_parameter
    type(extpar_input),intent(out)      :: extpar_io
    real(8),intent(in)                  :: areal(3,3)   !real-space lattice parameters
    integer,intent(in),optional         :: io_in
    character(*),intent(in),optional    :: fname_in
    !internal
    character(*),parameter              :: fname_default='input'
    character(:), allocatable           :: fname
    integer                             :: io_param

    real(8) :: lat_norm(3,3)
    real(8) :: ext_lat(4)
    integer :: i


    if(present(fname_in))then
        fname=fname_in
    else
        fname=fname_default
    endif
    if(present(io_in))then
        io_param=io_in
    else
        io_param=open_file_read(fname)
    endif
    call get_parameter(io_param,fname,'H_ext',extpar_io%H)
    call get_parameter(io_param,fname,'E_ext',extpar_io%E)
    call get_parameter(io_param,fname,'Tini',extpar_io%T(1))
    call get_parameter(io_param,fname,'Tfin',extpar_io%T(2))
    call get_parameter(io_param,fname,'enable_H',extpar_io%enable_H)
    call get_parameter(io_param,fname,'enable_E',extpar_io%enable_E)
    call get_parameter(io_param,fname,'enable_T',extpar_io%enable_T)
    call get_parameter(io_param,fname,'enable_M',extpar_io%enable_M)
    call get_parameter(io_param,fname,'enable_u',extpar_io%enable_u)
    call get_parameter(io_param,fname,'enable_phonon',extpar_io%enable_u)
    call get_parameter(io_param,fname,'enable_w',extpar_io%enable_w)

    lat_norm=transpose(areal)
    do i=1,3
        lat_norm(:,i)=lat_norm(:,i)/norm2(lat_norm(:,i))
    enddo
    ext_lat=0.0d0
    call get_parameter(io_param,fname,'H_ext_lat',ext_lat)
    if(any(ext_lat/=0.0d0))then
        if(norm2(ext_lat(1:3))==0.0d0)then
            write(error_unit,'(/A)') 'ERROR, found "H_ext_lat" input flag, but first 3 entries are zero so that no direction can be determined'
            write(error_unit,'(/A)') 'Cannot initialize magnetic field as direction in space of lattice vectors is undefined'
            STOP
        endif
        ext_lat(1:3)=ext_lat(1:3)/norm2(ext_lat(1:3))
        write(output_unit,'(A)') 'Found input "H_ext_lat", adding external magnetic field'
        write(output_unit,'(3X,A,3F16.8)') 'Direction is space on lattice vectors: ',ext_lat(1:3)
        write(output_unit,'(3X,A,F16.8,A)')  'Magnitude: ',ext_lat(4), ' T'
        ext_lat(1:3)=matmul(lat_norm,ext_lat(1:3))
        extpar_io%H=extpar_io%H+ext_lat(1:3)*ext_lat(4)
    endif

    ext_lat=0.0d0
    call get_parameter(io_param,fname,'E_ext_lat',ext_lat)
    if(any(ext_lat/=0.0d0))then
        if(norm2(ext_lat(1:3))==0.0d0)then
            write(error_unit,'(/A)') 'ERROR, found "E_ext_lat" input flag, but first 3 entries are zero so that no direction can be determined'
            write(error_unit,'(/A)') 'Cannot initialize electric field as direction in space of lattice vectors is undefined'
            STOP
        endif
        ext_lat(1:3)=ext_lat(1:3)/norm2(ext_lat(1:3))
        write(output_unit,'(A)') 'Found input "H_ext_lat", adding external electric field'
        write(output_unit,'(3X,A,3F16.8)') 'Direction is space on lattice vectors: ',ext_lat(1:3)
        write(output_unit,'(3X,A,F16.8,A)')  'Magnitude: ',ext_lat(4), ' V/nm'
        ext_lat(1:3)=matmul(lat_norm,ext_lat(1:3))
        extpar_io%E=extpar_io%E+ext_lat(1:3)*ext_lat(4)
    endif

    if(.not.present(io_in)) call close_file(fname,io_param)
end subroutine
end module
