module m_write_config
!module with write_config, 
!which writes out current state of the different order parameters in sererate files
!, if the respective order-parameter is associated

!use m_io_utils
use m_io_files_utils, only: open_file_write,close_file
use m_io_utils, only: dump_config
use m_type_lattice,only: lattice,order_parameter_name,order_par
use m_netcdf_routine
implicit none


interface write_config
 module procedure write_config_char_real
 module procedure write_config_char
 module procedure write_config_int
end interface

private
public write_config
#ifdef CPP_NETCDF
public write_netcdf     !is this really needed?
#endif 
contains

!write the configuration in the Netcdf file
!:========================================================================
#ifdef CPP_NETCDF
subroutine write_netcdf(fname,lat,time)
    ! write out all order parameters that are set with the respective name and simulations
    ! parameters in a netcdf
    character(len=*),intent(in) ::  fname
    type(lattice),intent(in)    ::  lat
    real(8), intent(in)         ::  time

    integer     :: i,xcoord(lat%dim_lat(1)),ycoord(lat%dim_lat(2)),zcoord(lat%dim_lat(3))
    character(len=100) :: filename
    real(8),pointer,contiguous ::  ord(:),ord_netcdf(:,:,:)
    logical         ::  used(size(order_parameter_name))

    Call lat%used_order(used)

    write(filename,'(2A)')  trim(adjustl(fname)),'.nc'

    xcoord=(/(i,i=1,lat%dim_lat(1))/)
    ycoord=(/(i,i=1,lat%dim_lat(2))/)
    zcoord=(/(i,i=1,lat%dim_lat(3))/)
    do i=1,size(order_parameter_name)
        if(used(i))then
            Call lat%set_order_point(i,ord)
            allocate(ord_netcdf(lat%dim_lat(1),lat%dim_lat(2),lat%dim_lat(3)))
            ord_netcdf=reshape(ord,(/lat%dim_lat(1),lat%dim_lat(2),lat%dim_lat(3)/))
            call writegrid(filename,order_parameter_name(i),xcoord,ycoord,zcoord,ord_netcdf,lat%dim_lat(1),lat%dim_lat(2),lat%dim_lat(3),time)
        endif
    enddo

end subroutine
#endif
!:========================================================================

subroutine write_config_char(fname,lat)
    !write out all order parameters that are set with the respective name in front
    !and with the respective dimmode of columns
    character(len=*),intent(in) ::  fname
    type(lattice),intent(in)    ::  lat

    integer     :: io,i
    character(len=100) :: filename
    real(8),pointer,contiguous ::  ord(:)
    logical         ::  used(size(order_parameter_name))

    Call lat%used_order(used)
    do i=1,size(order_parameter_name)
        if(used(i))then
            write(filename,'(4A)')  trim(order_parameter_name(i)),'_',trim(adjustl(fname)),'.dat'
            io=open_file_write(filename)
            Call lat%set_order_point(i,ord)
            call dump_config(io,lat%get_order_dim(i),ord)
            call close_file(filename,io)
        endif
    enddo
end subroutine

subroutine write_config_int(i,lat)
    integer,intent(in)          ::  i
    type(lattice),intent(in)    ::  lat

    character(len=100) :: fname

    write(fname,*) i
    Call write_config_char(trim(adjustl(fname)),lat)
end subroutine

subroutine write_config_char_real(fname,r,lat,acc_op)
    character(len=*),intent(in) :: fname
    real(8),intent(in)          :: r
    type(lattice),intent(in)    :: lat
    integer,intent(in),optional :: acc_op(2)

    integer     ::  i,j
    integer,parameter   :: acc_par(2)=[6,4]
    integer ::  acc(2)

    character(len=100) :: comb_name
    character(len=100) :: name_format

    if(present(acc_op))then
        acc=acc_op
    else
        acc=acc_par
    endif
    write(name_format,*) '(A,A,I0.',acc(1),',A,I0.',acc(2),')'
    i=int(r)
    j=nint((r-real(int(r)))*10**(acc(2)))
    write(comb_name,name_format) Trim(adjustl(fname)),'_',i,'_',j
    Call write_config_char(trim(adjustl(comb_name)),lat)
end subroutine


end module
