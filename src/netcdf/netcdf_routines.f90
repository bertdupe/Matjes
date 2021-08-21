module m_netcdf_routine
#ifdef CPP_NETCDF
use netcdf
#endif
use, intrinsic  ::  ISO_FORTRAN_ENV, only: error_unit, output_unit
IMPLICIT NONE

private
#ifdef CPP_NETCDF
public :: readgrid,griddims,writegrid,check
#endif
public :: netcdf_type

type netcdf_type
    integer     :: ncid=0    !netcdf id 
    integer     :: dims(3)  =[0,0,0]   !number of entries in each direction 
    integer     :: dim_id(3)=[0,0,0]   !netcdf-id for positions
    integer     :: var_id              !netcdf-id for actual variable 
    integer     :: N_entries=0  !number of times the state has been written (does not do anything as there is no time yet)
    
    logical     :: is_open=.false.  !file is open
contains
    procedure   :: file_open
    procedure   :: file_write
    procedure   :: file_close
end type

contains

subroutine file_open(this,fname,var_name,dims)
    class(netcdf_type),intent(inout)    :: this
    character(len=*),intent(in)         :: fname
    character(len=*),intent(in)         :: var_name
    integer,intent(in)                  :: dims(3)
#ifdef CPP_NETCDF    

    integer     :: xpos(dims(1)),ypos(dims(2)), zpos(dims(3))
    integer     :: pos_id(3), i

    if(this%is_open)then
        write(error_unit,'(2/,A)')  "Trying to open an opened netcdf-file"
        ERROR STOP
    endif

    CALL check(nf90_create(fname, NF90_CLOBBER, this%ncid))
    this%dims=dims
    this%is_open=.true.

    !Define the dimensions.
    CALL check(nf90_def_dim(this%ncid, "lon",   this%dims(1),   this%dim_id(1)))
    CALL check(nf90_def_dim(this%ncid, "lat",   this%dims(2),   this%dim_id(2)))
    CALL check(nf90_def_dim(this%ncid, "level", this%dims(3),   this%dim_id(3)))
!    CALL check(nf90_def_dim(this%ncid, "time",  NF90_UNLIMITED, this%dim_id(4)))

    !Define variables
    !!coordinate variables with regular integer position array
    CALL check(nf90_def_var(this%ncid, "lon"  , NF90_INT, this%dim_id(1), pos_id(1)))
    CALL check(nf90_def_var(this%ncid, "lat"  , NF90_INT, this%dim_id(2), pos_id(2)))
    CALL check(nf90_def_var(this%ncid, "level", NF90_INT, this%dim_id(3), pos_id(3)))
    !!coordinate variables with regular integer position array

    !TODO add actual variable

    CALL check(nf90_enddef(this%ncid)) !End Definitions

    !write position array
    xpos=[(i,i=1,this%dims(1))] 
    ypos=[(i,i=1,this%dims(2))] 
    zpos=[(i,i=1,this%dims(3))] 

    CALL check(nf90_put_var(this%ncid, this%dim_id(1), xpos))
    CALL check(nf90_put_var(this%ncid, this%dim_id(2), ypos))
    CALL check(nf90_put_var(this%ncid, this%dim_id(3), zpos))
#else
    write(error_unit,'(/A)') "WARNING, called netcdf routine, while Matjes is called without CPP_NETCDF"
    write(error_unit,'(A/)') "         this will have no effect"
#endif
end subroutine    

subroutine file_close(this)
    class(netcdf_type),intent(inout)    :: this
    if(.not.this%is_open)then
        write(error_unit,'(2/,A)')  "Trying to close to an unopened netcdf-file"
        ERROR STOP
    endif
    this%is_open=.false.
#ifdef CPP_NETCDF
    CALL check(nf90_close(this%ncid))
#else
    write(error_unit,'(/A)') "WARNING, called netcdf routine, while Matjes is called without CPP_NETCDF"
    write(error_unit,'(A/)') "         this will have no effect"
#endif
end subroutine

subroutine file_write(this,dat)
    class(netcdf_type),intent(inout)    :: this
    real(8),intent(in)                  :: dat(:)

    if(.not.this%is_open)then
        write(error_unit,'(2/,A)')  "Trying to write to an unopened netcdf-file"
        ERROR STOP
    endif

    !TODO, write to the actual data 
    STOP "NETCDF file write, and most other things, is not implemented"

    this%N_entries=this%N_entries+1
end subroutine


#ifdef CPP_NETCDF
!WRITEGRID - write a netCDF gridfile
!:=========================================================================
SUBROUTINE writegrid(outfile,var_name,xpos,ypos,zpos,idata,NX,NY,NZ,time)
integer, INTENT(IN) :: xpos(:),ypos(:),zpos(:)
integer, INTENT(IN) :: NX, NY, NZ
real(8), INTENT(IN) :: idata(:,:,:),time
CHARACTER(LEN=*), INTENT(IN) :: outfile,var_name

integer :: ncid, x_dimid, y_dimid, z_dimid
integer :: x_varid, y_varid, z_varid, varid
integer, DIMENSION(3) :: dimids

!Create the netCDF file.
CALL check(nf90_create(outfile, NF90_CLOBBER, ncid))
!Define the dimensions.
CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
CALL check(nf90_def_dim(ncid, "level", NZ, z_dimid))
CALL check(nf90_def_var(ncid, "lon", NF90_INT, x_dimid, x_varid))
CALL check(nf90_def_var(ncid, "lat", NF90_INT, y_dimid, y_varid))
CALL check(nf90_def_var(ncid, "level", NF90_INT, z_dimid, z_varid))
dimids = (/ x_dimid, y_dimid, z_dimid/)
!Define variable
CALL check(nf90_def_var(ncid, var_name, NF90_FLOAT, dimids, varid))
CALL check(nf90_enddef(ncid)) !End Definitions
!Write Data
CALL check(nf90_put_var(ncid, x_varid, xpos))
CALL check(nf90_put_var(ncid, y_varid, ypos))
CALL check(nf90_put_var(ncid, z_varid, zpos))
CALL check(nf90_put_var(ncid, varid, idata))
CALL check(nf90_close(ncid))

END SUBROUTINE writegrid
!:=========================================================================


!GRIDDIMS - Get dimensions of a netCDF 2D gridfile
!:======================================================================
SUBROUTINE griddims(infile,NX,NY)
INTEGER, INTENT(OUT) :: NX, NY
INTEGER :: ncid
CHARACTER(LEN=50), INTENT(IN) :: infile
CHARACTER(LEN=50) :: xname, yname

!Open netCDF file
!:-------:-------:-------:-------:-------:-------:-------:
CALL check(nf90_open(infile, nf90_nowrite, ncid))
!Inquire about the dimensions
!:-------:-------:-------:-------:-------:-------:-------:
CALL check(nf90_inquire_dimension(ncid,1,xname,NX))
CALL check(nf90_inquire_dimension(ncid,2,yname,NY))
!Close netCDF file
!:-------:-------:-------:-------:-------:-------:-------:
CALL check(nf90_close(ncid))
END SUBROUTINE griddims
!:======================================================================


!READGRID - read a netCDF gridfile
!:========================================================================
SUBROUTINE readgrid(infile,xpos,ypos,idata,NX,NY)
REAL(KIND=4), DIMENSION(NX), INTENT(OUT) :: xpos
REAL(KIND=4), DIMENSION(NY), INTENT(OUT) :: ypos
REAL(KIND=4), DIMENSION(NX,NY), INTENT(OUT) :: idata
INTEGER(KIND=4), INTENT(IN) :: NX, NY
INTEGER(KIND=4), DIMENSION(2) :: dimids
INTEGER(KIND=4) :: ncid, xtype, ndims, varid
CHARACTER(LEN=50), INTENT(IN) :: infile
CHARACTER(LEN=50) :: xname, yname, vname
!Open netCDF file
!:-------:-------:-------:-------:-------:-------:-------:-------:
CALL check(nf90_open(infile, nf90_nowrite, ncid))
!Get the values of the coordinates and put them in xpos & ypos
!:-------:-------:-------:-------:-------:-------:-------:-------:
CALL check(nf90_inquire_variable(ncid,1,vname,xtype,ndims,dimids))
CALL check(nf90_inq_varid(ncid,vname,varid))
CALL check(nf90_get_var(ncid,varid,xpos))
CALL check(nf90_inquire_variable(ncid,2,vname,xtype,ndims,dimids))
CALL check(nf90_inq_varid(ncid,vname,varid))
CALL check(nf90_get_var(ncid,varid,ypos))
!Get the values of the perturbations and put them in idata
!:-------:-------:-------:-------:-------:-------:-------:-------:
CALL check(nf90_inquire_variable(ncid,3,vname,xtype,ndims,dimids))
CALL check(nf90_inq_varid(ncid,vname,varid))
CALL check(nf90_get_var(ncid,varid,idata))
!Close netCDF file
!:-------:-------:-------:-------:-------:-------:-------:-------:
CALL check(nf90_close(ncid))
END SUBROUTINE readgrid
!:=========================================================================



!Check
!:======================================================================
SUBROUTINE check(istatus)
INTEGER, INTENT (IN) :: istatus
IF (istatus /= nf90_noerr) THEN
 write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
END IF
END SUBROUTINE check
!:======================================================================

#endif
end module m_netcdf_routine
