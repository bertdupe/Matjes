module m_netcdf_routine
use netcdf
IMPLICIT NONE

private
public :: readgrid,griddims,writegrid,check

contains

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
!Define coordinate variables
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

end module m_netcdf_routine
