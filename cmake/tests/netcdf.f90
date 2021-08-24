program test_netcdf
use netcdf
IMPLICIT NONE
integer :: ncid,x_dimid
integer :: stat

stat=nf90_create("test", NF90_CLOBBER, ncid)
stat=nf90_def_dim(ncid, "lon", 1, x_dimid)
stat=nf90_close(ncid)

!contains
!SUBROUTINE writegrid(outfile,var_name,xpos,ypos,zpos,idata,NX,NY,NZ,time)
!integer, INTENT(IN) :: xpos(:),ypos(:),zpos(:)
!integer, INTENT(IN) :: NX, NY, NZ
!real(8), INTENT(IN) :: idata(:,:,:),time
!CHARACTER(LEN=*), INTENT(IN) :: outfile,var_name
!
!integer :: ncid, x_dimid, y_dimid, z_dimid
!integer :: x_varid, y_varid, z_varid, varid
!integer, DIMENSION(3) :: dimids
!
!!Create the netCDF file.
!CALL check(nf90_create(outfile, NF90_CLOBBER, ncid))
!!Define the dimensions.
!CALL check(nf90_def_dim(ncid, "lon", NX, x_dimid))
!CALL check(nf90_def_dim(ncid, "lat", NY, y_dimid))
!CALL check(nf90_def_dim(ncid, "level", NZ, z_dimid))
!CALL check(nf90_def_var(ncid, "lon", NF90_INT, x_dimid, x_varid))
!CALL check(nf90_def_var(ncid, "lat", NF90_INT, y_dimid, y_varid))
!CALL check(nf90_def_var(ncid, "level", NF90_INT, z_dimid, z_varid))
!dimids = (/ x_dimid, y_dimid, z_dimid/)
!!Define variable
!CALL check(nf90_def_var(ncid, var_name, NF90_FLOAT, dimids, varid))
!CALL check(nf90_enddef(ncid)) !End Definitions
!!Write Data
!CALL check(nf90_put_var(ncid, x_varid, xpos))
!CALL check(nf90_put_var(ncid, y_varid, ypos))
!CALL check(nf90_put_var(ncid, z_varid, zpos))
!CALL check(nf90_put_var(ncid, varid, idata))
!CALL check(nf90_close(ncid))
!END SUBROUTINE 

end program
