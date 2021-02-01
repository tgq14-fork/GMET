subroutine save_rndnum (rndnum, file, error)
  ! only for saving random numbers using nf90_float to save space
  ! latitude/longitude is not saved to save space
  
  use netcdf
  use nrtype
  implicit none
  
  ! in/out
  real (sp), intent (in) :: rndnum (:, :, :)
  character (len=500), intent (in) :: file
  integer, intent (out) :: error

  ! Dimension names
  character (len=*), parameter :: y_name = "y"
  character (len=*), parameter :: x_name = "x"
  character (len=*), parameter :: time_name = "time"
  
  ! Variable Names (lat/lon is not saved)
  character (len=*), parameter :: lat_name = "latitude"
  character (len=*), parameter :: lon_name = "longitude"
  character (len=*), parameter :: rnd_name = "rndnum"
  
  ! others
  integer (i4b), allocatable :: matshp(:) ! shape of input matrix (i.e., rndnum)
  integer (i4b) :: nx, ny, ntimes
  integer :: ncid, x_dimid, y_dimid, time_dimid
  integer :: lat_varid, lon_varid, rnd_varid
  integer :: count2 (2), start2 (2), count3 (3), start3 (3), dimids2 (2), dimids3 (3)
	
  ! matrix shape
  matshp = shape(rndnum)
  nx = matshp(1)
  ny = matshp(2)
  ntimes = matshp(3)

  ! Create NEW output file
  call check (nf90_create(file, NF90_NETCDF4, ncid), "File creation error", error)
  if (error /= 0) return

  ! Define the dimensions.
  call check (nf90_def_dim(ncid, y_name, ny, y_dimid), "y dim def error", error)
  call check (nf90_def_dim(ncid, x_name, nx, x_dimid), "x dim def error", error)
  call check (nf90_def_dim(ncid, time_name, nf90_unlimited, time_dimid), "time dim def error", error)
  if (error /= 0) return

  ! Define the variables.
  ! dimids2 = (/ x_dimid, y_dimid /)
  ! call check (nf90_def_var(ncid, lat_name, nf90_float, dimids2, lat_varid, deflate_level=9), "lat var def error", error)
  ! call check (nf90_def_var(ncid, lon_name, nf90_float, dimids2, lon_varid, deflate_level=9), "lon var def error", error)
  dimids3 = (/ x_dimid, y_dimid, time_dimid /)
  call check (nf90_def_var(ncid, rnd_name, nf90_float, dimids3, rnd_varid, deflate_level=9), "pcp rndnum", error)
  if (error /= 0) return

  ! End define mode.
  call check (nf90_enddef(ncid), "end define mode error", error)
  if (error /= 0) return
  
  ! Write variables
  count2 = (/ nx, ny /) 
  start2 = (/ 1, 1 /)
  count3 = (/ nx, ny, ntimes /)
  start3 = (/ 1, 1, 1 /)
  ! call check (nf90_put_var(ncid, lat_varid, grdlat, start=start2, count=count2), "put lat error", error)
  ! call check (nf90_put_var(ncid, lon_varid, grdlon, start=start2, count=count2), "put lon error", error)
  call check (nf90_put_var(ncid, rnd_varid, rndnum, start=start3, count=count3), "put rndnum error", error)
  if (error /= 0) return

  ! Close netcdf file
  call check (nf90_close(ncid), "closing file error", error)


contains
  subroutine check (status, info, error)
    integer, intent (in) :: status
    character (len=*), intent (in) :: info
    integer, intent (out) :: error

    if (status /= nf90_noerr) then
      print *, trim (info) // ": " // trim (nf90_strerror(status))
      error = 1
    end if
  end subroutine check
end subroutine save_rndnum
