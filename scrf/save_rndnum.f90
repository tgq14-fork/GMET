subroutine save_rndnum (pcp_rndnum, tmean_rndnum, trange_rndnum, nx, ny, ntimes, grdlat, grdlon, grdalt, file, error)
  ! only for saving random numbers
  use netcdf
  use nrtype
  implicit none
!
  real (sp), intent (in) :: pcp_rndnum (:, :, :), tmean_rndnum (:, :, :), trange_rndnum (:, :, :)
  integer (i4b), intent (in) :: nx, ny, ntimes
  real (dp), intent (in) :: grdlat (:), grdlon (:), grdalt (:)
  character (len=500), intent (in) :: file
  integer, intent (out) :: error
!
!
  ! Dimension names
  character (len=*), parameter :: y_name = "y"
  character (len=*), parameter :: x_name = "x"
  character (len=*), parameter :: time_name = "time"
!
  ! Variable Names
  character (len=*), parameter :: lat_name = "latitude"
  character (len=*), parameter :: lon_name = "longitude"
  character (len=*), parameter :: alt_name = "elevation"
  character (len=*), parameter :: pcp_name = "pcp"
  character (len=*), parameter :: pop_name = "t_mean"
  character (len=*), parameter :: pcp_error_name = "t_range"
  character (len=*), parameter :: pcprnd_name = "pcp_rndnum"
  character (len=*), parameter :: pcpcprob_name = "pcp_cprob"
  character (len=*), parameter :: tmeanrnd_name = "tmean_rndnum"
  character (len=*), parameter :: trangernd_name = "trange_rndnum"
!
  character (len=*), parameter :: long_name = "long_name"
  character (len=*), parameter :: pcp_long_name = "estimated precip in mm/day"
  character (len=*), parameter :: pop_long_name = "estimated daily mean temperature"
  character (len=*), parameter :: pcp_error_long_name = "estimated diurnal range"
!
  ! Units
  character (len=*), parameter :: units = "units"
  character (len=*), parameter :: pcp_units = "mm"
  character (len=*), parameter :: pop_units = "deg_C"
  character (len=*), parameter :: pcp_error_units = "deg_C"
  character (len=*), parameter :: lat_units = "degrees_north"
  character (len=*), parameter :: lon_units = "degrees_east"
  character (len=*), parameter :: alt_units = "meters"
  character (len=*), parameter :: time_units = "seconds since 1970-01-01 00:00:00.0 0:00"
  character (len=*), parameter :: fill = "_FillValue"
!
  real (dp), allocatable :: file_times (:)
!
  integer :: n_chars, n_times, inx, iny
  integer :: ncid, x_dimid, y_dimid, time_dimid
  integer :: lat_varid, lon_varid, alt_varid, time_varid, pcp_varid, pop_varid, pcp_error_varid, &
 & pcprnd_varid, pcpcprob_varid, tmeanrnd_varid, trangernd_varid
  integer :: count1 (1), start1 (1), count2 (2), start2 (2), count3 (3), start3 (3), dimids2 (2), &
 & dimids3 (3)
  integer :: trec, nrecs, file_nx, file_ny, file_ntimes, i
!
  trec = 0
  n_times = ntimes
  n_chars = 100
  inx = nx
  iny = ny
!
  if (size(grdlat) /= inx*iny) then
    print *, "Error "
  end if

  ! AW initial version allowed appending to file if it existed. 
  ! AW for now desired behavior is to just create a new file  
  !    so a lot of the following code (after the 'else') is commented out

!a  error = nf90_open (file, nf90_write, ncid)
!a  if (error /= nf90_noerr) then
!a    error = 0

    ! Create NEW output file
    call check (nf90_create(file, NF90_NETCDF4, ncid), "File creation error", error)
    if (error /= 0) return

    ! Define the dimensions.
    call check (nf90_def_dim(ncid, y_name, iny, y_dimid), "y dim def error", error)
    call check (nf90_def_dim(ncid, x_name, inx, x_dimid), "x dim def error", error)
    call check (nf90_def_dim(ncid, time_name, nf90_unlimited, time_dimid), "time dim def error", &
   & error)
    if (error /= 0) return

    ! Define the variables.
    dimids2 = (/ x_dimid, y_dimid /)
    call check (nf90_def_var(ncid, lat_name, nf90_double, dimids2, lat_varid, deflate_level=9), "lat var def error", &
   & error)
    call check (nf90_def_var(ncid, lon_name, nf90_double, dimids2, lon_varid, deflate_level=9), "lon var def error", &
   & error)
    call check (nf90_def_var(ncid, alt_name, nf90_double, dimids2, alt_varid, deflate_level=9), "alt var def error", &
   & error)
    if (error /= 0) return
   
   ! Add by TGQ
   dimids3 = (/ x_dimid, y_dimid, time_dimid /)
    call check (nf90_def_var(ncid, pcprnd_name, nf90_float, dimids3, pcprnd_varid, deflate_level=9), "pcp rndnum", &
   & error)
    call check (nf90_def_var(ncid, tmeanrnd_name, nf90_float, dimids3, tmeanrnd_varid, deflate_level=9), "tmean rndnum", &
   & error)
   call check (nf90_def_var(ncid, trangernd_name, nf90_float, dimids3, trangernd_varid, deflate_level=9), "trange rndnum", &
   & error)
   ! Add by TGQ
    if (error /= 0) return



    ! End define mode.
    call check (nf90_enddef(ncid), "end define mode error", error)
    if (error /= 0) return



    trec = 1
    nrecs = 0

  count3 = (/ inx, iny, n_times /)
  start3 = (/ 1, 1, trec /)
  ! Add by TGQ
  call check (nf90_put_var(ncid, pcprnd_varid, pcp_rndnum, start=start3, count=count3), "put pcp rndnum", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, tmeanrnd_varid, tmean_rndnum, start=start3, count=count3), "put tmean rndnum", error)
  if (error /= 0) return
  call check (nf90_put_var(ncid, trangernd_varid, trange_rndnum, start=start3, count=count3), "put trange rndnum", error)
  if (error /= 0) return
  ! Add by TGQ
  
!
  call check (nf90_close(ncid), "closing file error", error)
!
contains
  subroutine check (status, info, error)
    integer, intent (in) :: status
    character (len=*), intent (in) :: info
    integer, intent (out) :: error
!
    if (status /= nf90_noerr) then
      print *, trim (info) // ": " // trim (nf90_strerror(status))
      error = 1
    end if
  end subroutine check
end subroutine save_rndnum
