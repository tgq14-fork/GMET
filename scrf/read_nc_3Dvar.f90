subroutine read_nc_3Dvar (file_name, var_name, var_out, error)
  ! read 3D variable from netcdf file
 
  use netcdf
  use nrtype
  implicit none

  character (len=*), intent (in) :: file_name
  character (len=*), intent (in) :: var_name
  real (sp), allocatable, intent (out) :: var_out (:, :, :) ! lag-1 auto-correlation 
  integer, intent (out) :: error

  ! local variables
  integer :: ncid !netcdf file id
  integer :: i
  integer :: varid !variable id
  integer (i4b), dimension (nf90_max_var_dims) :: dimids !dimension ids for dimensions of grid file
  integer (i4b) :: ndims, nlat, nlon, ntime


  ! open netcdf grid file
  call check (nf90_open(trim(file_name), nf90_nowrite, ncid), "File open error", error)
  if (error /= 0) return
  
  !!! read var_out
  ! inquire variable
  call check (nf90_inq_varid(ncid, var_name, varid), "var_out name error", error)
  if (error /= 0) return
  ! check dimensions (only check for the first variable)
  call check (nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), "Dimension inq error", error)
  if (error /= 0) return
  ! get x,y dimensions (only check for the first variable)
  call check (nf90_inquire_dimension(ncid, dimids(1), len=nlon), "x dim error", error)
  call check (nf90_inquire_dimension(ncid, dimids(2), len=nlat), "y dim error", error)
  call check (nf90_inquire_dimension(ncid, dimids(3), len=ntime), "z dim error", error)
  if (error /= 0) return
  ! get the variable
  allocate (var_out(nlon, nlat, ntime))
  call check (nf90_get_var(ncid, varid, var_out), "var_out read error", error)
  if (error /= 0) return


  !!! Close the file.
  i = nf90_close (ncid)

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
!
end subroutine read_nc_3Dvar
