subroutine read_nc_exp2p (file_name, c0, s0, error)
 
  use netcdf
  use nrtype
  implicit none

  character (len=*), intent (in) :: file_name

  ! Exponential 2-parameter: r = exp(-(d/c0)**s0)
  real (dp), allocatable, intent (out) :: c0 (:, :, :) ! correlation length parameter
  real (dp), allocatable, intent (out) :: s0 (:, :, :) ! shape parameter
  
  integer, intent (out) :: error

  ! local variables
  integer :: ncid !netcdf file id
  integer :: i
  integer :: varid !variable id
  integer (i4b), dimension (nf90_max_var_dims) :: dimids !dimension ids for dimensions of grid file
  integer (i4b) :: ndims, nlat, nlon, ntime
  character (len=*), parameter :: c0_name = "c0"
  character (len=*), parameter :: s0_name = "s0"


  ! open netcdf grid file
  call check (nf90_open(trim(file_name), nf90_nowrite, ncid), "File open error", error)
  if (error /= 0) return
  
  !!! read c0
  ! inquire variable
  call check (nf90_inq_varid(ncid, c0_name, varid), "c0 name error", error)
  if (error /= 0) return
  ! check dimensions (only check for the first variable)
  call check (nf90_inquire_variable(ncid, varid, ndims=ndims, dimids=dimids), "Dimension inq error", error)
  if (error /= 0 .or. ndims /= 3) return
  ! get x,y dimensions (only check for the first variable)
  call check (nf90_inquire_dimension(ncid, dimids(1), len=nlon), "x dim error", error)
  call check (nf90_inquire_dimension(ncid, dimids(2), len=nlat), "y dim error", error)
  call check (nf90_inquire_dimension(ncid, dimids(3), len=ntime), "z dim error", error)
  if (error /= 0) return
  ! get the variable
  allocate (c0(nlon, nlat, ntime))
  call check (nf90_get_var(ncid, varid, c0), "c0 read error", error)
  if (error /= 0) return
  
  
  !!! read s0
  ! inquire variable
  call check (nf90_inq_varid(ncid, s0_name, varid), "s0 name error", error)
  if (error /= 0) return
  ! get the variable
  allocate (s0(nlon, nlat, ntime))
  call check (nf90_get_var(ncid, varid, s0), "s0 read error", error)
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
end subroutine read_nc_exp2p
