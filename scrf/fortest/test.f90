program test

  use netcdf !netcdf
  use utim   !time utility routines
  use nrtype ! Numerical recipies types
  use linkstruct !structure from topnet model for grid information
  use gridweight !grid structure used by spcorr
  use nr, only: erf, erfcc ! Numerical Recipies error function
  use namelist_module_rndnum, only: read_namelist_rndnum !namelist module
  use namelist_module_rndnum, only: start_time, ntimes, start_ens, stop_ens, exp2p_file, grid_name, out_spcorr_prefix, out_rndnum_prefix
  
  implicit none
  
 
  ! ######################################################################################
  ! start interface
  interface

 	subroutine save_rndnum (pcp_rndnum, tmean_rndnum, trange_rndnum, nx, ny, ntimes, grdlat, grdlon, grdalt, file, error)
      use netcdf
      use nrtype
      real (sp), intent (in) :: pcp_rndnum (:, :, :), tmean_rndnum (:, :, :), trange_rndnum (:, :, :)
      integer (i4b), intent (in) :: nx, ny, ntimes
      real (dp), intent (in) :: grdlat (:), grdlon (:), grdalt (:)
      character (len=500), intent (in) :: file
      integer, intent (out) :: error
    end subroutine save_rndnum
 	
 	subroutine field_rand_nopointer (nspl1, nspl2, sp_wght, sp_sdev, sp_ipos, sp_jpos, sp_num, cfield)
  		use nrtype ! variable types (DP, I4B, etc.)
  		use nr, only: gasdev ! Num. Recipies
  		use inputdat2d ! use to relate basins to gridpoints
  		implicit none
  		integer (i4b), intent (in) :: nspl1 ! # points (1st spatial dimension)
  		integer (i4b), intent (in) :: nspl2 ! # points (2nd spatial dimension)
  		real (dp), intent (in) :: sp_wght(:, :, :)
  		real (dp), intent (in) :: sp_sdev(:, :)
  		integer (i4b), intent (in) :: sp_ipos(:,:,:), sp_jpos(:,:,:)
  		integer (i4b), intent (in) ::sp_num(:,:)
  		real (dp), dimension (nspl1, nspl2), intent (out) :: cfield ! correlated random field
    end subroutine field_rand_nopointer
 
    function erfinv (x)
      use nrtype
      real (sp), intent (in) :: x
      real (sp) :: erfinv
    end function erfinv
    
    subroutine read_nc_grid (file_name, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
      use netcdf
      use nrtype
      character (len=*), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), grad_e (:, :), mask (:, :)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_nc_grid
 
    subroutine read_nc_exp2p (file_name, c0, s0, error)
      use netcdf
      use nrtype
      character (len=*), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: c0 (:, :, :), s0 (:, :, :)
      integer, intent (out) :: error
    end subroutine read_nc_exp2p

  end interface
  ! end interface
  ! ######################################################################################
 
  ! Local variables
  integer (i4b) :: i, j, k, igrd, istep, iens, mm !counter variables
  character(10) ::  mmstr
  integer (i4b), dimension (1:2) :: order1 = (/ 2, 1 /)!order for reshape array
  integer (i4b) :: ierr, jerr !error variables for various error checks
  integer (i4b) :: nspl1 ! # points (1st spatial dimension)
  integer (i4b) :: nspl2 ! # points (2nd spatial dimension)
  integer (i4b) :: isp1  ! first grid dimension location
  integer (i4b) :: isp2  ! second grid dimension location
  real (dp), dimension (:, :), allocatable :: rho ! temporal correlation parameter
  real (dp), dimension (:, :), allocatable :: old_random ! previous correlated random field
  real (dp), dimension (:, :), allocatable :: pcp_random ! new correlated random field for pcp
  real (dp), dimension (:, :), allocatable :: tmean_random ! new correlated rand field, tmean
  real (dp), dimension (:, :), allocatable :: trange_random ! new correlated rand field, trange
 
  real (sp) :: acorr !value from scrf
  real (sp) :: aprob !probability from scrf
  real (sp) :: a_ra
  real (sp) :: aprob_ra
 
  real (dp) :: cprob !cdf value from scrf
  real (dp) :: amult !multiplier value to get actual precip from normalized value
  real (dp) :: rn
  real (dp) :: ra
  real (dp) :: ra_err
  real (dp) :: cs
  real (dp) :: cprob_ra
  integer(I4B)   :: cs_percentile !for climo precip distribution
 
  real (dp) :: transform
 
  integer :: f ! AWW for command line argument read
  character (len=200) :: namelist_filename !AWW now an argument to the program
  character (len=1024) :: arg ! AWW command line arg for configuration file
  character (len=1024) :: out_name !base output name for netcdf files
  character (len=128) :: suffix !suffix for ensemble member output
  character (len=1024) :: var_name !name of netcdf variable grabbed from jason's netcdf file
  real (dp), allocatable :: lon_out (:)! lon output to netcdf
  real (dp), allocatable :: lat_out (:)! lat output to netcdf
  real (dp), allocatable :: hgt_out (:)! hgt output to netcdf
  real (dp), allocatable :: lat (:, :), lon (:, :)
  real (dp), allocatable :: c0 (:, :, :), s0 (:, :, :), c0m (:, :), s0m (:, :)
  real (dp), allocatable :: hgt (:, :)
  real (dp), allocatable :: slp_e (:, :)
  real (dp), allocatable :: slp_n (:, :)
  real (dp), allocatable :: mask (:, :)
  real (dp), allocatable :: weight (:, :)!weights from spcorr
  real (dp), allocatable :: std (:, :)!std from spcorr
  real (dp), allocatable :: var (:, :, :)!generic variable
  real (dp), allocatable :: pcp (:, :, :)!output from qpe code, normalized precip
  real (dp), allocatable :: pop (:, :, :)!output from qpe code, normalized pop
  real (dp), allocatable :: pcp_error (:, :, :)!error from ols regression in qpe code
  real (dp), allocatable :: tmean (:, :, :)
  real (dp), allocatable :: tmean_error (:, :, :)
  real (dp), allocatable :: trange (:, :, :)
  real (dp), allocatable :: trange_error (:, :, :)
  real(DP)               :: obs_max !precipitation limit using observations (used as cap for ensemble members)

  real (dp), allocatable :: lons (:, :)!lons array from qpe code
  real (dp), allocatable :: lats (:, :)!lats array from qpe code
  real (dp), allocatable :: times (:)!time vector from qpe code
  real (dp), allocatable :: auto_corr (:)!lag-1 autocorrelation vector from qpe code
  real (dp), allocatable :: tpc_corr (:)!temp-precip correlation vector from qpe code
  real (dp), allocatable :: obs_max_pcp (:, :, :) !max of non-0 pcp (each tstep)
  ! Add by TGQ
  real (sp), allocatable :: pcp_rndnum (:, :, :)!
  real (sp), allocatable :: pcp_cprob (:, :, :)!
  real (sp), allocatable :: tmean_rndnum (:, :, :)!
  real (sp), allocatable :: trange_rndnum (:, :, :)!
  ! Add by TGQ
  integer (i4b) :: nx, ny !grid size
  integer (i4b) :: spl1_start, spl2_start !starting point of x,y grid
  integer (i4b) :: spl1_count, spl2_count !length of x,y grid
  integer (i4b) :: tot_times
  integer (i4b) :: ncid, dimid, varid, error
  integer (i4b) :: nTimesRegression 

  !climo grid variables
  real(SP),allocatable  :: climo_tmin(:,:,:)        !monthly climo tmax grids
  real(SP),allocatable  :: climo_tmax(:,:,:)        !monthly climo tmin grids
  real(SP),allocatable  :: climo_precip(:,:,:)        !climo precip grid for current day
  real(SP),allocatable  :: climo_tmean(:,:)         !climo tmean grid for current day
  real(SP),allocatable  :: climo_trange(:,:)        !climo trange grid for current day
  real(SP),allocatable  :: uncert_tmean(:,:)         !uncertainty tmean grid for current month
  real(SP),allocatable  :: uncert_trange(:,:)        !uncertainty trange grid for current month
  logical               :: first_climo = .FALSE.    !logical for first read
  integer               :: prev_month = -999        !previous month
  integer               :: current_month            !current month
  integer               :: prev_delta               !number of days from previous month for temporal interpolation
  integer               :: next_delta               !number of days to next month for temporal interpolation
  integer               :: year, day, hour, minute, second  !dates
  integer,dimension(12) :: month_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  character(len=2)      :: mnth_str    !string to contain current month we're in
  character(len=1024)   :: climo_file
  real(DP)              :: combined_error           !total error of daily anomaly uncertainty and climo uncertainty
  real(SP)              :: max_pcp                  !maximum allowable precip for a grid cell
  ! add by TGQ. Regressed using station data for North America (km).
  ! fit using Pearson correlation coefficient using (cij = exp(-dij/Clen))
!   real(DP),dimension(12) :: clen_daily_prcp= (/290.25, 278.77, 246.90, 219.92, 183.17, 143.70, 111.21, 119.47, 183.84, 241.29, 263.62, 284.26/)
  real(DP),dimension(12) :: clen_daily_tmean= (/1058.53, 1216.52, 1137.26, 962.64, 895.36, 763.74, 629.39, 650.82, 996.21, 1129.24, 1250.76, 1152.05/)
  real(DP),dimension(12) :: clen_daily_trange= (/344.97, 352.29, 349.27, 378.64, 366.32, 346.70, 295.97, 291.39, 388.96, 483.01, 403.84, 333.49/)

  ! prcp only: fit using Spearman correlation coefficient using (cij = exp(-dij/Clen))
  real(DP),dimension(12) :: clen_daily_prcp= (/350.65, 339.90, 320.07, 306.85, 268.04, 234.21, 185.38, 185.22, 273.86, 344.70, 351.25, 346.14/)
  
  ! Auto corr
  real(DP),dimension(12) :: auto_corr_daily= (/0.629, 0.636, 0.648, 0.634, 0.650, 0.620, 0.559, 0.561, 0.633, 0.654, 0.652, 0.643/)
  real(DP),dimension(12) :: tp_corr_daily= (/-0.167, -0.221, -0.246, -0.277, -0.279, -0.261, -0.229, -0.249, -0.287, -0.266, -0.199, -0.151/)
  ! add by TGQ

  type (coords), pointer :: grid !coordinate structure for grid
  type (splnum), dimension (:, :), pointer :: sp_pcp, sp_temp, sp_trange ! structures of spatially correlated random field weights
  
  real (dp), dimension (:, :, :), allocatable :: sp_wght_prcp, sp_wght_tmean, sp_wght_trange
  real (dp), dimension (:, :), allocatable :: sp_sdev_prcp, sp_sdev_tmean, sp_sdev_trange
  integer (i4b), dimension (:,:,:), allocatable :: sp_ipos_prcp, sp_jpos_prcp, sp_ipos_tmean, sp_jpos_tmean, sp_ipos_trange, sp_jpos_trange
  integer (i4b), dimension (:,:), allocatable ::sp_num_prcp, sp_num_tmean, sp_num_trange
  logical :: file_exists
  integer (i4b) :: initflag, shp1, shp2, shp3
  integer (i4b), allocatable :: shp(:)
  ! ========== code starts below ==============================
   
  f = 0
  do
    call get_command_argument (f, arg)
    if (f .eq. 1) namelist_filename = arg
    if (len_trim(arg) == 0) exit
    f = f + 1
  end do

  ! read namelist in
  call read_namelist_rndnum (namelist_filename)
  exp2p_file = trim(exp2p_file)
  grid_name = trim(grid_name)
  out_spcorr_prefix = trim(out_spcorr_prefix)
  out_rndnum_prefix = trim(out_rndnum_prefix)
  
  print *, start_time
  print *, ntimes
  print *, start_ens
  print *, stop_ens
  print *, exp2p_file
  print *, grid_name
  print *, out_spcorr_prefix
  print *, out_rndnum_prefix
  
  error = 0
  call read_nc_exp2p (exp2p_file, c0, s0, error)
  if (error .ne. 0) call exit_scrf (1, 'problem in read_nc_grid ')
  shp = shape(c0)
  shp2 = shp(2)
  print *, 'shape', shp2, shp
  
 
end program test
