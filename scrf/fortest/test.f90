program test


  use netcdf !netcdf
  use utim   !time utility routines
  use nrtype ! Numerical recipies types
  use linkstruct !structure from topnet model for grid information
  use gridweight !grid structure used by spcorr
  use nr, only: erf, erfcc ! Numerical Recipies error function
  use namelist_module_rndnum, only: read_namelist_rndnum !namelist module
  use namelist_module_rndnum, only: start_date, stop_date, start_ens, stop_ens, cross_cc_flag,  weight_judge, &
    & exp2p_file, cross_file_prefix, cc_file, grid_name, out_spcorr_prefix, out_rndnum_prefix


  implicit none
  
  ! ######################################################################################
  ! start interface
  interface

 	subroutine save_rndnum (rndnum, file, error)
      use netcdf
      use nrtype
      real (sp), intent (in) :: rndnum (:, :, :)
      character (len=500), intent (in) :: file
      integer, intent (out) :: error
    end subroutine save_rndnum
 	
 	subroutine field_rand_nopointer (nspl1, nspl2, sp_wght, sp_sdev, sp_ipos, sp_jpos, sp_num, iorder1d, jorder1d, cfield)
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
  		integer (i4b), intent (in) :: iorder1d(:), jorder1d(:)
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
    
    subroutine read_nc_3Dvar (file_name, var_name, var_out, error)
      use netcdf
      use nrtype
      character (len=*), intent (in) :: file_name
  	  character (len=*), intent (in) :: var_name
      real (sp), allocatable, intent (out) :: var_out (:, :, :)
      integer, intent (out) :: error
    end subroutine read_nc_3Dvar
    
    subroutine spcorr_grd_exp2p (nspl1, nspl2, c0m, s0m, grid, weight_judge, sp_wght_var, sp_sdev_var, sp_ipos_var, sp_jpos_var, sp_num_var, iorder1d, jorder1d)
      use nrtype 
  	  use nr, only: ludcmp, lubksb 
      use nrutil, only: arth
      use trig_degrees, only: sind, cosd
      use linkstruct
      use gridweight
      implicit none
      integer (i4b), intent (in) :: nspl1 
  	  integer (i4b), intent (in) :: nspl2
  	  real (dp), intent (in) :: c0m(:, :)
	  real (dp), intent (in) :: s0m(:, :)
	  real (dp), intent (in) :: weight_judge
      type (coords), intent (in) :: grid
      real (dp), intent (out) :: sp_wght_var(:,:,:)
      real (dp), intent (out) :: sp_sdev_var(:,:)
      integer (i4b), intent (out) :: sp_ipos_var(:,:,:), sp_jpos_var(:,:,:)
      integer (i4b), intent (out) :: sp_num_var(:,:)
      integer (i4b), intent (out) :: iorder1d(:), jorder1d(:)
    end subroutine spcorr_grd_exp2p
    
    subroutine generate_date_series (start_date, stop_date, date_series)
      implicit none
      integer, intent(in):: start_date, stop_date ! start/stop date (yyyymmdd)
      integer, intent(out), allocatable :: date_series(:, :) ! (month numbers, year/month/start_day/stop_day)
    end subroutine generate_date_series

  end interface
  ! end interface
  ! ######################################################################################
 
  ! Local variables
  integer (i4b) :: i, j, k, mstep, istep, iens, mm, yy, ntimes !counter variables
  character(10) ::  mmstr, yymmstr
  integer (i4b) :: ierr, jerr !error variables for various error checks
  integer (i4b) :: nspl1 ! # points (1st spatial dimension)
  integer (i4b) :: nspl2 ! # points (2nd spatial dimension)
 
  integer :: f ! AWW for command line argument read
  character (len=200) :: namelist_filename !AWW now an argument to the program
  character (len=1024) :: arg ! AWW command line arg for configuration file
  character (len=1024) :: file_spcc_struct ! file name of spatial correlation structure
  character (len=1024) :: file_scrf ! file name of output SCRF (random number)
  character (len=1024) :: file_scrf_forcross ! file name of outside scrf for cross correlation purpose
  character (len=128) :: suffix !suffix for ensemble member output
  real (dp), allocatable :: lon_out (:), lat_out (:), hgt_out (:) ! lon/lat/height output to netcdf
  real (dp), allocatable :: lat (:, :), lon (:, :), hgt (:, :), slp_e (:, :), slp_n (:, :), mask (:, :) ! grid information read from gridfile
  real (dp), allocatable :: c0 (:, :, :), s0 (:, :, :), c0m (:, :), s0m (:, :) ! correlation model parameters read from exp2p file
  real (sp), allocatable :: cc_lc (:, :, :), cc_lcm (:, :) ! lag1 or cross correlation
  
  real (sp), allocatable :: rndnum_3D (:, :, :) ! scrf (x, y, time)
  real (sp), allocatable :: rndnum_3D_forcross (:, :, :) ! scrf from outside files for cross correlation purpose
  real (dp), allocatable :: rndnum_2D (:, :), old_random(:, :) ! scrf (x, y)
  
  integer (i4b) :: nx, ny !grid size
  integer (i4b) :: spl1_start, spl2_start !starting point of x,y grid
  integer (i4b) :: spl1_count, spl2_count !length of x,y grid
  integer (i4b) :: error

  integer               :: year, day, hour, minute, second, monnum  !dates
  character(len=2)      :: mnth_str    !string to contain current month we're in
  integer, allocatable :: date_series(:, :), datesize(:)

  ! spatial correlation structure
  type (coords), pointer :: grid 
  real (dp), dimension (:, :, :), allocatable :: sp_wght_var
  real (dp), dimension (:, :), allocatable :: sp_sdev_var
  integer (i4b), dimension (:,:,:), allocatable :: sp_ipos_var, sp_jpos_var
  integer (i4b), dimension (:,:), allocatable ::sp_num_var
  integer (i4b), dimension (:), allocatable :: iorder1d, jorder1d
  integer (i4b) :: maxprev = 49 ! hard coded following spcorr_grd
  
  logical :: file_flag
  integer (i4b) :: initflag ! whether it is the first time step when generating scrf

  real :: a, b
  
  a=0.348
  print *, 'erfinv ',erfinv(a)
   
end program test
