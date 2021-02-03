program generate_rndnum_exp2p
! Creator: Guoqiang Tang, 2020
! Purpose: Produce spatiotemporally correlated random numbers (SCRF) for all days and 
! all grids in a single run
! spcorr is no longer pointer to enable explicit memory control

! -----------------------------------------------------------------------------
! Creator(s):
!   Andy Newman, 2013
! Modified:
!   Andy Wood, 2016 -- adding namelist file as argument
!                   -- no longer hardwired; clean formatting
!                   -- adding documentation
!                   -- altered namelist args to specify ens range
!   Guoqiang Tang, 2021 -- adopts spatially varied parameters from Exponential 2-parameter
! -----------------------------------------------------------------------------
! Purpose:
!   Driver for spatially correlated random field code from Martyn Clark
!   Generates ensebles of precipitation and temperature from regression step
!   For version 0 of CONUS ensemble product.  See Newman et al. 2015 JHM
! -----------------------------------------------------------------------------
 
  use netcdf !netcdf
  use utim   !time utility routines
  use nrtype ! Numerical recipies types
  use linkstruct !structure from topnet model for grid information
  use gridweight !grid structure used by spcorr
  use nr, only: erf, erfcc ! Numerical Recipies error function
  use namelist_module_rndnum, only: read_namelist_rndnum !namelist module
  use namelist_module_rndnum, only: start_date, stop_date, start_ens, stop_ens, cross_cc_flag,  &
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
    
    subroutine spcorr_grd_exp2p (nspl1, nspl2, c0m, s0m, grid,  sp_wght_var, sp_sdev_var, sp_ipos_var, sp_jpos_var, sp_num_var, iorder1d, jorder1d)
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
  real (dp) :: ctime_start, ctime_end ! computation time

  
  ! ========== code starts below ==============================
  

  ! ######################################################################################
  ! read file: part-1
  ! read parameter/file settings from input the namelist file
  f = 0
  do
    call get_command_argument (f, arg)
    if (f .eq. 1) namelist_filename = arg
    if (len_trim(arg) == 0) exit
    f = f + 1
  end do

  call read_namelist_rndnum (namelist_filename)
  exp2p_file = trim(exp2p_file)
  grid_name = trim(grid_name)
  out_spcorr_prefix = trim(out_spcorr_prefix)
  out_rndnum_prefix = trim(out_rndnum_prefix)
  
  print *, 'start_date is ', start_date
  print *, 'stop_date is ', stop_date
  print *, 'start_ens is ', start_ens
  print *, 'stop_ens is ', stop_ens
  print *, 'cross_cc_flag is', cross_cc_flag
!   print *, 'exp2p_file is ', exp2p_file
!   print *, 'cross_file_prefix is ', cross_file_prefix
!   print *, 'cc_file is ', cc_file
!   print *, 'grid_name is ', grid_name
!   print *, 'out_spcorr_prefix is ', out_spcorr_prefix
!   print *, 'out_rndnum_prefix is ', out_rndnum_prefix
  
  ! Generate date series (four columns: year/month/start_day/stop_day)
  call generate_date_series(start_date, stop_date, date_series)
  datesize = shape(date_series)
  monnum = datesize(1)
  print *, 'start date', date_series(1, :)
  print *, 'stop date', date_series(monnum, :)

  ! ######################################################################################
  ! read file: part-2
  ! read grid file which contains information such as lat/lon/elevation etc (i.e., lat, 
  ! lon, hgt, slp_n, slp_e, mask, nx, ny)

  ! initialize error status
  error = 0
  ierr = 0
  jerr = 0

  !read in netcdf grid file
  call read_nc_grid (grid_name, lat, lon, hgt, slp_n, slp_e, mask, nx, ny, error)
  if (error .ne. 0) call exit_scrf (1, 'problem in read_nc_grid ')
  
  ! generate lat/lon vector only used for output purpose
  allocate (lat_out(nx*ny), lon_out(nx*ny), hgt_out(nx*ny), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating for 1-d output variables')
  lon_out = pack (lon, .true.)
  lat_out = pack (lat, .true.)
  hgt_out = pack (hgt, .true.)
  
  ! set/name a few variables for generating spatial correlation structure (may be simplified in future versions)
  nspl1 = nx
  nspl2 = ny
  spl1_start = 1
  spl2_start = 1
  spl1_count = nx
  spl2_count = ny
  
  ! allocate space for a pointer that stores the grid information (only lat/lon/elevation) for later use
  nullify (grid)
  allocate (grid, stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating structure grid')
  grid%idx%spl1_start = spl1_start
  grid%idx%spl2_start = spl2_start
  grid%idx%spl1_count = spl1_count
  grid%idx%spl2_count = spl2_count
  allocate (grid%lat(spl1_count, spl2_count), grid%lon(spl1_count, spl2_count), grid%elv(spl1_count, spl2_count), stat=jerr)
  if (ierr .ne. 0 .or. jerr .ne. 0) call exit_scrf (1, ' problem allocating space for lat-lon-elev coordinates ')
  grid%lat = lat
  grid%lon = lon
  grid%elv = hgt
  
  ! ######################################################################################
  ! read file: part-3
  ! read exp2p parameters
  error = 0
  call read_nc_exp2p (exp2p_file, c0, s0, error)
  if (error .ne. 0) call exit_scrf (1, 'problem in read_nc_grid ')
  
  ! read auto correlation CC_lag1 or cross correlation CC_cross
  if (cross_cc_flag .lt. 0) then
      call read_nc_3Dvar (cc_file, 'CC_lag1', cc_lc, error)
  else
  	  call read_nc_3Dvar (cc_file, 'CC_cross', cc_lc, error)
  end if
  if (error .ne. 0) call exit_scrf (1, 'problem in reading cc using read_nc_3Dvar')
  
  ! ######################################################################################
  ! Produce spatial correlation structure based on the correlation model (e.g., Exponential 2parameter)
  ! The structure will be used in random number generation
  print *, 'Generating spatial correlation structure for 12 months'
  do mm = 1, 12
  	  call cpu_time(ctime_start) ! computation time: start 
  
      print *,'Processing month', mm
      ! check if output file exists
      write( mmstr, '(i2)' )  mm
	  file_spcc_struct =  trim(out_spcorr_prefix) // 'month_'// trim(mmstr)
	  
	  INQUIRE(FILE=file_spcc_struct, EXIST=file_flag)
	  if (file_flag) then
	  	  print *, 'Outfile exists. Continue to next month.'
	  else
		  ! generate structure: spcorr
		  if (associated(spcorr)) deallocate (spcorr, stat=ierr)
		  c0m = c0(:,:,mm)
		  s0m = s0(:,:,mm)
		  call spcorr_grd_exp2p (nspl1, nspl2, c0m, s0m, grid, sp_wght_var, sp_sdev_var, sp_ipos_var, sp_jpos_var, sp_num_var, iorder1d, jorder1d)
  		  
		  ! save structure information to files
		  open(unit=34,file= file_spcc_struct,form='unformatted',iostat=error)
		  if(error .ne. 0) then; print *, 'Error opening station weight file', error; stop; end if
		  write(unit=34,iostat=error) sp_wght_var, sp_ipos_var, sp_jpos_var, sp_num_var, sp_sdev_var, iorder1d, jorder1d
		  if(error .ne. 0) then; print *, 'Error writing station weight file ', error; stop; end if
		  close(unit=34)
		  
		  deallocate (sp_wght_var, sp_ipos_var, sp_jpos_var, sp_num_var, sp_sdev_var, iorder1d, jorder1d)
	  end if
	  
	  call cpu_time(ctime_end) ! computation time: end 
      print *, 'computation time (seconds): ', ctime_end - ctime_start
  end do ! end loop mm
 
  
  ! ######################################################################################
  ! Generate SCRF
  
  ! deallocate/allocate
  if (allocated(rndnum_2D))  deallocate (rndnum_2D, stat=ierr)
  allocate (rndnum_2D(nspl1, nspl2), stat=ierr)

  if (allocated(old_random)) deallocate (old_random, stat=ierr)
  allocate (old_random(nspl1, nspl2), stat=ierr)
  
  if (allocated(sp_wght_var)) deallocate (sp_wght_var)
  if (allocated(sp_ipos_var)) deallocate (sp_ipos_var)
  if (allocated(sp_jpos_var)) deallocate (sp_jpos_var)
  if (allocated(sp_num_var)) deallocate (sp_num_var)
  if (allocated(sp_sdev_var)) deallocate (sp_sdev_var) 
  if (allocated(iorder1d)) deallocate (iorder1d) 
  if (allocated(jorder1d)) deallocate (jorder1d) 
  allocate (sp_wght_var(nspl1,nspl2,maxprev), sp_ipos_var(nspl1,nspl2,maxprev), & 
    & sp_jpos_var(nspl1,nspl2,maxprev), sp_num_var(nspl1,nspl2), sp_sdev_var(nspl1,nspl2), & 
    & iorder1d(nspl1*nspl2), jorder1d(nspl1*nspl2), stat=ierr)

  
  ! start generation
  print *, 'Generating SCRF for every ensemble member and every month'
  do iens = start_ens, stop_ens
    initflag = 1 ! different ensemble members are independent with each other
    do mstep = 1, monnum
        call cpu_time(ctime_start) ! computation time: start 
        
        yy = date_series(mstep, 1)
        mm = date_series(mstep, 2)
        ntimes = date_series(mstep, 4) - date_series(mstep, 3) + 1
        
        print *, 'Generating random numbers: member/year/month', iens, yy, mm
        write( mmstr, '(i2)' )  mm
        write( yymmstr, '(i6)' )  yy*100+mm
		write (suffix, '(I3.3)') iens
		cc_lcm = cc_lc(:, :, mm)
		
        ! load spatial correlation structure
        print *, 'load spatial correlation structure'
        file_spcc_struct = trim(out_spcorr_prefix) // 'month_'// trim(mmstr)
		open(unit=34,file=file_spcc_struct,form='unformatted',iostat=error)
	    read(unit=34,iostat=error) sp_wght_var, sp_ipos_var, sp_jpos_var, sp_num_var, sp_sdev_var, iorder1d, jorder1d
		close(unit=34)
		
		! if using cross correlation is true, reading reference scrf from files
		if (cross_cc_flag .gt. 0) then
			file_scrf_forcross = trim(cross_file_prefix) // trim(yymmstr) // '_' // trim(suffix) // '.nc'
			call read_nc_3Dvar (file_scrf_forcross, 'rndnum', rndnum_3D_forcross, error)
        	if (error .ne. 0) call exit_scrf (1, 'problem in reading random numbers for cross correlation')
        end if
		
		! generate random numbers for every day
		print *, 'generate random numbers'	
	    if (allocated(rndnum_3D))  deallocate (rndnum_3D, stat=ierr)
		allocate (rndnum_3D(nx, ny, ntimes), stat=ierr)
		rndnum_3D = 0.0
		
		do istep = 1, ntimes   ! loop for all days in one month
		   print *, 'current time step', istep, '    ///    total steps', ntimes 
		   call field_rand_nopointer (nspl1, nspl2, sp_wght_var, sp_sdev_var, sp_ipos_var, sp_jpos_var, sp_num_var, iorder1d, jorder1d, rndnum_2D)
		   if (initflag .eq. 1) then
		   	  rndnum_3D(:, :, istep) = rndnum_2D
			  initflag = 0
		   else
              old_random = rndnum_2D
      		  ! update temporally isolated random number (rndnum_2D) using random number from previous time step or from another variable
      		  do i = 1, nspl1
			    do j = 1, nspl2
			      if (cross_cc_flag .lt. 0) then
      		      	rndnum_3D(i,j,istep)=old_random(i,j)*cc_lcm(i,j) + sqrt(1 - cc_lcm(i,j)*cc_lcm(i,j))*rndnum_2D(i,j)
      		      else
      		      	rndnum_3D(i,j,istep)=rndnum_3D_forcross(i,j,istep)*cc_lcm(i,j) + sqrt(1 - cc_lcm(i,j)*cc_lcm(i,j))*rndnum_2D(i,j)
      		      end if
      		    end do
      		  end do
		   end if
		end do ! end loop days

		! save random numbers to netcdf files
		print *, 'save random numbers to nc file'
		file_scrf = trim(out_rndnum_prefix) // trim(yymmstr) // '_' // trim(suffix) // '.nc'
		call save_rndnum (rndnum_3D, file_scrf, ierr)
        
        call cpu_time(ctime_end) ! computation time: end 
        print *, 'computation time (seconds): ', ctime_end - ctime_start
    end do !end month loop
  end do !end ensemble member loop
   
 
end program generate_rndnum_exp2p
