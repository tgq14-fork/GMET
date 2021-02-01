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
  use namelist_module_rndnum, only: start_time, ntimes, start_ens, stop_ens, exp2p_file, grid_name, out_spcorr_prefix, out_rndnum_prefix


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
  integer (i4b) :: i, j, k, istep, iens, mm !counter variables
  character(10) ::  mmstr
  integer (i4b) :: ierr, jerr !error variables for various error checks
  integer (i4b) :: nspl1 ! # points (1st spatial dimension)
  integer (i4b) :: nspl2 ! # points (2nd spatial dimension)
 
  integer :: f ! AWW for command line argument read
  character (len=200) :: namelist_filename !AWW now an argument to the program
  character (len=1024) :: arg ! AWW command line arg for configuration file
  character (len=1024) :: file_spcc_struct ! file name of spatial correlation structure
  character (len=1024) :: file_scrf ! file name of SCRF (random number)
  character (len=128) :: suffix !suffix for ensemble member output
  real (dp), allocatable :: lon_out (:), lat_out (:), hgt_out (:) ! lon/lat/height output to netcdf
  real (dp), allocatable :: lat (:, :), lon (:, :), hgt (:, :), slp_e (:, :), slp_n (:, :), mask (:, :) ! grid information read from gridfile
  real (dp), allocatable :: c0 (:, :, :), s0 (:, :, :), c0m (:, :), s0m (:, :) ! correlation model parameters read from exp2p file
  
  real (dp), allocatable :: auto_corr (:)!lag-1 autocorrelation vector from qpe code
  real (sp), allocatable :: rndnum_3D (:, :, :) ! scrf (x, y, time)
  real (sp), allocatable :: rndnum_2D (:, :), old_random(:, :) ! scrf (x, y)
  
  integer (i4b) :: nx, ny !grid size
  integer (i4b) :: spl1_start, spl2_start !starting point of x,y grid
  integer (i4b) :: spl1_count, spl2_count !length of x,y grid
  integer (i4b) :: error

  integer               :: year, day, hour, minute, second  !dates
  integer,dimension(12) :: month_days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  character(len=2)      :: mnth_str    !string to contain current month we're in

  ! spatial correlation structure
  type (coords), pointer :: grid 
  real (dp), dimension (:, :, :), allocatable :: sp_wght_var
  real (dp), dimension (:, :), allocatable :: sp_sdev_var
  integer (i4b), dimension (:,:,:), allocatable :: sp_ipos_var, sp_jpos_var
  integer (i4b), dimension (:,:), allocatable ::sp_num_var
  integer (i4b) :: maxprev = 49 ! max number of previous generated points. Hard coded following spcorr_grd.f90
  
  logical :: file_flag
  integer (i4b) :: initflag ! whether it is the first time step when generating scrf

  
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

  call read_namelist (namelist_filename)
  exp2p_file = trim(exp2p_file)
  grid_name = trim(grid_name)
  out_spcorr_prefix = trim(out_spcorr_prefix)
  out_rndnum_prefix = trim(out_rndnum_prefix)
  
  print *, 'start_time is ', start_time
  print *, 'ntimes is ', ntimes
  print *, 'start_ens is ', start_ens
  print *, 'stop_ens is ', stop_ens
  print *, 'exp2p_file is ', exp2p_file
  print *, 'grid_name is ', grid_name
  print *, 'out_spcorr_prefix is ', out_spcorr_prefix
  print *, 'out_rndnum_prefix is ', out_rndnum_prefix
  
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
  ! read exp2p parameters. this file must be generated following the lat/lon format of grid file
  error = 0
  call read_nc_exp2p (exp2p_file, c0, s0, error)
  if (error .ne. 0) call exit_scrf (1, 'problem in read_nc_grid ')
  
  ! ######################################################################################
  ! Produce spatial correlation structure based on the correlation model (e.g., Exponential 2parameter)
  ! The structure will be used in random number generation
  print *, 'Generating spatial correlation structure for every month'
  do mm = 1, 12
      print *,'Processing month', mm
      
      write( mmstr, '(i2)' )  mm
	  file_spcc_struct =  trim(out_spcorr_prefix) // trim(mmstr)
	  
	  INQUIRE(FILE=file_spcc_struct, EXIST=file_flag)
	  if (file_flag) then
	  	  print *, 'Outfile exists. Continue to next month.'
	  else
	      ! deallocate/allocate spatial correlation structure variables
	      if (allocated(sp_wght_var)) deallocate (sp_wght_var)
	      if (allocated(sp_ipos_var)) deallocate (sp_ipos_var)
	      if (allocated(sp_jpos_var)) deallocate (sp_jpos_var)
	      if (allocated(sp_num_var)) deallocate (sp_num_var)
	      if (allocated(sp_sdev_var)) deallocate (sp_sdev_var)
		  allocate (sp_wght_var(nspl1,nspl2,maxprev), sp_ipos_var(nspl1,nspl2,maxprev), & 
		  	& sp_jpos_var(nspl1,nspl2,maxprev),sp_num_var(nspl1,nspl2), sp_sdev_var(nspl1,nspl2), stat=ierr)
		  
		  ! generate structure: spcorr
		  if (associated(spcorr)) deallocate (spcorr, stat=ierr)
		  c0m = c0(:,:,mm)
		  s0m = s0(:,:,mm)
		  call spcorr_grd_exp2p (nspl1, nspl2, c0m, s0m, grid)
  
		  ! assign values from pointers to matrices
		  sp_wght_var = 0.0
		  sp_ipos_var = 0
		  sp_jpos_var = 0
		  do i = 1, nspl1
			do j = 1, nspl2
			  sp_num_var(i, j)=size(spcorr(i, j)%wght)
			  sp_ipos_var(i, j, 1:sp_num_var(i, j))=spcorr(i, j)%ipos
			  sp_jpos_var(i, j, 1:sp_num_var(i, j))=spcorr(i, j)%jpos
			  sp_wght_var(i, j, 1:sp_num_var(i, j))=spcorr(i, j)%wght
			  sp_sdev_var(i, j) = spcorr(i, j)%sdev
			end do
		  end do
		  
		  ! save structure information to files
		  open(unit=34,file= file_spcc_struct,form='unformatted',iostat=error)
		  if(error .ne. 0) then; print *, 'Error opening station weight file', error; stop; end if
		  write(unit=34,iostat=error) sp_wght_var, sp_ipos_var, sp_jpos_var, sp_num_var, sp_sdev_var
		  if(error .ne. 0) then; print *, 'Error writing station weight file ', error; stop; end if
		  close(unit=34)
		  
	  end if
	  
  end do ! end loop mm
 
	
  clen = 150.0
  call spcorr_grd (nspl1, nspl2, grid)
  
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
  allocate (sp_wght_var(nspl1,nspl2,maxprev), sp_ipos_var(nspl1,nspl2,maxprev), & 
    & sp_jpos_var(nspl1,nspl2,maxprev), sp_num_var(nspl1,nspl2), sp_sdev_var(nspl1,nspl2), stat=ierr)

  
  ! start generation
  print *, 'Generating SCRF for every ensemble member and every month'
  do iens = start_ens, stop_ens
    initflag = 1 ! different ensemble members are independent with each other
    do i = 1979, 2018
      do j = 1, 12
        print *, 'Generating random numbers: member/year/month', iens, i, j
        write( mmstr, '(i2)' )  j

        ! load spatial correlation structure
        print *, 'load spatial correlation structure'
        file_spcc_struct = trim(out_spcorr_prefix) // trim(mmstr)
		open(unit=34,file=file_spcc_struct,form='unformatted',iostat=error)
	    read(unit=34,iostat=error) sp_wght_var, sp_ipos_var, sp_jpos_var, sp_num_var, sp_sdev_var
		close(unit=34)

		! estimate the number of days in the month
		if ((j .eq. 1) .or. (j .eq. 3) .or. (j .eq. 5) .or. (j .eq. 7) .or. (j .eq. 8) .or. (j .eq. 10) .or. (j .eq. 12)) ntimes = 31
		if ((j .eq. 4) .or. (j .eq. 6) .or. (j .eq. 9) .or. (j .eq. 11)) ntimes = 30
		if (j .eq. 2) then
			if (((mod(i,4) .eq. 0) .and. (mod(i,100) .ne. 0)) .or. (mod(i,400) .eq. 0)) then
		   		ntimes = 29
			else
		   		ntimes = 28
		    end if 
		end if  
		
		! generate random numbers for every day
		print *, 'generate random numbers'	
	    if (allocated(rndnum_3D))  deallocate (rndnum_3D, stat=ierr)
		allocate (rndnum_3D(nx, ny, ntimes), stat=ierr)
		rndnum_3D = 0.0
		
		do istep = 1, ntimes   ! loop for all days in one month
		   print *, 'current time step', istep, '    ///    total steps', ntimes 
		   if (initflag .eq. 1) then
		   	  call field_rand_nopointer (nspl1, nspl2, sp_wght_var, sp_sdev_var, sp_ipos_var, sp_jpos_var, sp_num_var, rndnum_2D)
		   	  rndnum_3D(:, :, istep) = rndnum_2D
			  initflag = 0
		   else
              old_random = rndnum_2D
      		  call field_rand_nopointer (nspl1, nspl2, sp_wght_var, sp_sdev_var, sp_ipos_var, sp_jpos_var, sp_num_var, rndnum_2D)
      		  rndnum_3D(:, :, istep) = old_random * auto_corr(j) + sqrt (1-auto_corr(j)*auto_corr(j)) * rndnum_2D
		   end if
		end do ! end loop days

		! save random numbers to netcdf files
		print *, 'save random numbers to nc file'
		write( mmstr, '(i6)' )  i*100+j
		write (suffix, '(I3.3)') iens
		file_scrf = trim(out_forc_name_base) // '/scrf.' // trim(mmstr) // '.' // trim(suffix) // '.nc'
		call save_rndnum (rndnum_3D, file_scrf, ierr)
          
      end do !end month
    end do !end year
  end do !end ensemble member loop
   
 
end program generate_rndnum_exp2p
