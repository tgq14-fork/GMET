program generate_rndnum
! Creator: Guoqiang Tang, 2020
! Purpose: Produce spatiotemporally correlated random numbers (SCRF) for all days and 
! all grids in a single run

! -----------------------------------------------------------------------------
! Creator(s):
!   Andy Newman, 2013
! Modified:
!   Andy Wood, 2016 -- adding namelist file as argument
!                   -- no longer hardwired; clean formatting
!                   -- adding documentation
!                   -- altered namelist args to specify ens range
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
  use namelist_module, only: read_namelist !namelist module
  ! use namelist_module, only: nens, ntimes, start_time AW edited
  use namelist_module, only: start_ens, stop_ens, ntimes, start_time
  use namelist_module, only: out_forc_name_base, in_regr_name, grid_name, clen
  use namelist_module, only: time_mode
  use namelist_module, only: climo_path

  implicit none
 
  interface

 	! Add by TGQ
 	subroutine save_rndnum (pcp_rndnum, tmean_rndnum, trange_rndnum, nx, ny, ntimes, grdlat, grdlon, grdalt, file, error)
      use netcdf
      use nrtype
 
      real (sp), intent (in) :: pcp_rndnum (:, :, :), tmean_rndnum (:, :, :), trange_rndnum (:, :, :)
      integer (i4b), intent (in) :: nx, ny, ntimes
      real (dp), intent (in) :: grdlat (:), grdlon (:), grdalt (:)
      character (len=500), intent (in) :: file
      integer, intent (out) :: error
    end subroutine save_rndnum
 	! Add by TGQ
 	
 
    function erfinv (x)
      use nrtype
 
      real (sp), intent (in) :: x
      real (sp) :: erfinv
    end function erfinv
 
    subroutine read_nc_grid (file_name, lat, lon, elev, grad_n, grad_e, mask, nx, ny, error)
      use netcdf
      use nrtype
 
      character (len=*), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: lat (:, :), lon (:, :), elev (:, :), grad_n (:, :), &
     & grad_e (:, :), mask (:, :)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_nc_grid


  end interface
  ! ================== END of INTERFACES ===============
 
  ! Local variables
  integer (i4b) :: i, j, k, igrd, istep, iens, mm !counter variables
  character(10) ::  mmstr
  integer (i4b) :: nens  ! AW number of ensemble members to generate
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
  real (dp), allocatable :: lat (:, :)
  real (dp), allocatable :: lon (:, :)
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
  real (sp), allocatable :: pcp_out (:, :, :)!
  real (sp), allocatable :: tmean_out (:, :, :)!
  real (sp), allocatable :: trange_out (:, :, :)!
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
  real(DP),dimension(12) :: clen_daily_tmean= (/1153.0, 1322.0, 1237.0, 1018.0, 961.0, 809.0, 600.0, 615.0, 1049.0, 1217.0, 1383.0, 1250.0/)
  real(DP),dimension(12) :: clen_month_tmean= (/1395.0, 1493.0, 1472.0, 1302.0, 1145.0, 862.0, 967.0, 914.0, 991.0, 1207.0, 1643.0, 1348.0/)
  real(DP),dimension(12) :: clen_daily_prcp= (/303.0, 277.0, 236.0, 189.0, 125.0, 76.0, 47.0, 52.0, 126.0, 215.0, 250.0, 281.0/)
  real(DP),dimension(12) :: clen_month_prcp= (/501.0, 384.0, 325.0, 357.0, 349.0, 283.0, 204.0, 180.0, 280.0, 422.0, 453.0, 449.0/)
  real(DP),dimension(12) :: clen_daily_trange= (/200.0, 191.0, 185.0, 229.0, 231.0, 211.0, 128.0, 121.0, 226.0, 430.0, 285.0, 171.0/)
  real(DP),dimension(12) :: clen_month_trange= (/573.0, 537.0, 379.0, 361.0, 368.0, 423.0, 308.0, 198.0, 399.0, 805.0, 696.0, 515.0/)
  real(DP),dimension(12) :: auto_corr_daily= (/0.629, 0.636, 0.648, 0.634, 0.650, 0.620, 0.559, 0.561, 0.633, 0.654, 0.652, 0.643/)
  real(DP),dimension(12) :: tp_corr_daily= (/-0.167, -0.221, -0.246, -0.277, -0.279, -0.261, -0.229, -0.249, -0.287, -0.266, -0.199, -0.151/)
  ! add by TGQ

  type (coords), pointer :: grid !coordinate structure for grid
  type (splnum), dimension (:, :), pointer :: sp_pcp, sp_temp, sp_trange ! structures of spatially correlated random field weights
  
  real (dp), dimension (:, :, :), allocatable :: sp_wght
  real (dp), dimension (:, :), allocatable :: sp_sdev
  integer (i4b), dimension (:,:,:), allocatable :: sp_ipos, sp_jpos
  integer (i4b), dimension (:,:), allocatable ::sp_num
  logical :: file_exists
  integer (i4b) :: initflag
  ! ========== code starts below ==============================
 
  ! AWW: get namelist filename from command line (no longer hardwired)
  f = 0
  do
    call get_command_argument (f, arg)
    if (f .eq. 1) namelist_filename = arg
    if (len_trim(arg) == 0) exit
    f = f + 1
  end do

  ! read namelist in
  call read_namelist (namelist_filename)
 
  ! set output file name from namelist
  out_name = out_forc_name_base
  error = 0
  ierr = 0
  jerr = 0

  ! AW set number of ensembles to generate
  !TGQ: the codes are repeated?
  nens = stop_ens - start_ens + 1
  print*, 'Generating ',nens,' ensembles from ',start_ens,' to ',stop_ens

  !read in netcdf grid file
  call read_nc_grid (grid_name, lat, lon, hgt, slp_n, slp_e, mask, nx, ny, error)
 
  if (error .ne. 0) call exit_scrf (1, 'problem in read_nc_grid ')
 
  allocate (lat_out(nx*ny), lon_out(nx*ny), hgt_out(nx*ny), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating for 1-d output variables')

  
  ! set up a few variables for spcorr structure
  nspl1 = nx
  nspl2 = ny
  spl1_start = 1
  spl2_start = 1
  spl1_count = nx
  spl2_count = ny
 
  ! allocate space for scrfs
  allocate (sp_pcp(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) then
    call exit_scrf (1, 'problem deallocating space for sp_pcp ')
  end if
 
  allocate (sp_temp(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) then
    call exit_scrf (1, 'problem deallocating space for sp_temp ')
  end if
  
  ! TGQ
  allocate (sp_trange(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) then
    call exit_scrf (1, 'problem deallocating space for sp_trange ')
  end if
  ! TGQ
 
  if (allocated(pcp_random)) then
    deallocate (pcp_random, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating space for pcp_random ')
  end if
  allocate (pcp_random(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating space for pcp_random ')
 
  if (allocated(tmean_random)) then
    deallocate (tmean_random, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating space for tmean_random ')
  end if
  allocate (tmean_random(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating space for tmean_random ')
 
  if (allocated(trange_random)) then
    deallocate (trange_random, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating space for trange_random ')
  end if
  allocate (trange_random(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating space for trange_random ')
 
  if (allocated(old_random)) deallocate (old_random, stat=ierr)
  allocate (old_random(nspl1, nspl2), stat=ierr)
 
  nullify (grid)
  allocate (grid, stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating structure grid')
 
  ! --------------------------------------------------------------------------------------
  ! allocate space for spatial arrays in grid structure
  grid%idx%spl1_start = spl1_start
  grid%idx%spl2_start = spl2_start
  grid%idx%spl1_count = spl1_count
  grid%idx%spl2_count = spl2_count
 
  allocate (grid%lat(spl1_count, spl2_count), grid%lon(spl1_count, spl2_count), &
 & grid%elv(spl1_count, spl2_count), stat=jerr)
  if (ierr .ne. 0 .or. jerr .ne. 0) call exit_scrf (1, ' problem allocating space for&
 & lat-lon-elev coordinates ')
 
  allocate (pcp_out(nx, ny, ntimes), tmean_out(nx, ny, ntimes), trange_out(nx, ny, ntimes), &
 & stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating for 2-d output variables')

  lon_out = pack (lon, .true.)
  lat_out = pack (lat, .true.)
  hgt_out = pack (hgt, .true.)
 
  grid%lat = lat
  grid%lon = lon
  grid%elv = hgt
  
  ! --------------------------------------------------------------------------------------
  print *, 'Generating spatial correlation structure for every month'
  do mm = 1, 12
      print *,'Processing month', mm
      write( mmstr, '(i2)' )  mm
      ! prcp
	  clen = clen_daily_prcp(mm)
	  out_name =  trim(out_forc_name_base) // '/' // trim('spcorr_prcp_month_') // trim(mmstr)
	  INQUIRE(FILE=out_name, EXIST=file_exists)
	  if (.NOT. file_exists) then
	      if (allocated(sp_wght)) deallocate (sp_wght)
	      if (allocated(sp_ipos)) deallocate (sp_ipos)
	      if (allocated(sp_jpos)) deallocate (sp_jpos)
	      if (allocated(sp_num)) deallocate (sp_num)
	      if (allocated(sp_sdev)) deallocate (sp_sdev)
		  allocate (sp_wght(nspl1,nspl2,49), sp_ipos(nspl1,nspl2,49), sp_jpos(nspl1,nspl2,49),sp_num(nspl1,nspl2), sp_sdev(nspl1,nspl2), stat=ierr)
		  if (associated(spcorr)) then
    		deallocate (spcorr, stat=ierr)
    		if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for the spatial correlation structure ')
  		  end if
  		  
		  call spcorr_grd (nspl1, nspl2, grid)
  
		  !output spcorr_grd to file
		  sp_wght = 0.0
		  sp_ipos = 0
		  sp_jpos = 0
		  do i = 1, nspl1
			do j = 1, nspl2
			  sp_num(i, j)=size(spcorr(i, j)%wght)
			  sp_ipos(i, j, 1:sp_num(i, j))=spcorr(i, j)%ipos
			  sp_jpos(i, j, 1:sp_num(i, j))=spcorr(i, j)%jpos
			  sp_wght(i, j, 1:sp_num(i, j))=spcorr(i, j)%wght
			  sp_sdev(i, j) = spcorr(i, j)%sdev
			end do
		  end do
		  open(unit=34,file= out_name,form='unformatted',iostat=error)
		  if(error .ne. 0) then; print *, 'Error opening station weight file', error; stop; end if
		  write(unit=34,iostat=error) sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev
		  if(error .ne. 0) then; print *, 'Error writing station weight file ', error; stop; end if
		  close(unit=34)
	  end if
	  
	  ! tmean
	  clen = clen_daily_tmean(mm)
	  out_name = trim(out_forc_name_base) // '/' // trim('spcorr_tmean_month_') // trim(mmstr)
	  INQUIRE(FILE=out_name, EXIST=file_exists)
	  if (.NOT. file_exists) then
	      if (allocated(sp_wght)) deallocate (sp_wght)
	      if (allocated(sp_ipos)) deallocate (sp_ipos)
	      if (allocated(sp_jpos)) deallocate (sp_jpos)
	      if (allocated(sp_num)) deallocate (sp_num)
	      if (allocated(sp_sdev)) deallocate (sp_sdev)
		  allocate (sp_wght(nspl1,nspl2,49), sp_ipos(nspl1,nspl2,49), sp_jpos(nspl1,nspl2,49),sp_num(nspl1,nspl2), sp_sdev(nspl1,nspl2), stat=ierr)
		  if (associated(spcorr)) then
    		deallocate (spcorr, stat=ierr)
    		if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for the spatial correlation structure ')
  		  end if
  		  
		  call spcorr_grd (nspl1, nspl2, grid)
  
		  !output spcorr_grd to file
		  sp_wght = 0.0
		  sp_ipos = 0
		  sp_jpos = 0
		  do i = 1, nspl1
			do j = 1, nspl2
			  sp_num(i, j)=size(spcorr(i, j)%wght)
			  sp_ipos(i, j, 1:sp_num(i, j))=spcorr(i, j)%ipos
			  sp_jpos(i, j, 1:sp_num(i, j))=spcorr(i, j)%jpos
			  sp_wght(i, j, 1:sp_num(i, j))=spcorr(i, j)%wght
			  sp_sdev(i, j) = spcorr(i, j)%sdev
			end do
		  end do
		  open(unit=34,file= out_name,form='unformatted',iostat=error)
		  if(error .ne. 0) then; print *, 'Error opening station weight file', error; stop; end if
		  write(unit=34,iostat=error) sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev
		  if(error .ne. 0) then; print *, 'Error writing station weight file ', error; stop; end if
		  close(unit=34)
	  end if
	  
	  ! trange
	  clen = clen_daily_trange(mm)
	  out_name = trim(out_forc_name_base) // '/' // trim('spcorr_trange_month_') // trim(mmstr)
	  INQUIRE(FILE=out_name, EXIST=file_exists)
	  if (.NOT. file_exists) then
	      if (allocated(sp_wght)) deallocate (sp_wght)
	      if (allocated(sp_ipos)) deallocate (sp_ipos)
	      if (allocated(sp_jpos)) deallocate (sp_jpos)
	      if (allocated(sp_num)) deallocate (sp_num)
	      if (allocated(sp_sdev)) deallocate (sp_sdev)
		  allocate (sp_wght(nspl1,nspl2,49), sp_ipos(nspl1,nspl2,49), sp_jpos(nspl1,nspl2,49),sp_num(nspl1,nspl2), sp_sdev(nspl1,nspl2), stat=ierr)
		  if (associated(spcorr)) then
    		deallocate (spcorr, stat=ierr)
    		if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for the spatial correlation structure ')
  		  end if
  		  
		  call spcorr_grd (nspl1, nspl2, grid)
  
		  !output spcorr_grd to file
		  sp_wght = 0.0
		  sp_ipos = 0
		  sp_jpos = 0
		  do i = 1, nspl1
			do j = 1, nspl2
			  sp_num(i, j)=size(spcorr(i, j)%wght)
			  sp_ipos(i, j, 1:sp_num(i, j))=spcorr(i, j)%ipos
			  sp_jpos(i, j, 1:sp_num(i, j))=spcorr(i, j)%jpos
			  sp_wght(i, j, 1:sp_num(i, j))=spcorr(i, j)%wght
			  sp_sdev(i, j) = spcorr(i, j)%sdev
			end do
		  end do
		  open(unit=34,file= out_name,form='unformatted',iostat=error)
		  if(error .ne. 0) then; print *, 'Error opening station weight file', error; stop; end if
		  write(unit=34,iostat=error) sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev
		  if(error .ne. 0) then; print *, 'Error writing station weight file ', error; stop; end if
		  close(unit=34)
	  end if
  end do ! end loop mm
  
  print *, 'Generating a useless SCRF'
  call spcorr_grd (nspl1, nspl2, grid) ! necessary to produce iorder and jorder
  
  if (allocated(pcp_rndnum))  deallocate (pcp_rndnum, stat=ierr)
		  if (allocated(tmean_rndnum))  deallocate (tmean_rndnum, stat=ierr)
		  if (allocated(trange_rndnum))  deallocate (trange_rndnum, stat=ierr)
		  allocate (pcp_rndnum(nx, ny, ntimes), tmean_rndnum(nx, ny, ntimes), trange_rndnum(nx, ny, ntimes), stat=ierr)
		  pcp_rndnum = 0.0
  		  tmean_rndnum = 0.0
  		  trange_rndnum = 0.0
  		  
  if (allocated(old_random)) deallocate (old_random, stat=ierr)
		  if (allocated(pcp_random)) deallocate (pcp_random, stat=ierr)
		  if (allocated(tmean_random)) deallocate (tmean_random, stat=ierr)
		  if (allocated(trange_random)) deallocate (trange_random, stat=ierr)
      	  allocate (old_random(nspl1, nspl2),pcp_random(nspl1, nspl2),tmean_random(nspl1, nspl2), trange_random(nspl1, nspl2), stat=ierr)
  ! --------------------------------------------------------------------------------------
  print *, 'Generating SCRF for every month'
  do iens = start_ens, stop_ens
    initflag = 1 ! different ensemble members are independent with each other
    do i = 1979, 2018
      do j = 1, 12
        print *, 'Processing ens/year/month', iens, i, j
        write( mmstr, '(i2)' )  j
        ! --------------------------------------------------------------------------------
        ! load spatial correlation structure
        print *, 'load spatial correlation structure'
        if (associated(sp_pcp)) deallocate (sp_pcp, stat=ierr)
        if (associated(sp_temp)) deallocate (sp_temp, stat=ierr)
        if (associated(sp_trange)) deallocate (sp_trange, stat=ierr)
        allocate (sp_pcp(nspl1, nspl2),sp_temp(nspl1, nspl2),sp_trange(nspl1, nspl2), stat=ierr)
        ! prcp
          out_name = trim(out_forc_name_base) // '/' // trim('spcorr_prcp_month_') // trim(mmstr)
          if (allocated(sp_wght)) deallocate (sp_wght)
	      if (allocated(sp_ipos)) deallocate (sp_ipos)
	      if (allocated(sp_jpos)) deallocate (sp_jpos)
	      if (allocated(sp_num)) deallocate (sp_num)
	      if (allocated(sp_sdev)) deallocate (sp_sdev)
		  allocate (sp_wght(nspl1,nspl2,49), sp_ipos(nspl1,nspl2,49), sp_jpos(nspl1,nspl2,49),sp_num(nspl1,nspl2), sp_sdev(nspl1,nspl2), stat=ierr)
		  open(unit=34,file=out_name,form='unformatted',iostat=error)
		  read(unit=34,iostat=error) sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev
		  close(unit=34)

		  if (associated(spcorr)) deallocate (spcorr, stat=ierr)
		  allocate (spcorr(nspl1, nspl2), stat=ierr)
		  do isp1 = 1, nspl1
			do isp2 = 1, nspl2
			  nullify (spcorr(isp1, isp2)%ipos, spcorr(isp1, isp2)%jpos, spcorr(isp1, isp2)%wght)
			end do
		  end do
  
		  do isp1 = 1, nspl1
			do isp2 = 1, nspl2
			  if (associated(spcorr(isp1, isp2)%ipos))  deallocate (spcorr(isp1, isp2)%ipos, stat=ierr)
			  if (associated(spcorr(isp1, isp2)%jpos))  deallocate (spcorr(isp1, isp2)%jpos, stat=ierr)
			  if (associated(spcorr(isp1, isp2)%wght)) deallocate (spcorr(isp1, isp2)%wght, stat=ierr)			  
			  allocate (spcorr(isp1, isp2)%ipos(sp_num(isp1, isp2)), spcorr(isp1, isp2)%jpos(sp_num(isp1, isp2)), stat=ierr)
			  allocate (spcorr(isp1, isp2)%wght(sp_num(isp1, isp2)), stat=ierr)
			  spcorr(isp1, isp2)%wght=sp_wght(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%ipos=sp_ipos(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%jpos=sp_jpos(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%sdev = sp_sdev(isp1, isp2)
			end do
		  end do
		  sp_pcp = spcorr
		! tmean
          out_name = trim(out_forc_name_base) // '/' // trim('spcorr_tmean_month_') // trim(mmstr)
          if (allocated(sp_wght)) deallocate (sp_wght)
	      if (allocated(sp_ipos)) deallocate (sp_ipos)
	      if (allocated(sp_jpos)) deallocate (sp_jpos)
	      if (allocated(sp_num)) deallocate (sp_num)
	      if (allocated(sp_sdev)) deallocate (sp_sdev)
		  allocate (sp_wght(nspl1,nspl2,49), sp_ipos(nspl1,nspl2,49), sp_jpos(nspl1,nspl2,49),sp_num(nspl1,nspl2), sp_sdev(nspl1,nspl2), stat=ierr)
		  open(unit=34,file=out_name,form='unformatted',iostat=error)
		  read(unit=34,iostat=error) sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev
		  close(unit=34)

		  if (associated(spcorr)) deallocate (spcorr, stat=ierr)
		  allocate (spcorr(nspl1, nspl2), stat=ierr)
		  do isp1 = 1, nspl1
			do isp2 = 1, nspl2
			  nullify (spcorr(isp1, isp2)%ipos, spcorr(isp1, isp2)%jpos, spcorr(isp1, isp2)%wght)
			end do
		  end do
  
		  do isp1 = 1, nspl1
			do isp2 = 1, nspl2
			  if (associated(spcorr(isp1, isp2)%ipos))  deallocate (spcorr(isp1, isp2)%ipos, stat=ierr)
			  if (associated(spcorr(isp1, isp2)%jpos))  deallocate (spcorr(isp1, isp2)%jpos, stat=ierr)
			  if (associated(spcorr(isp1, isp2)%wght)) deallocate (spcorr(isp1, isp2)%wght, stat=ierr)			  
			  allocate (spcorr(isp1, isp2)%ipos(sp_num(isp1, isp2)), spcorr(isp1, isp2)%jpos(sp_num(isp1, isp2)), stat=ierr)
			  allocate (spcorr(isp1, isp2)%wght(sp_num(isp1, isp2)), stat=ierr)
			  spcorr(isp1, isp2)%wght=sp_wght(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%ipos=sp_ipos(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%jpos=sp_jpos(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%sdev = sp_sdev(isp1, isp2)
			end do
		  end do
		  sp_temp = spcorr
		! trange
          out_name = trim(out_forc_name_base) // '/' // trim('spcorr_trange_month_') // trim(mmstr)
          if (allocated(sp_wght)) deallocate (sp_wght)
	      if (allocated(sp_ipos)) deallocate (sp_ipos)
	      if (allocated(sp_jpos)) deallocate (sp_jpos)
	      if (allocated(sp_num)) deallocate (sp_num)
	      if (allocated(sp_sdev)) deallocate (sp_sdev)
		  allocate (sp_wght(nspl1,nspl2,49), sp_ipos(nspl1,nspl2,49), sp_jpos(nspl1,nspl2,49),sp_num(nspl1,nspl2), sp_sdev(nspl1,nspl2), stat=ierr)
		  open(unit=34,file=out_name,form='unformatted',iostat=error)
		  read(unit=34,iostat=error) sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev
		  close(unit=34)

		  if (associated(spcorr)) deallocate (spcorr, stat=ierr)
		  allocate (spcorr(nspl1, nspl2), stat=ierr)
		  do isp1 = 1, nspl1
			do isp2 = 1, nspl2
			  nullify (spcorr(isp1, isp2)%ipos, spcorr(isp1, isp2)%jpos, spcorr(isp1, isp2)%wght)
			end do
		  end do
  
		  do isp1 = 1, nspl1
			do isp2 = 1, nspl2
			  if (associated(spcorr(isp1, isp2)%ipos))  deallocate (spcorr(isp1, isp2)%ipos, stat=ierr)
			  if (associated(spcorr(isp1, isp2)%jpos))  deallocate (spcorr(isp1, isp2)%jpos, stat=ierr)
			  if (associated(spcorr(isp1, isp2)%wght)) deallocate (spcorr(isp1, isp2)%wght, stat=ierr)			  
			  allocate (spcorr(isp1, isp2)%ipos(sp_num(isp1, isp2)), spcorr(isp1, isp2)%jpos(sp_num(isp1, isp2)), stat=ierr)
			  allocate (spcorr(isp1, isp2)%wght(sp_num(isp1, isp2)), stat=ierr)
			  spcorr(isp1, isp2)%wght=sp_wght(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%ipos=sp_ipos(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%jpos=sp_jpos(isp1, isp2,1:sp_num(isp1, isp2))
			  spcorr(isp1, isp2)%sdev = sp_sdev(isp1, isp2)
			end do
		  end do
		  sp_trange = spcorr
		  deallocate (sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev)
! 		  deallocate (spcorr)
		  ! --------------------------------------------------------------------------------
		  ! days of the month
		  if ((j .eq. 1) .or. (j .eq. 3) .or. (j .eq. 5) .or. (j .eq. 7) .or. (j .eq. 8) .or. (j .eq. 10) .or. (j .eq. 12)) then
		    ntimes = 31
		  end if 
		  if ((j .eq. 4) .or. (j .eq. 6) .or. (j .eq. 9) .or. (j .eq. 11)) then
		    ntimes = 30
		  end if 
		  if (j .eq. 2) then
		    if (((mod(j,4) .eq. 0) .and. (mod(j,100) .ne. 0)) .or. (mod(j,400) .eq. 0)) then
		      ntimes = 29
		    else
		      ntimes = 28
		    end if 
		  end if  
		  print *, 'day number', ntimes
		  ! if (allocated(pcp_rndnum))  deallocate (pcp_rndnum, stat=ierr)
! 		  if (allocated(tmean_rndnum))  deallocate (tmean_rndnum, stat=ierr)
! 		  if (allocated(trange_rndnum))  deallocate (trange_rndnum, stat=ierr)
! 		  allocate (pcp_rndnum(nx, ny, ntimes), tmean_rndnum(nx, ny, ntimes), trange_rndnum(nx, ny, ntimes), stat=ierr)
! 		  pcp_rndnum = 0.0
!   		  tmean_rndnum = 0.0
!   		  trange_rndnum = 0.0
		  ! --------------------------------------------------------------------------------
		  print *, 'generate random numbers'
		  ! if (allocated(old_random)) deallocate (old_random, stat=ierr)
! 		  if (allocated(pcp_random)) deallocate (pcp_random, stat=ierr)
! 		  if (allocated(tmean_random)) deallocate (tmean_random, stat=ierr)
! 		  if (allocated(trange_random)) deallocate (trange_random, stat=ierr)
!       	  allocate (old_random(nspl1, nspl2),pcp_random(nspl1, nspl2),tmean_random(nspl1, nspl2), trange_random(nspl1, nspl2), stat=ierr)
		  
		  do istep = 1, ntimes
			  if (initflag .eq. 1) then
			    print *, 'init time step', istep
			  	spcorr = sp_pcp
			    call field_rand (nspl1, nspl2, pcp_random)
			    pcp_rndnum(:, :, istep) = pcp_random
			    
			    spcorr = sp_temp
			    call field_rand (nspl1, nspl2, tmean_random)
			    tmean_rndnum(:, :, istep) = tmean_random
			    
			    spcorr = sp_trange
			    call field_rand (nspl1, nspl2, trange_random)
				trange_rndnum(:, :, istep) = trange_random
				initflag = 0
			  else
			    print *, 'time step', istep
			    
			    spcorr = sp_temp
                old_random = tmean_random
      			call field_rand (nspl1, nspl2, tmean_random)
      			tmean_rndnum(:, :, istep) = old_random * auto_corr_daily(j) + sqrt (1-auto_corr_daily(j)*auto_corr_daily(j)) * tmean_random
      			
      			spcorr = sp_trange
      			old_random = trange_random
      			call field_rand (nspl1, nspl2, trange_random)
      			trange_rndnum(:, :, istep) = old_random * auto_corr_daily(j) + sqrt (1-auto_corr_daily(j)*auto_corr_daily(j)) * trange_random
 
      			spcorr = sp_pcp
      			call field_rand (nspl1, nspl2, pcp_random)
      			pcp_rndnum(:, :, istep) = trange_random * tp_corr_daily (j) + sqrt (1-tp_corr_daily (j)*tp_corr_daily (j)) * pcp_random
			  end if
		  end do ! end days
		  ! --------------------------------------------------------------------------------
		  ! save random numbers to netcdf files
		  print *, 'save random numbers to nc file'
		  write( mmstr, '(i6)' )  i*100+j
		  write (suffix, '(I3.3)') iens
		  out_name = trim(out_forc_name_base) // '/scrf.' // trim(mmstr) // '.' // trim(suffix) // '.nc'
! 		  call save_rndnum (pcp_rndnum, tmean_rndnum, trange_rndnum, nx, ny, ntimes, lat_out, lon_out, hgt_out, out_name, ierr)
          deallocate (pcp_rndnum, tmean_rndnum, trange_rndnum)
      end do !end month
    end do !end year
  end do !end ensemble member loop
   
 
end program generate_rndnum
