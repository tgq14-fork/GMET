program generate_ensembles
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
    subroutine read_grid_list (file_name, lats, lons, alts, slp_n, slp_e, nx, ny, error)
      use nrtype
      character (len=500), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: lats (:), lons (:), alts (:), slp_n (:), slp_e (:)
      integer (i4b), intent (out) :: nx, ny
      integer, intent (out) :: error
    end subroutine read_grid_list
 
    subroutine save_vars (pcp, tmean, trange, nx, ny, grdlat, grdlon, grdalt, times, file, error)
      use netcdf
      use nrtype
 
      real (sp), intent (in) :: pcp (:, :, :), tmean (:, :, :), trange (:, :, :)
      integer (i4b), intent (in) :: nx, ny
      real (dp), intent (in) :: grdlat (:), grdlon (:), grdalt (:)
      real (dp), intent (in) :: times (:)
      character (len=500), intent (in) :: file
      integer, intent (out) :: error
    end subroutine save_vars
 	
 	! Add by TGQ
 	subroutine save_vars_rndnum (pcp, tmean, trange, pcp_rndnum, pcp_cprob, tmean_rndnum, trange_rndnum, &
   & nx, ny, grdlat, grdlon, grdalt, times, file, error)
      use netcdf
      use nrtype
 
      real (sp), intent (in) :: pcp (:, :, :), tmean (:, :, :), trange (:, :, :)
      real (sp), intent (in) :: pcp_rndnum (:, :, :), pcp_cprob(:, :, :), tmean_rndnum (:, :, :), trange_rndnum (:, :, :)
      integer (i4b), intent (in) :: nx, ny
      real (dp), intent (in) :: grdlat (:), grdlon (:), grdalt (:)
      real (dp), intent (in) :: times (:)
      character (len=500), intent (in) :: file
      integer, intent (out) :: error
    end subroutine save_vars_rndnum
 	! Add by TGQ
 	
    subroutine read_grid_qpe_nc_ens (file_name, var_name, var, lats, lons, auto_corr, tp_corr, &
   & times, tot_times, error)
      use netcdf
      use nrtype
 
      character (len=*), intent (in) :: file_name
      character (len=*), intent (in) :: var_name
      real (dp), allocatable, intent (out) :: var (:, :, :)
      real (dp), allocatable, intent (out) :: lats (:, :), lons (:, :)
      real (dp), allocatable, intent (out) :: times (:)
      real (dp), allocatable, intent (out) :: auto_corr (:)
      real (dp), allocatable, intent (out) :: tp_corr (:)
 
      integer, intent (out) :: error
      integer (i4b), intent (out) :: tot_times
    end subroutine read_grid_qpe_nc_ens
 
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

    subroutine read_climo_precip (climo_path,current_month,climo_field,climo_out,error)
      use netcdf
      use nrtype
  
      character (len = *), intent(in)     :: climo_path      !path to climo netcdf files
      integer, intent(in)                 :: current_month   !integer month
      character (len = *), intent(in)     :: climo_field     !name of climo variable to read
      real(SP), intent(out)               :: climo_out(:,:,:)    !climo variable grid
      integer, intent(out)                :: error           !error code
    end subroutine read_climo_precip

    subroutine read_climo_temp(climo_path,current_month,climo_field,climo_out,error)
      use netcdf
      use nrtype

      character (len = *), intent(in)     :: climo_path      !path to climo netcdf files
      integer, intent(in)                 :: current_month   !integer month
      character (len = *), intent(in)     :: climo_field     !name of climo variable to read
      real(SP), intent(out)               :: climo_out(:,:)    !climo variable grid 
      integer, intent(out)                :: error           !error code
    end subroutine read_climo_temp

    subroutine read_climo_uncertainty(uncert_path,current_month,uncert_field,uncert_out,error)
      use netcdf
      use nrtype

      character (len = *), intent(in)     :: uncert_path      !path to uncert netcdf files
      integer, intent(in)                 :: current_month   !integer month
      character (len = *), intent(in)     :: uncert_field     !name of uncert variable to read
      real(SP), intent(out)               :: uncert_out(:,:)    !uncertainty variable grid
      integer, intent(out)                :: error           !error code
    end subroutine read_climo_uncertainty

  end interface
  ! ================== END of INTERFACES ===============
 
  ! Local variables
  integer (i4b) :: i, j, k, igrd, istep, iens, mm !counter variables
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
  ! add by TGQ

  type (coords), pointer :: grid !coordinate structure for grid
  type (splnum), dimension (:, :), pointer :: sp_pcp, sp_temp, sp_trange ! structures of spatially correlated random field weights
  
  real (dp), dimension (:, :, :), allocatable :: sp_wght
  real (dp), dimension (:, :), allocatable :: sp_sdev
  integer (i4b), dimension (:,:,:), allocatable :: sp_ipos, sp_jpos
  integer (i4b), dimension (:,:), allocatable ::sp_num
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
  if(nens <= 0) call exit_scrf (1, 'number of ensembles to generate is 0 or less')
  print*, 'Generating ',nens,' ensembles from ',start_ens,' to ',stop_ens
  nens = stop_ens - start_ens + 1
  if(stop_ens <= start_ens) call exit_scrf (1, 'stop_ens equal or before start_ens')
  print*, 'Generating ',nens,' ensembles from ',start_ens,' to ',stop_ens
 
  !read in netcdf grid file
  call read_nc_grid (grid_name, lat, lon, hgt, slp_n, slp_e, mask, nx, ny, error)
 
  if (error .ne. 0) call exit_scrf (1, 'problem in read_nc_grid ')
 
  allocate (lat_out(nx*ny), lon_out(nx*ny), hgt_out(nx*ny), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating for 1-d output variables')
 

 
  allocate (times(ntimes))

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
  
  print *, 'Generating SCRF spatial correlation structure for every month'
  ! spcc structure of prcp and first random number
  call unix_to_date(times(ntimes/2+1),year,current_month,day,hour,minute,second)
  if(trim(time_mode) .eq. 'daily_anom' .or. trim(time_mode) .eq. 'DAILY_ANOM' .or. trim(time_mode) .eq. 'daily' .or. trim(time_mode) .eq. 'DAILY') then
    clen = clen_daily_prcp(current_month)
  elseif(trim(time_mode) .eq. 'climo' .or. trim(time_mode) .eq. 'CLIMO') then
    clen = clen_month_prcp(current_month)
  end if
  
  allocate (sp_wght(nspl1,nspl2,49), sp_ipos(nspl1,nspl2,49), sp_jpos(nspl1,nspl2,49),sp_num(nspl1,nspl2), sp_sdev(nspl1,nspl2), stat=ierr)
  call spcorr_grd (nspl1, nspl2, grid)
  sp_pcp = spcorr !this is location, weigth, and std of previously generated points. it won't be changed.
  
  print *, 'spcorr value', spcorr(800, 400)%ipos
  print *, 'spcorr value', spcorr(800, 400)%wght
  
  !output spcorr_grd to file
  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating spcorr variables')
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
  
  
  print *, 'start output sp_pcp'
  open(unit=34,file='spcorr_grd',form='unformatted',iostat=error)
  if(error .ne. 0) then; print *, 'Error opening station weight file', error; stop; end if
  write(unit=34,iostat=error) sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev
  if(error .ne. 0) then; print *, 'Error writing station weight file ', error; stop; end if
  close(unit=34)
   
  print *, 'start read sp_pcp'
  open(unit=34,file='spcorr_grd',form='unformatted',iostat=error)
  if(error .ne. 0) then; print *, 'Error opening station weight file ', error; stop; end if
  read(unit=34,iostat=error) sp_wght, sp_ipos, sp_jpos, sp_num, sp_sdev
  if(error .ne. 0) then; print *, 'Error reading station weight file ', error; stop; end if
  close(unit=34)
  
  ! allocate space for the correlation structure (SPCORR has a pointer structure)
  if (associated(spcorr)) then
    deallocate (spcorr, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for the spatial correlation st&
   &ructure ')
  end if
  allocate (spcorr(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, ' problem allocating spatial corr struct ')
    do isp1 = 1, nspl1
      do isp2 = 1, nspl2
        nullify (spcorr(isp1, isp2)%ipos, spcorr(isp1, isp2)%jpos, spcorr(isp1, isp2)%wght)
      end do
    end do
  
  
  print *, 'sp_num size', size(sp_num)
  
  do isp1 = 1, nspl1
    do isp2 = 1, nspl2
    ! allocate space for the (i,j) position of previously generated points
          if (associated(spcorr(isp1, isp2)%ipos)) then
            deallocate (spcorr(isp1, isp2)%ipos, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for ipos, jpos ')
          end if
          if (associated(spcorr(isp1, isp2)%jpos)) then
            deallocate (spcorr(isp1, isp2)%jpos, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for ipos, jpos ')
          end if
          allocate (spcorr(isp1, isp2)%ipos(sp_num(isp1, isp2)), spcorr(isp1, isp2)%jpos(sp_num(isp1, isp2)), stat=ierr)
          if (ierr .ne. 0) call exit_scrf (1, ' problem allocating space for the (i,j) position of &
         &previously generated points ')
    ! allocate space for the weights assigned to previously generated points
          if (associated(spcorr(isp1, isp2)%wght)) then
            deallocate (spcorr(isp1, isp2)%wght, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for weights assigned t&
           &o previously generated points ')
          end if
          allocate (spcorr(isp1, isp2)%wght(sp_num(isp1, isp2)), stat=ierr)
          if (ierr .ne. 0) call exit_scrf (1, ' problem allocating space for weights assigned to pr&
         &eviously generated points ')
    ! populate the structures (-1 excludes the current (i,j) point)
    
      spcorr(isp1, isp2)%wght=sp_wght(isp1, isp2,1:sp_num(isp1, isp2))
      spcorr(isp1, isp2)%ipos=sp_ipos(isp1, isp2,1:sp_num(isp1, isp2))
      spcorr(isp1, isp2)%jpos=sp_jpos(isp1, isp2,1:sp_num(isp1, isp2))
      spcorr(isp1, isp2)%sdev = sp_sdev(isp1, isp2)
    end do
  end do
  
  print *, 'spcorr value', spcorr(800, 400)%ipos
  print *, 'spcorr value', spcorr(800, 400)%wght
  
  
 
end program generate_ensembles
