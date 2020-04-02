! AWW-2016Jan, modifications to handle time subsetting and reduce mem alloc, and clean up
!   renamed from estimate_precip; add also 'directory' var, changed some var names

subroutine estimate_forcing_regression (gen_sta_weights, sta_weight_name, x, z, ngrid, maxdistance, times, st_rec, end_rec, &
  & stnid, stnvar, directory, pcp, pop, pcperr, tmean, tmean_err, &
  & trange, trange_err, mean_autocorr, mean_tp_corr, y_max, error)

  ! ==============================================================================================
  ! This routine is called during MODE 2 usage:  creates gridded ensembles from station/point data
  ! ==============================================================================================

  use string_mod
  use utim
  use type
  implicit none

  ! ===== start interfaces =======
  interface
    subroutine read_transform_exp (ntimes, file_name, texp)
      use type
      integer (i4b), intent (in) :: ntimes
      character (len=*), intent (in) :: file_name
      real (dp), allocatable, intent (out) :: texp (:)
    end subroutine read_transform_exp

    subroutine read_station (stnvar, stnid, directory, st_rec, end_rec, vals, tair_vals, &
   & vals_miss, vals_miss_t, error)
      use type
      character (len=100), intent (in) :: stnvar
      character (len=100), intent (in) :: stnid
      character (len=500), intent (in) :: directory
      integer (i4b), intent (in) :: st_rec, end_rec
      real (dp), allocatable, intent (out) :: vals (:), tair_vals (:, :)
      logical, allocatable, intent (out) :: vals_miss (:), vals_miss_t (:)
      integer, intent (out) :: error
    end subroutine read_station

    subroutine compute_station_weights(sta_weight_name,ngrid,nstns,X,Z,search_distance, &
                                   sta_limit,sta_data,tair_data, &
                                   close_meta,close_meta_t,close_loc,close_loc_t, &
                                   close_count,close_count_t,close_weights,close_weights_t,error)
      use type
      !inputs
      character(len=500),intent(in) :: sta_weight_name   !name of station weight binary file
      integer(I4B), intent(in)      :: ngrid               !number of grid points
      integer(I4B), intent(in)      :: nstns               !number of stations
      real(DP), intent(in)          :: Z(:,:)              !grid metadata array
      real(DP), intent(in)          :: X(:,:)              !station metadata array
      real(DP), intent(in)          :: search_distance     !default station search distance
      integer(I4B), intent(in)      :: sta_limit           !maximum number of stations for a grid point
      real(DP), intent(in)          :: sta_data(:,:)       !station data values for precipitation
      real(DP), intent(in)          :: tair_data(:,:,:)    !station air temperature data
      !in/out
      real(DP), intent(inout)     :: close_meta(:,:,:)
      real(DP), intent(inout)     :: close_meta_t(:,:,:)
      integer(I4B), intent(inout) :: close_loc(:,:)
      integer(I4B), intent(inout) :: close_loc_t(:,:)
      integer(I4B), intent(inout) :: close_count(:)
      integer(I4B), intent(inout) :: close_count_t(:)
      real(DP), intent(inout)     :: close_weights(:,:)
      real(DP), intent(inout)     :: close_weights_t(:,:) 
      integer(I4B),intent(inout)  :: error
    end subroutine compute_station_weights

    subroutine write_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !input
                                close_weights_t,close_count,close_count_t,error) !input
      use type
      !inputs
      character(len=500), intent(in)    :: sta_weight_name
      real(DP), intent(in)      :: close_meta(:,:,:)
      real(DP), intent(in)      :: close_meta_t(:,:,:)
      integer(I4B), intent(in)  :: close_loc(:,:)
      integer(I4B), intent(in)  :: close_loc_t(:,:)
      real(DP), intent(in)      :: close_weights(:,:)
      real(DP), intent(in)      :: close_weights_t(:,:)
      integer(I4B), intent(in)  :: close_count(:)
      integer(I4B), intent(in)  :: close_count_t(:)
      integer(I4B), intent(inout):: error
    end subroutine write_station_weights

    subroutine read_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !output
                                close_weights_t,close_count,close_count_t,error) !output
      use type
      !input
      character(len=500), intent(in)    :: sta_weight_name
      !output
      real(DP), intent(out)      :: close_meta(:,:,:)
      real(DP), intent(out)      :: close_meta_t(:,:,:)
      integer(I4B), intent(out)  :: close_loc(:,:)
      integer(I4B), intent(out)  :: close_loc_t(:,:)
      real(DP), intent(out)      :: close_weights(:,:)
      real(DP), intent(out)      :: close_weights_t(:,:)
      integer(I4B), intent(out)  :: close_count(:)
      integer(I4B), intent(out)  :: close_count_t(:)
      integer(I4B), intent(inout):: error
    end subroutine read_station_weights

    subroutine normalize_x (x)
      use type
      real (dp), intent (inout) :: x (:, :)
    end subroutine normalize_x

    ! added AJN Sept 2013
    subroutine normalize_xv (x, weight, yp, smax)
      use type
      real (dp), intent (in) :: x (:)
      integer (i4b), intent (in) :: yp (:)
      real (dp), intent (in) :: weight (:)
      real (dp), intent (out) :: smax
    end subroutine normalize_xv

    subroutine normalize_y (texp, y)
      use type
      real (dp), intent (in) :: texp !transform exponent
      real (dp), intent (inout) :: y (:)
    end subroutine normalize_y

    subroutine least_squares (x, y, tx, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine least_squares

    subroutine logistic_regression (x, y, tx, yp, b)
      use type
      real (dp), intent (in) :: x (:, :)
      real (dp), intent (in) :: y (:)
      real (dp), intent (in) :: tx (:, :)
      integer (i4b), intent (in) :: yp (:)
      real (dp), allocatable, intent (out) :: b (:)
    end subroutine logistic_regression

    ! added AJN Sept 2013
    subroutine generic_corr (prcp_data, tair_data, lag, window, auto_corr, t_p_corr)
      use type
      real (dp), intent (in) :: prcp_data (:)
      real (dp), intent (in) :: tair_data (:, :)
      integer (i4b), intent (in) :: lag
      integer (i4b), intent (in) :: window
      real (dp), intent (out) :: auto_corr
      real (dp), intent (out) :: t_p_corr
    end subroutine generic_corr

    subroutine calc_distance (lat1, lon1, lat2, lon2, dist)
      use type
      implicit none
      real (dp), intent (in) :: lat1, lon1, lat2, lon2
      real (dp), intent (out) :: dist
    end subroutine calc_distance

    subroutine heapsort (n, ra, rn)
      use type
      implicit none
      integer (i4b), intent (in) :: n
      integer (i4b), dimension (:), intent (inout) :: rn
      real (dp), dimension (:), intent (inout) :: ra
    end subroutine heapsort

  end interface
  ! =========== end interfaces, start code =============

  real (dp), intent (in) :: x (:, :), z (:, :)  ! station and grid point description arrays
  real (dp), intent (in) :: maxdistance         ! max distance for weight function. !TGQ: This is not used? 
  integer (i4b), intent (in) :: ngrid           ! number of grid points
  real (dp), intent (in) :: times (:)!time step array

  ! AWW added next
  integer (i4b), intent (in) :: st_rec, end_rec

  character (len=100), intent (in) :: stnid (:)!station id array
  character (len=100), intent (in) :: stnvar !control file variables
  character (len=500), intent (in) :: directory
  character (len = 500),intent(in) :: gen_sta_weights     ! flag for generating station weight file
  character (len = 500),intent(in) :: sta_weight_name     ! station weight file name

  real (sp), allocatable, intent (out) :: pcp (:, :), pop (:, :), pcperr (:, :)!output variables for precipitation
  real (sp), allocatable, intent (out) :: tmean (:, :), tmean_err (:, :)!OLS tmean estimate and error
  real (sp), allocatable, intent (out) :: trange (:, :), trange_err (:, :)!OLS trange estimate and error

  integer, intent (out) :: error ! integer error flag
  real (dp), intent (out) :: mean_autocorr (:)!mean autocorrelation from all stations over entire time period
  real (dp), intent (out) :: mean_tp_corr (:)!mean correlation for mean temp and precip

  ! vary at each grid point and time step
  real (dp), intent (out) :: y_max (:, :) !max time step precipitation

  real (dp), allocatable :: y (:), b (:)

  real (dp), allocatable :: y_red (:)! reduced matrix for ...
  real (dp), allocatable :: x_red (:, :)! reduced matrix for ...
  real (dp), allocatable :: x_red_t (:, :)! reduced matrix for ...

  !reduced size matrices for quicker matrix operations
  real (dp), allocatable :: twx_red (:, :), tx_red (:, :)!reduced matricies

  real (dp), allocatable :: w_base (:, :)!initial distance weight matrix
  real (dp), allocatable :: w_pcp_1d (:), w_temp_1d (:)
  integer (i4b), allocatable :: w_pcp_1d_loc (:), w_temp_1d_loc (:)
  !real(DP), allocatable :: w_pcp(:,:), w_temp(:,:) !distance weight matrices for a specific grid point
  real (dp), allocatable :: w_pcp_red (:, :), w_temp_red (:, :)!reduced distance weigth matricies

  real (dp), allocatable :: y_tmean (:), y_trange (:)!transformed station data arrays
  real (dp), allocatable :: y_tmean_red (:), y_trange_red (:)!transformed station data arrays
  real (dp), allocatable :: stn_prcp (:), prcp_data (:, :), tair_data (:, :, :), stn_tair (:, :)! orig stn data arrays
  real (dp), allocatable :: auto_corr (:)   ! lag-1 autocorrelation for stations over entire time period used
  real (dp), allocatable :: t_p_corr (:)    ! correlation between temp and precip
  integer (i4b), allocatable :: yp (:)      ! binary for logistic regression
  integer (i4b), allocatable :: yp_red (:)  ! reduced binary for logistic regression

  real(DP), allocatable :: y_red_loocv(:)                                                 !size reduced obs precip matrix for loocv error
  real(DP), allocatable :: twx_red_loocv(:,:), tx_red_loocv(:,:), x_red_loocv(:,:)        !size reduced station matrices for loocv error estimation
  real(DP), allocatable :: y_tmean_red_loocv(:), y_trange_red_loocv(:)                    !size reduced precip matrices for loocv error estimation
  real(DP), allocatable :: x_red_t_loocv(:,:)                                             !size reduced temp matrices for loocv error estimation
  real(DP), allocatable :: w_pcp_red_loocv(:,:), w_temp_red_loocv(:,:)                    !size reduced weight matrices for loocv error estimation

  real(SP)    :: pcp_tmp    !loocv precipitation prediction
  real(SP)    :: temp_pred  !loocv temperature prediction

  logical, allocatable :: stn_miss (:), stn_miss_t (:) ! missing value logical arrays

  real (dp) :: errsum, wgtsum, sta_temp
  real (dp) :: auto_corr_sum, tp_corr_sum
  real (dp) :: step_max ! timestep statistics

  integer (i4b) :: xsize !size of second dimension of input X array
  integer (i4b) :: ntimes, nstns
  integer (i4b) :: t, i, j, g, ndata, nodata
  integer (i4b) :: ndata_t, nodata_t
  integer (i4b) :: lag, window
  integer (i4b) :: auto_cnt, tp_cnt

  ! variables for tracking closest N stations for precipitation
  integer (i4b), parameter :: sta_limit = 30
  integer (i4b), allocatable :: close_loc (:, :)
  integer (i4b), allocatable :: close_count (:)
  real (dp), allocatable :: close_weights (:, :)
  real (dp), allocatable :: close_meta (:, :, :)
  real (dp) :: max_distance
  real (dp), parameter :: search_distance = 1000.0

  ! variables for tracking closest N stations for temperature
  integer (i4b), allocatable :: close_loc_t (:, :)
  integer (i4b), allocatable :: close_count_t (:)
  real (dp), allocatable :: close_weights_t (:, :)
  real (dp), allocatable :: close_meta_t (:, :, :)
  real (dp) :: max_distance_t
  real (dp) :: tmp_weight
  real(DP),allocatable  :: tmp_weight_arr(:,:)

  integer (i4b) :: slope_flag_pcp
  integer (i4b) :: slope_flag_temp

  ! variables to check for singular matrix
  real (dp), allocatable :: mat_test (:, :)
  real (dp), allocatable :: tmp (:, :)
  real (dp), allocatable :: vv (:)

  ! variables for timing code AJN
  integer (i4b) :: t1, t2, count_rate
  integer (i4b) :: tg1, tg2

  !==============================================================!
  !                     code starts below here                   !
  !==============================================================!

  nstns = size (stnid)
  ntimes = size (times)
  xsize = size (x, 2)

  ! allocate variables
  allocate (y(nstns))
  allocate (y_tmean(nstns), y_trange(nstns))
  allocate (prcp_data(nstns, ntimes))
  allocate (tair_data(2, nstns, ntimes))
  allocate (w_pcp_red(sta_limit, sta_limit))
  allocate (w_temp_red(sta_limit, sta_limit))
  allocate (y_red(sta_limit))
  allocate (y_tmean_red(sta_limit), y_trange_red(sta_limit))
  allocate (x_red(sta_limit, xsize))
  allocate (x_red_t(sta_limit, xsize))
  allocate (w_pcp_1d(sta_limit))
  allocate (w_temp_1d(sta_limit))
  allocate (w_pcp_1d_loc(sta_limit))
  allocate (w_temp_1d_loc(sta_limit))
  allocate (tmp(6, 6))
  allocate (vv(6))
  allocate (pcp(ngrid, ntimes))
  allocate (pop(ngrid, ntimes))
  allocate (pcperr(ngrid, ntimes))
  allocate (tmean(ngrid, ntimes))
  allocate (tmean_err(ngrid, ntimes))
  allocate (trange(ngrid, ntimes))
  allocate (trange_err(ngrid, ntimes))
  allocate (auto_corr(nstns))
  allocate (t_p_corr(nstns))
  allocate (yp(nstns))
  allocate (yp_red(sta_limit))
  allocate (mat_test(6,sta_limit))

  ! station limit arrays (precip)
  allocate (close_weights(ngrid, sta_limit))
  allocate (close_loc(ngrid, sta_limit))
  allocate (close_meta(5, ngrid, sta_limit))
  allocate (close_count(ngrid))

  ! station limit arrays (temp)
  allocate (close_weights_t(ngrid, sta_limit))
  allocate (close_loc_t(ngrid, sta_limit))
  allocate (close_meta_t(5, ngrid, sta_limit))
  allocate (close_count_t(ngrid))

  ! base weight array
  allocate (w_base(ngrid, nstns))

 !loocv error estimation vectors and matrices
  allocate(y_red_loocv(sta_limit-1))
  allocate(x_red_loocv(sta_limit-1,xsize))
  allocate(w_pcp_red_loocv(sta_limit-1,sta_limit-1))
  allocate(x_red_t_loocv(sta_limit-1,xsize))
  allocate(w_temp_red_loocv(sta_limit-1,sta_limit-1))
  allocate(Y_tmean_red_loocv(sta_limit-1),Y_trange_red_loocv(sta_limit-1)) ! revise by TGQ
  allocate(tmp_weight_arr(sta_limit,sta_limit))


  ! initializations
  pcp = 0.0d0
  pop = 0.0d0
  pcperr = 0.0d0
  auto_corr = 0.0d0

  tmean = 0.0d0
  trange = 0.0d0
  tmean_err = 0.0d0
  trange_err = 0.0d0

  auto_corr_sum = 0.0d0
  auto_cnt = 0
  tp_corr_sum = 0.0d0
  tp_cnt = 0

  w_base = 0.0d0

  ! ================= LOOP OVER STATIONS ============
  ! this part calls subroutines that calc various  correlations
  ! can do autocorrelations and correlation between temperature and precipitation
  ! uses an n-day moving average (window) to remove "monthly" cycle from temp
  ! and computes autocorrelation on the anomalies

  do i = 1, nstns, 1

    call read_station (stnvar, stnid(i), directory, st_rec, end_rec, stn_prcp, stn_tair, &
   & stn_miss, stn_miss_t, error)
	
	! tgq: those output may be unnecessary because when there are too many stations, it is impossible to identify what is being output in screen
!     print*, "first value check: "
!     print*, "   prcp(1): ",stn_prcp(1)
!     print*, "   tavg(1): ",stn_tair(1,1)
!     print*, "   trng(1): ",stn_tair(2,1)
!     print*

    prcp_data (i, :) = stn_prcp
    tair_data (1, i, :) = stn_tair (1, :)
    tair_data (2, i, :) = stn_tair (2, :)

    ! check data if needed
    ! print *, '-----------'
    ! print *, 'precip',prcp_data(i,:),'MISS',stn_miss
    ! print *, '-----------'
    ! print *,'tmin',tair_data(1,i,:),'MISS',stn_miss_t
    ! print *, '-----------'

    ! compute mean autocorrelation for all stations and all times
    lag = 1
    window = 31 ! AWW:  hardwired parameter; should bring out
    call generic_corr (prcp_data(i, :), tair_data(:, i, :), lag, window, auto_corr(i), t_p_corr(i))
    ! print *,auto_corr(i)

    ! check for values outside of -1 to 1
    ! stations with incomplete data are set to -999
    print *, 'station id & auto_corr'
    print *, stnid(i), auto_corr(i)
    if (auto_corr(i) .ge.-1.0 .and. auto_corr(i) .le. 1.0) then
      auto_corr_sum = auto_corr_sum + auto_corr (i)
      auto_cnt = auto_cnt + 1
    end if
    if (t_p_corr(i) .ge.-1.0 .and. t_p_corr(i) .le. 1.0) then
      tp_corr_sum = tp_corr_sum + t_p_corr (i)
      tp_cnt = tp_cnt + 1
    end if

    deallocate (stn_miss_t)  ! must be allocated within read_station
    deallocate (stn_miss)
    deallocate (stn_prcp)
    deallocate (stn_tair)

  end do
  ! =========== end station read loop ============

  error = 0 ! AWW:  why is this set?  not used again in subroutine

  ! AWW: some checks
  print *, 'auto_cnt, tp_cnt=', auto_cnt, tp_cnt
  if (auto_cnt == 0 .or. tp_cnt == 0) then
    print *, 'ERROR:  autocorr or crosscorr (TxP) could not be calculated due to lack of matching p&
   &airs'
    stop
  end if
  mean_autocorr = auto_corr_sum / real (auto_cnt, kind(dp))
  mean_tp_corr = tp_corr_sum / real (tp_cnt, kind(dp))

  print *, ' '
  print *, '===================================================='
  print *, 'Temp lag-1 autocorrelation: ', mean_autocorr (1)
  print *, 'Temp-precip correlation: ', mean_tp_corr (1)
  print *, '===================================================='
  print *, ' '

  call system_clock (t1, count_rate)

  ! ========= LOOP OVER GRID CELLS ==================
  !Create station-grid cell weight matrices before time stepping
  print *, 'Generating base weight matrix and finding nearest stations for each gridpoint'

  if(gen_sta_weights .eq. "TRUE" .or. gen_sta_weights .eq. "true") then
    call compute_station_weights(sta_weight_name,ngrid,nstns,X,Z,search_distance, & !input
                                 sta_limit,prcp_data,tair_data, & !input
                                 close_meta,close_meta_t,close_loc,close_loc_t, &  !output
                                 close_count,close_count_t,close_weights,close_weights_t,error) !output

    if(error /= 0) then
       return
    endif
    call write_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !input
                                close_weights_t,close_count,close_count_t,error) !input
    if(error /= 0) then
       return
    endif
  else
    call read_station_weights(sta_weight_name, & !input
                                close_meta,close_meta_t,close_loc,close_loc_t,close_weights,& !output
                                close_weights_t,close_count,close_count_t,error) !output
    if(error /= 0) then
       return
    endif
  endif

  call system_clock (t2, count_rate)
  print *, 'Elapsed time for weight generation: ', real (t2-t1) / real (count_rate)

  !for LOOCV, keep constant weights for stations
  do j = 1,sta_limit
    tmp_weight_arr(j,j) = 1.0
  enddo

  ! =========== now LOOP through all TIME steps and populate grids ===============
  do t = 1, ntimes, 1

    call system_clock (tg1, count_rate)
    print *, "TIME STEP = ", times (t), " (", t, "/", ntimes, ")"

    ! --- assign vectors of station values for prcp, temp, for current time step
    do i = 1, nstns, 1
      y (i) = prcp_data (i, t)
      y_tmean (i) = (tair_data(1, i, t)+tair_data(2, i, t)) / 2.0d0
      y_trange (i) = abs (tair_data(2, i, t)-tair_data(1, i, t))
    end do

    ! do power transformation on precip vector (AWW: consider alternate transforms)
    call normalize_y (4.0d0, y)    ! SHOULD NOT BE HARDWIRED

    ! -------- loop through all grid cells for a given time step --------
    do g = 1, ngrid, 1
      ! call system_clock(tg1,count_rate)

      ! IF the elevation is valid for this grid cell
      ! (this starts a long section working first on precip, then temp)
      if (z(g, 4) .gt.-200) then
        ! call system_clock(t1,count_rate)

        ! want to reset weights for closest sta_limit stations...
        ! recalc calc_distance_weight function for selected stations
        ! set max_distance equal to the farthest station distance

        ! ---- first, PRECIP ----

        ! set data count integers and initialize reduced arrays to zero
        ndata = 0
        nodata = 0
        w_pcp_red = 0.0
        y_red = 0.0
        x_red = 0.0
        yp_red = 0

        max_distance = 0.0d0
        do i = 1, (close_count(g)-1) ! close_count is sta_limit+1. "close_count(g)-1" should be done in station_weights.f90
          if (close_meta(5, g, i) .gt. max_distance) then
            max_distance = close_meta (5, g, i)
          end if
        end do

        ! reduced matrices for precip
        slope_flag_pcp = 0
        do i = 1, (close_count(g)-1)
          call calc_distance_weight (max_distance, close_meta(1, g, i), close_meta(2, g, i), &
         & close_meta(3, g, i), close_meta(4, g, i), tmp_weight)

          w_pcp_red (i, i) = tmp_weight ! diagonal matrix: temporal weight based on max_distance
          w_pcp_1d (i) = tmp_weight
          w_pcp_1d_loc (i) = close_loc (g, i) ! close station location/id/number
          y_red (i) = y (close_loc(g, i)) ! close station prcp
          x_red (i, :) = x (close_loc(g, i), :) ! close station basic information [1, lat, lon, altitude, slp_n, slp_e]

          if (prcp_data(close_loc(g, i), t) .gt. 0.0) then
            ndata = ndata + 1    ! count data points with non-zero precip
            yp_red (i) = 1
          else
            nodata = nodata + 1  ! count data points with zero precip
          end if
        end do

        call normalize_xv (y_red, w_pcp_1d, yp_red, step_max) ! TGQ: why calling a complex function normalize_xv just to get the max value of y_red?

        y_max (g, t) = step_max

        ! ---- second, TEMPERATURES ----

        ! start with initializations
        ndata_t = 0
        nodata_t = 0
        w_temp_red = 0.0
        y_tmean_red = 0.0
        y_trange_red = 0.0
        x_red_t = 0.0

        ! max_distance_t = maxval(close_meta_t(5,g,:))
        max_distance_t = 0.0d0
        do i = 1, (close_count_t(g)-1)
          if (close_meta_t(5, g, i) .gt. max_distance_t) then
            max_distance_t = close_meta_t (5, g, i)
          end if
        end do

        ! reduced matrices for temperature
        slope_flag_temp = 0
        do i = 1, (close_count_t(g)-1)
          call calc_distance_weight (max_distance_t, close_meta_t(1, g, i), close_meta_t(2, g, i), &
         & close_meta_t(3, g, i), close_meta_t(4, g, i), tmp_weight)

          w_temp_red (i, i) = tmp_weight
          w_temp_1d (i) = tmp_weight
          y_tmean_red (i) = y_tmean (close_loc_t(g, i))
          y_trange_red (i) = y_trange (close_loc_t(g, i))
          x_red_t (i, :) = x (close_loc_t(g, i), :)

          if (y_tmean(close_loc_t(g, i)) .gt.-100.0) then
            ndata_t = ndata_t + 1    ! count data points with valid temperature
          else
            nodata_t = nodata_t + 1  ! count data points with invalid temperature
          end if
        end do

        ! ---- checks on station availability for precip and temp

        if (ndata == 0 .and. nodata == 0) then
          !print *, "WARNING:  No stations with data within max distance of grid cell!"
          ! this should not happen if station data are filled
          pop (g, t) = 0.0
          pcp (g, t) = 0.0
          pcperr (g, t) = 0.0
        end if

        ! added AJN Sept 2013
        if (ndata_t == 0 .and. nodata_t == 0) then
          if (t .gt. 1) then
            tmean (g, t) = tmean (g, t-1)
            trange (g, t) = trange (g, t-1)
            tmean_err (g, t) = tmean_err (g, t-1)
            trange_err (g, t) = trange_err (g, t-1)
          else
            tmean (g, t) = -999.
            trange (g, t) = -999.
            tmean_err (g, t) = 0.0
            trange_err (g, t) = 0.0
          end if
        end if
        !print *,ndata

        ! ========= Precip & temp are processed sequentially, again ===========
        ! this is the start of the PRECIP processing block ---

        if (ndata >= 1) then  ! at least one station close by has pcp > 0

          ! tmp needs to be matmul(TX,X) where TX = TWX_red and X = X_red
          mat_test = matmul (transpose(x_red), w_pcp_red)
          tmp = matmul (mat_test, x_red)
          vv = maxval (abs(tmp), dim=2)
          
          if (any(vv == 0.0) .or. (abs(Z(g,5)) .lt. 3.6 .and. abs(Z(g,6)) .lt. 3.6)) then
            slope_flag_pcp = 0
          else
            slope_flag_pcp = 1
          end if

          ! -------------- 1. CALCULATING POP -----------------
          if (nodata == 0) then
            print *, "All stations have precip>0, POP = 1.0"
            pop (g, t) = 1.0
          else
            ! some stations don't have precip > 0
            if (slope_flag_pcp .eq. 0) then
              ! --- regression without slope ---
              ! AWW note that these now use the 2nd set of T* variables (different dimension)

              allocate(tx_red(4, sta_limit))
              allocate(twx_red(4, sta_limit))

              tx_red = transpose(x_red(:, 1:4))
              twx_red = matmul(tx_red, w_pcp_red)
              call logistic_regression (x_red(:, 1:4), y_red, twx_red, yp_red, b)!AJN

              if(-dot_product(Z(g,1:4),B) < 25.) then
                pop (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, 1:4), b))), kind(sp))
              else
                POP(g,t) = 0.0
              endif

              deallocate (b)! B must get allocated in logistic reg.; could this also be allocated just once?
              deallocate(tx_red)
              deallocate(twx_red)
            else
              ! --- regression with slope ---
              allocate(tx_red(6, sta_limit))
              allocate(twx_red(6, sta_limit))
              tx_red = transpose (x_red)
              twx_red = matmul (tx_red, w_pcp_red)
              call logistic_regression (x_red, y_red, twx_red, yp_red, b) ! AJN

              !pop (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, :), b))), kind(sp))
              if(-dot_product(Z(g,:),B) < 25.) then
                pop (g, t) = real (1.0/(1.0+exp(-dot_product(z(g, :), b))), kind(sp))
              else
                POP(g,t) = 0.0
              endif

              deallocate (b)
              deallocate(tx_red)
              deallocate(twx_red)
            end if

          end if
          ! print *, "POP: ", POP(g,t)

          ! -------------- 2. NOW CALCULATING PCP -----------------
          if(slope_flag_pcp .eq. 0) then
            !regression without slope terms
            allocate(twx_red(4, sta_limit))
            allocate(tx_red(4, sta_limit))
            allocate(tx_red_loocv(4,sta_limit-1))
            allocate(twx_red_loocv(4,sta_limit-1))

            tx_red = transpose (x_red(:, 1:4))
            twx_red = matmul (tx_red, w_pcp_red)
            call least_squares (x_red(:, 1:4), y_red, twx_red, b)

            pcp (g, t) = real (dot_product(z(g, 1:4), b), kind(sp))

            wgtsum = 0.0
            errsum = 0.0
            do i = 1, (close_count(g)-1), 1
              wgtsum = wgtsum + w_pcp_red (i, i)
              if(i .gt. 1) then
                x_red_loocv(1:i-1,:) = x_red(1:i-1,:)
                x_red_loocv(i:close_count(g)-2,:) = x_red(i+1:close_count(g)-1,:)
                w_pcp_red_loocv(1:i-1,1:i-1) = tmp_weight_arr(1:i-1,1:i-1)
                w_pcp_red_loocv(i:close_count(g)-2,i:close_count(g)-2) = tmp_weight_arr(i+1:close_count(g)-1,i+1:close_count(g)-1)
                y_red_loocv(1:i-1) = y_red(1:i-1)
                y_red_loocv(i:close_count(g)-2) = y_red(i+1:close_count(g)-1)
              else
                x_red_loocv(i:close_count(g)-2,:) = x_red(i+1:close_count(g)-1,:)
                w_pcp_red_loocv(i:close_count(g)-2,i:close_count(g)-2) = tmp_weight_arr(i+1:close_count(g)-1,i+1:close_count(g)-1)
                y_red_loocv(i:close_count(g)-2) = y_red(i+1:close_count(g)-1)
              end if

              TX_red_loocv = transpose(X_red_loocv(:,1:4))
              TWX_red_loocv = matmul(TX_red_loocv,w_pcp_red_loocv)
              call least_squares(X_red_loocv(:,1:4), Y_red_loocv, TWX_red_loocv, B)

              pcp_tmp = real(dot_product(X_red(i,1:4),B),kind(sp))

              errsum = errsum + (w_pcp_red(i,i) * (pcp_tmp - Y_red(i))**2 )
            end do

            pcperr (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(sp))

            deallocate(twx_red)
            deallocate(tx_red)
            deallocate(twx_red_loocv)
            deallocate(tx_red_loocv)
            deallocate(B)
          else
            ! regression with slope terms
            allocate(twx_red(6, sta_limit))
            allocate(tx_red(6, sta_limit))
            allocate(tx_red_loocv(6,sta_limit-1))
            allocate(twx_red_loocv(6,sta_limit-1))
            tx_red = transpose (x_red)
            twx_red = matmul (tx_red, w_pcp_red)

            call least_squares (x_red, y_red, twx_red, b)
            pcp (g, t) = real (dot_product(z(g, :), b), kind(sp))
            deallocate (b)  !AWW-seems to be missing

            wgtsum = 0.0
            errsum = 0.0
            do i = 1, (close_count(g)-1), 1
              wgtsum = wgtsum + w_pcp_red (i, i)
              if(i .gt. 1) then
                x_red_loocv(1:i-1,:) = x_red(1:i-1,:)
                x_red_loocv(i:close_count(g)-2,:) = x_red(i+1:close_count(g)-1,:)
                w_pcp_red_loocv(1:i-1,1:i-1) = tmp_weight_arr(1:i-1,1:i-1)
                w_pcp_red_loocv(i:close_count(g)-2,i:close_count(g)-2) = tmp_weight_arr(i+1:close_count(g)-1,i+1:close_count(g)-1)
                y_red_loocv(1:i-1) = y_red(1:i-1)
                y_red_loocv(i:close_count(g)-2) = y_red(i+1:close_count(g)-1)
              else
                x_red_loocv(i:close_count(g)-2,:) = x_red(i+1:close_count(g)-1,:)
                w_pcp_red_loocv(i:close_count(g)-2,i:close_count(g)-2) = tmp_weight_arr(i+1:close_count(g)-1,i+1:close_count(g)-1)
                y_red_loocv(i:close_count(g)-2) = y_red(i+1:close_count(g)-1)
              end if

              TX_red_loocv = transpose(X_red_loocv)
              TWX_red_loocv = matmul(TX_red_loocv,w_pcp_red_loocv)
              call least_squares(X_red_loocv, Y_red_loocv, TWX_red_loocv, B)

              pcp_tmp = real(dot_product(X_red(i,:),B),kind(sp))

              errsum = errsum + (w_pcp_red(i,i) * (pcp_tmp - Y_red(i))**2 )
            end do

            pcperr (g, t) = real ((errsum/wgtsum)**(1.0/2.0), kind(sp))

            deallocate(twx_red)
            deallocate(tx_red)
            deallocate(twx_red_loocv)
            deallocate(tx_red_loocv)
            deallocate(b)
          end if

        else
          ! this means ndata = 0 for this grid cell and timestep
          ! print *, "INFO:  No stations nearby have pcp > 0, so precip for this cell being set to zero"
          pop (g, t) = 0.0
          pcp (g, t) = 0.0
          pcperr (g, t) = 0.0

        end if ! done with precip if (ndata>=1) block

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !  Temperature OLS
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (ndata_t .ge. 1) then !
          ! regression without slope terms  only for temperature
          allocate(tx_red(4,sta_limit))
          allocate(twx_red(4,sta_limit))
          allocate(tx_red_loocv(4,sta_limit-1))
          allocate(twx_red_loocv(4,sta_limit-1))
		  
          tx_red = transpose (x_red_t(:, 1:4))
          twx_red = matmul (tx_red, w_temp_red)
          call least_squares (x_red_t(:,1:4), y_tmean_red, twx_red, b)

          tmean (g, t) = real (dot_product(z(g, 1:4), b), kind(sp))

          errsum = 0.0
          wgtsum = 0.0
          do i = 1, (close_count_t(g)-1), 1
            wgtsum = wgtsum + w_temp_red (i, i)
            if(i .gt. 1) then
              x_red_t_loocv(1:i-1,:) = x_red_t(1:i-1,:)
              x_red_t_loocv(i:close_count_t(g)-2,:) = x_red_t(i+1:close_count_t(g)-1,:)
              w_temp_red_loocv(1:i-1,1:i-1) = tmp_weight_arr(1:i-1,1:i-1)
              w_temp_red_loocv(i:close_count_t(g)-2,i:close_count_t(g)-2) = tmp_weight_arr(i+1:close_count_t(g)-1,i+1:close_count_t(g)-1)
              y_tmean_red_loocv(1:i-1) = y_tmean_red(1:i-1)
              y_tmean_red_loocv(i:close_count_t(g)-2) = y_tmean_red(i+1:close_count_t(g)-1)
            else
              x_red_t_loocv(i:close_count_t(g)-2,:) = x_red_t(i+1:close_count_t(g)-1,:)
              w_temp_red_loocv(i:close_count_t(g)-2,i:close_count_t(g)-2) = tmp_weight_arr(i+1:close_count_t(g)-1,i+1:close_count_t(g)-1)
              y_tmean_red_loocv(i:close_count_t(g)-2) = y_tmean_red(i+1:close_count_t(g)-1)
            endif
            TX_red_loocv = transpose(X_red_t_loocv(:, 1:4))
            TWX_red_loocv = matmul(TX_red_loocv,w_temp_red_loocv)
            call least_squares(X_red_t_loocv(:, 1:4), Y_tmean_red_loocv,TWX_red_loocv, B)
			
            temp_pred = dot_product(X_red_t(i,1:4),B)
            errsum = errsum + (w_temp_red(i,i) * (temp_pred - Y_tmean_red(i))**2)
          enddo
          tmean_err(g,t) = real((errsum / wgtsum)**(1.0/2.0),kind(sp))

          deallocate (b)

          ! ===== NOW do TRANGE ============
          !regression without slope terms
          TX_red = transpose(X_red_t(:,1:4))
          TWX_red = matmul(TX_red, w_temp_red)

          call least_squares(X_red_t(:,1:4), y_trange_red, TWX_red, B)
          trange(g,t) = real(dot_product(Z(g,1:4),B),kind(sp))

          errsum = 0.0
          wgtsum = 0.0
          !do i = 1, nstns, 1
          do i = 1,(close_count_t(g)-1),1
              wgtsum = wgtsum + w_temp_red(i,i)

              if(i .gt. 1) then
                x_red_t_loocv(1:i-1,:) = x_red_t(1:i-1,:)
                x_red_t_loocv(i:close_count_t(g)-2,:) = x_red_t(i+1:close_count_t(g)-1,:)
                w_temp_red_loocv(1:i-1,1:i-1) = tmp_weight_arr(1:i-1,1:i-1)
                w_temp_red_loocv(i:close_count_t(g)-2,i:close_count_t(g)-2) = tmp_weight_arr(i+1:close_count_t(g)-1,i+1:close_count_t(g)-1)
                y_trange_red_loocv(1:i-1) = y_trange_red(1:i-1)
                y_trange_red_loocv(i:close_count_t(g)-2) = y_trange_red(i+1:close_count_t(g)-1)
              else
                x_red_t_loocv(i:close_count_t(g)-2,:) = x_red_t(i+1:close_count_t(g)-1,:)
                w_temp_red_loocv(i:close_count_t(g)-2,i:close_count_t(g)-2) = tmp_weight_arr(i+1:close_count_t(g)-1,i+1:close_count_t(g)-1)
                y_trange_red_loocv(i:close_count_t(g)-2) = y_trange_red(i+1:close_count_t(g)-1)
              endif

              TX_red_loocv = transpose(X_red_t_loocv(:,1:4))
              TWX_red_loocv = matmul(TX_red_loocv,w_temp_red_loocv)
              call least_squares(X_red_t_loocv(:,1:4), Y_trange_red_loocv,TWX_red_loocv, B)
              temp_pred = dot_product(X_red_t(i,1:4),B)

              errsum = errsum + (w_temp_red(i,i) * (temp_pred - Y_trange_red(i))**2)

          enddo
          trange_err(g,t) = real((errsum / wgtsum)**(1.0/2.0),kind(sp))

          deallocate (b)  !AWW-seems to be missing
          deallocate(tx_red)
          deallocate(twx_red)
          deallocate(tx_red_loocv)
          deallocate(twx_red_loocv)

        else ! alternative to having (ndata_t <= 1)

          ! if not enough stations with data
          ! just use value from previous grid point for now AJN
          print *, 'WARNING:  not enough data stations for current point for temperature; using las&
         &t grid point'

          if (g .gt. 1) then
            trange (g, t) = trange (g-1, t)
            trange_err (g, t) = trange_err (g-1, t)
            tmean (g, t) = tmean (g-1, t)
            tmean_err (g, t) = tmean_err (g-1, t)
          else
            trange (g, t) = trange (g, t-1)
            trange_err (g, t) = trange_err (g-1, t-1)
            tmean (g, t) = tmean (g, t-1)
            tmean_err (g, t) = tmean_err (g, t-1)
          end if
        end if ! end data check if statement for temperature

      end if ! end check for valid elevation

    end do ! end grid loop

    call system_clock (tg2, count_rate)
    print *, 'Elapsed time for one time step: ', real (tg2-tg1) / real (count_rate)

  end do ! end time record loop

  ! AWW -- just deallocate once at end of subroutine
!  deallocate (twx_red)
!  deallocate (tx_red)
!  deallocate (tx_red_loocv)
!  deallocate (twx_red_loocv)

end subroutine estimate_forcing_regression
