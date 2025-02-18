! This source code file contains a number of subroutines related to reading files

subroutine read_transform_exp (ntimes, file_name, texp)
  use type
  implicit none
 
  integer (i4b), intent (in) :: ntimes
  character (len=*), intent (in) :: file_name
  real (dp), allocatable, intent (out) :: texp (:)
 
  integer (i4b) :: i
 
  allocate (texp(ntimes))
 
  print *, 'Reading transform file: ', trim (file_name), ' ', ntimes, ' times'
  open (55, file=file_name, status='old')
 
  do i = 1, ntimes
    read (55, "(F3.1)") texp (i)
  end do
 
  close (55)
 
end subroutine read_transform_exp
 
subroutine read_refcst (startdate, enddate, file_var, perturbation, var_name, forecast, v, x, y, t, &
& error)
  use string_mod
  use utim
  use type
  implicit none
 
  interface
    subroutine read_nc_file (file_name, valid_time, var_name, var, lats, lons, error)
      use type
      character (len=*), intent (in) :: file_name
      character (len=*), intent (in) :: var_name
      real (dp), intent (in) :: valid_time
      real (dp), allocatable, intent (out) :: var (:, :, :)
      real (dp), allocatable, intent (out) :: lats (:), lons (:)
      integer, intent (out) :: error
    end subroutine read_nc_file
  end interface
 
  character (len=100), intent (in) :: startdate, enddate, file_var, perturbation
  character (len=*), intent (in) :: var_name
  integer (i4b), intent (in) :: forecast
  real (dp), allocatable, intent (out) :: v (:, :), x (:), y (:)
  real (dp), allocatable, intent (out) :: t (:)
  integer, intent (out) :: error
 
  real (dp), allocatable :: lats (:), lons (:), var (:, :, :)
  real (dp) :: valid_time, utime
  character (len=500) :: file
  character (len=100) :: date
  integer (i4b) :: vshape (3)
  integer (i4b) :: sday, eday, ntimes, ngrids
  integer (i4b) :: sec, min, hour, day, month, year
  integer (i4b) :: i, j
 
  error = 0
 
  call parse_date (startdate, year, month, day, hour, min, sec, error)
  sday = julian_date (day, month, year)
  call parse_date (enddate, year, month, day, hour, min, sec, error)
  eday = julian_date (day, month, year)
  ntimes = eday - sday + 1
 
  print *, startdate, 'sday=', sday, enddate, 'eday=', eday, 'ntimes=', ntimes
 
  allocate (t(ntimes))
 
  i = 1
  utime = date_to_unix (startdate)
  print *, 'utime=', utime
 
  do
    if (utime > date_to_unix(enddate)) exit
    call unix_to_date (utime, year, month, day, hour, min, sec)
 
    write (date, "(I4.4I2.2I2.2)") year, month, day
    if (forecast >= 190) then
      file = trim (date) // "/" // trim (file_var) // "_" // trim (date) // "00_" // trim &
     & (perturbation) // "_t190.nc"
    else
      file = trim (date) // "/" // trim (file_var) // "_" // trim (date) // "00_" // trim &
     & (perturbation) // ".nc"
    end if
 
    valid_time = utime + (forecast*3600.0)
 
    print *, 'utime ', utime, 'valid ', valid_time, year, month, day, file
 
    call read_nc_file (file, valid_time, var_name, var, lats, lons, error)
 
    if (trim(var_name) == 'APCP_ens_mean_surface') then
      var = var ** (1.0/4.0)
      print *, 'normalizing GEFS precip'
    end if
 
    if (error == 1) then
      print *, "Failed reading file ", trim (file)
      print *, 'var ', trim (var_name)
      exit
    else
      print ("(3A,F11.0)"), "Success reading file ", trim(file), " Valid at ", valid_time
    end if
 
    vshape = shape (var)
    if (vshape(1) /= size(lons) .or. vshape(2) /= size(lats) .or. vshape(3) /= 1) then
      print *, "Dimensions from file do not match"
      error = 1
      exit
    end if
 
    if (i == 1) then
      ngrids = vshape (1) * vshape (2)
      allocate (v(ngrids, ntimes))
      allocate (x(vshape(2)))
      allocate (y(vshape(1)))
      x (:) = lats (:)
      do j = 1, vshape (1), 1
        y (j) = mod (lons(j), 360.0) - 360.0 ! generates warning
      end do
        !do j = 1, vshape(1), 1
        !   do k = 1, vshape(2), 1
        !      X(((j-1)*vshape(2)) + k) = lats(k)
        !      Y(((j-1)*vshape(2)) + k) = mod(lons(j),360.0)-360.0
        !   enddo
        !enddo
    else
      if (vshape(2) /= size(x) .or. vshape(1) /= size(y)) then
        print *, "Dimensions from file do not match previous files"
        error = 1
        exit
      end if
    end if
 
    t (i) = valid_time
    v (1:ngrids, i) = reshape (var, (/ ngrids /))
 
    i = i + 1
    utime = utime + 86400
 
  end do
 
end subroutine read_refcst

 
subroutine read_station_list (file_name, id, lat, lon, elev, sslp_n, sslp_e, n_stations)
  use string_mod
  use type
  implicit none
 
  character (len=500), intent (in) :: file_name
  character (len=100), allocatable, intent (out) :: id (:)
  real (dp), allocatable, intent (out) :: lat (:), lon (:), elev (:), sslp_n (:), sslp_e (:)
  integer (i4b), intent (out) :: n_stations
 
  ! local variables
  integer (i4b) s, stat
  character (500) line

  print *, "Reading stations list: ", trim (file_name)

  ! open file and find number of stations (rows-1), then rewind it
  open (11, file=file_name, status='old', iostat=stat)
  if (stat /= 0) then
    print *, "Cannot find file: ", trim(file_name); stop    ! program cannot continue
  end if
  n_stations = -1   ! set to skip header in count
  do 
    read (11, *, end=10) 
    n_stations = n_stations + 1 
  end do
10 rewind(11)    ! rewind file to starting position
  print*, 'Found ', n_stations, ' station records'
  
  ! allocate variables to hold station records
  allocate (id(n_stations))
  allocate (lat(n_stations))
  allocate (lon(n_stations))
  allocate (elev(n_stations))
  allocate (sslp_n(n_stations))
  allocate (sslp_e(n_stations))  
 
  ! now read the command separated records (skip header)
  ! fields are ordered:  station ID, lat, lon, elev, slope_n, slope_e
  read (11, *, iostat=stat) line    ! header
  do s = 1, n_stations
    read (11, *, iostat=stat) id(s), lat(s), lon(s), elev(s), sslp_n(s), sslp_e(s)
    !print*, s,trim(id(s)),lat(s),lon(s)
    if (stat /= 0) then
      print*, 'STOP:  problem reading records from station file'; stop
    end if
  end do
  close (11)
  print *, 'Done reading ', (s-1), ' station info records'; print *, ' '

end subroutine read_station_list
 
 
! moved out of subroutine read_station
subroutine check (status)
  use netcdf
  integer, intent (in) :: status
 
  if (status /= nf90_noerr) then
    print *, trim (nf90_strerror(status))
    stop "Stopped"
  end if
end subroutine check
 
 
! ======== subroutine to read station data in netcdf format ===========
!modified EAC Dec 2015 adapted from AJN version for ASCII, replacing it
!modified AWW Dec 2015, Feb 2016 -- added directory, st_rec, end_rec

subroutine read_station (prcp_varname, stnid, directory, st_rec, end_rec, vals, tair_vals, &
                       & vals_miss, vals_miss_t, error)
  use string_mod
  use utim
  use type
  use netcdf
  implicit none
 
  character (len=100), intent (in) :: prcp_varname !! prcp variable name
  character (len=100), intent (in) :: stnid
  character (len=500), intent (in) :: directory ! AWW-feb2016 data dir separated from site list
 
  integer (i4b), intent (in) :: st_rec, end_rec ! AWW
 
  real, allocatable :: vals1 (:), vals2 (:), valsp (:)
  real (dp), allocatable, intent (out) :: vals (:), tair_vals (:, :)! AWW should call vals prcp_vals
  logical, allocatable, intent (out) :: vals_miss (:), vals_miss_t (:)
  integer, intent (out) :: error
  integer :: ncid ! the file_id of the netcdf file
  integer :: varid ! the integer ID of the variable in the netcdf file
  integer :: recorddimid, nctimes
  character (len=20) :: recorddimname
 
  character (len=500) :: file_name
  !character (len=500) :: directory  AWW
  character (len=100) :: stnvar_t1, stnvar_t2
  integer (i4b) :: t
  logical :: fexist
 
  integer :: status, status_tn, status_tx
 
  !!!! === first read in station precipitation data
  error = 0
 
  !!!! Get directory and file name based on station_list and stnid
 
  file_name = trim (directory) // "/" // trim (stnid) // ".nc"
 
  print *, " "
  print *, " ----- Reading station file: ", trim (file_name)
  !!!! length of file and initialization of vals and vals_miss
 
  !!!! check for file
  inquire (file=file_name, exist=fexist)
 
  if ( .not. fexist) then
    print *, "ERROR: Cannot find file: ", trim (file_name)
    print *, "Setting station precipitation timeseries to missing"
    stop! quit and fix problem -- station should not be missing
 
  else ! file exists -- read it
 
    call check (nf90_open(file_name, nf90_nowrite, ncid))
    call check (nf90_inquire(ncid, unlimiteddimid=recorddimid))
    call check (nf90_inquire_dimension(ncid, recorddimid, name=recorddimname, len=nctimes))
    print *, "netcdf file has ", nctimes, " records (nctimes)"
    ! check to make sure desired period is valid
    if (end_rec .gt. nctimes) then
      print *, 'ERROR:  trying to access station data from record greater than those in station'
      print *, 'rec wanted=', end_rec, 'rec available=', nctimes
      stop
    end if
    print *, "storing ", (end_rec-st_rec+1), " output period records from [st_rec, end_rec]: ",st_rec, end_rec
 
    !!! ----- Get precipitation values if they exist -----
    allocate (valsp(nctimes))                        ! nctimes = all the records in the netcdf file
    allocate (vals(end_rec-st_rec+1))
    allocate (vals_miss(end_rec-st_rec+1))
!    allocate(vals(nctimes))
!    allocate(vals_miss(nctimes))
    status = nf90_inq_varid (ncid, prcp_varname, varid)   
    if (status /= nf90_noerr) then
      vals_miss = .false.
      vals = - 999.0
    else
      call check (nf90_get_var(ncid, varid, valsp))  ! read the precip
      print *, "reading in prcp data"
 
      ! AWW here could make it so that if the station var not specified (PT), don't use

      do t = st_rec, end_rec, 1 ! AWW store only within desired period
        if (valsp(t) == nf90_fill_float) then
          vals_miss = .false.   ! this is the reverse of what the name suggests -AWW need to FIX
          vals = - 999.0
          print *, 'found missing precip at timestep', t, ' -- all values set to void'
          exit
        else
          vals (t-st_rec+1) = dble (valsp(t))
          vals_miss (t-st_rec+1) = .true.
        end if
      end do
 
    end if
    deallocate (valsp)
 
    !!! ----- Get temperature values -----
    stnvar_t1 = 'tmin'
    stnvar_t2 = 'tmax'
!    allocate(tair_vals(2, nctimes))
    allocate (tair_vals(2, end_rec-st_rec+1))       ! AW - get all these allocation / deallocates out of loop; inefficient
    allocate (vals1(nctimes))
    allocate (vals2(nctimes))
!    allocate(vals_miss_t(nctimes))
    allocate (vals_miss_t(end_rec-st_rec+1))
    vals_miss_t = .false.
    tair_vals = - 999.0
 
    ! check for tmin and tmax
    ! AW:  check on both variables at once -- if either is missing, don't use temperatures
    !      from this station
    status_tx = nf90_inq_varid (ncid, stnvar_t2, varid)
    status_tn = nf90_inq_varid (ncid, stnvar_t1, varid)! leaves varid setting for TMIN
    if (status_tn /= nf90_noerr .or. status_tx /= nf90_noerr) then
      ! one or more temp variable is missing
      vals_miss_t = .false.
      tair_vals = - 999.0
 
    else
      ! both variables are supposed to be present -- try to read them
      ! check for tmax
      print *, "reading in tmin"
      call check (nf90_get_var(ncid, varid, vals1))! get tmin data
      print *, "reading in tmax"
      call check (nf90_inq_varid(ncid, stnvar_t2, varid))! reset varid for TMAX
      call check (nf90_get_var(ncid, varid, vals2))! get tmax data
 
      do t = st_rec, end_rec, 1 ! AWW store only within desired period
        if (vals1(t) == nf90_fill_float .or. vals2(t) == nf90_fill_float) then
          tair_vals = - 999.0
          vals_miss_t = .false. ! quit loop on any void and assign whole ts to void
          print *, 'found missing temperature(s) at timestep', t, '-- all values set to void'
          exit
        else
          tair_vals (1, t-st_rec+1) = dble (vals1(t))
          tair_vals (2, t-st_rec+1) = dble (vals2(t))
          vals_miss_t = .true.
        end if
      end do
 
      deallocate (vals1)
      deallocate (vals2)
 
    end if ! end of if 'temperature is present' case
 
    ! close netcdf file
    call check (nf90_close(ncid))
    print *, "closed netcdf file"
 
  end if ! end of 'if station file exists' case
 
end subroutine read_station
 
 
subroutine read_grid_list (file_name, lats, lons, elevs, slp_n, slp_e, nx, ny, error)
  use string_mod
  use type
  implicit none
 
  character (len=500), intent (in) :: file_name
  real (dp), allocatable, intent (out) :: lats (:), lons (:), elevs (:), slp_n (:), slp_e (:)
  integer (i4b), intent (out) :: nx, ny
  integer, intent (out) :: error
 
  character (len=100) :: settings (5)
  integer :: ngrid, i, nsettings, stat, err, ipos
  real (dp) :: dx, dy, startx, starty
  character (300) :: line
  logical fexist
 
  ngrid = 0
  nx = 0
  ny = 0
  dx = 0.0
  dy = 0.0
  startx = 0.0
  starty = 0.0
 
  error = 0
  print *, "Reading grid file: ", trim (file_name)
 
  inquire (file=file_name, exist=fexist)
  if ( .not. fexist) then
    print *, "Cannot find file: ", trim (file_name)
    error = 1
    return
  end if
 
  open (13, file=file_name, status='old')
 
  i = 1
  do
    read (13, "(A)", iostat=stat) line
    if (stat < 0) exit
    line = adjustl (line)
    if (line(1:1) .ne. "#" .and. len(trim(line)) /= 0) then
 
      ipos = index (line, "NX")
      if (ipos > 0) then
        ipos = ipos + 2
        call value (line(ipos:), nx, err)
        if (nx > 0 .and. ny > 0) then
          ngrid = nx * ny
          allocate (lats(ngrid))
          allocate (lons(ngrid))
          allocate (elevs(ngrid))
          allocate (slp_n(ngrid))
          allocate (slp_e(ngrid))
        end if
      end if
      ipos = index (line, "NY")
      if (ipos > 0) then
        ipos = ipos + 2
        call value (line(ipos:), ny, err)
        if (nx > 0 .and. ny > 0) then
          ngrid = nx * ny
          allocate (lats(ngrid))
          allocate (lons(ngrid))
          allocate (elevs(ngrid))
          allocate (slp_n(ngrid))
          allocate (slp_e(ngrid))
        end if
      end if
      ipos = index (line, "DX")
      if (ipos > 0) then
        ipos = ipos + 2
        call value (line(ipos:), dx, err)
      end if
      ipos = index (line, "DY")
      if (ipos > 0) then
        ipos = ipos + 2
        call value (line(ipos:), dy, err)
      end if
      ipos = index (line, "STARTX")
      if (ipos > 0) then
        ipos = ipos + 6
        call value (line(ipos:), startx, err)
      end if
      ipos = index (line, "STARTY")
      if (ipos > 0) then
        ipos = ipos + 6
        call value (line(ipos:), starty, err)
      end if
    end if
 
    ipos = index (line, ',')
    if (ipos > 0 .and. ngrid > 0) then
      call parse (line, ",", settings, nsettings)
      if (nsettings == 5) then
        call value (settings(1), lats(i), err)
        if (err /= 0) lats (i) = - 999.99
        call value (settings(2), lons(i), err)
        if (err /= 0) lons (i) = - 999.99
        call value (settings(3), elevs(i), err)
        if (err /= 0) elevs (i) = - 999.99
        call value (settings(4), slp_n(i), err)
        if (err /= 0) slp_n (i) = - 999.99
        call value (settings(5), slp_e(i), err)
        if (err /= 0) slp_e (i) = - 999.99
           !print *, lats(i), lons(i), elevs(i)
        i = i + 1
      end if
    end if
 
  end do
  if (ngrid == 0) then
    print *, "Failed to find NX, NY in grid list: ", trim (file_name)
    error = 1
  else
    if (i /= ngrid+1) then
      print *, "Found only ", i, " out of ", ngrid, " points from: ", trim (file_name)
      error = 1
    end if
  end if
 
end subroutine read_grid_list



! previous version of routine w/ mixed format station metadata read ... can be deleted if not used
subroutine read_station_list_prev (file_name, id, name, lat, lon, elev, sslp_n, sslp_e, n_stations, &
                              & error, vars)
  use string_mod
  use type
  implicit none
 
  character (len=500), intent (in) :: file_name
  character (len=100), allocatable, intent (out) :: id (:), name (:)
  character (len=2), allocatable, intent (out) :: vars (:) ! AWW-feb2016 holds P/T flags
  real (dp), allocatable, intent (out) :: lat (:), lon (:), elev (:), sslp_n (:), sslp_e (:)
  integer (i4b), intent (out) :: n_stations
  integer, intent (out) :: error
 
  character (len=100) :: settings (8)
  integer (i4b) i, nsettings, stat, err, ipos
  character (500) line
  logical fexist
 
  error = 0
  print *, "Reading stations list: ", trim (file_name)
 
  inquire (file=file_name, exist=fexist)
  if ( .not. fexist) then
    print *, "Cannot find file: ", trim (file_name)
    stop
  end if
 
  open (11, file=file_name, status='old')
 
  n_stations = 0
  i = 1
  do
    read (11, "(A)", iostat=stat) line
    if (stat < 0) exit
    line = adjustl (line)
 
    if (line(1:1) .ne. "#" .and. len(trim(line)) /= 0) then
      ipos = index (line, "NSITES")
      if (ipos > 0) then
        ipos = ipos + 6
        call value (line(ipos:), n_stations, err)
        if (err /= 0) then
          n_stations = 0
        else
          print *, "Stations: ", n_stations
          allocate (id(n_stations))
          allocate (name(n_stations))
          allocate (lat(n_stations))
          allocate (lon(n_stations))
          allocate (elev(n_stations))
          allocate (sslp_n(n_stations))
          allocate (sslp_e(n_stations))
          allocate (vars(n_stations))  !AWW-feb2016 holds P/T flags
        end if
      end if
      ipos = index (line, ',')
      if (ipos > 0 .and. n_stations > 0) then
        call parse (line, ",", settings, nsettings)
        if (nsettings == 7 .or. nsettings == 8) then  ! orig number of fields
          id(i) = settings(1)
          name(i) = settings(7)
          call delall (name(i), '"')
          call value (settings(2), lat(i), err)
          if (err /= 0) lat (i) = - 999.99
          call value (settings(3), lon(i), err)
          if (err /= 0) lon (i) = - 999.99
          call value (settings(4), elev(i), err)
          if (err /= 0) elev (i) = - 999.99
          call value (settings(5), sslp_n(i), err)
          if (err /= 0) sslp_n (i) = - 999.99
          call value (settings(6), sslp_e(i), err)
          if (err /= 0) sslp_e (i) = - 999.99

          if (nsettings == 8) then ! AWW-feb2016
            vars(i) = settings(8)  ! need to upgrade station list
          end if
         
          ! print *, trim(id(i)), "  ", trim(name(i)), lat(i), lon(i), elev(i), vars[i]
          i = i + 1
        end if
 
     end if
    end if
 
  end do
  if (n_stations == 0) then
    print *, "Failed to find NSITES in station list: ", trim (file_name)
    error = 1
  else
    if (i /= n_stations+1) then
      print *, "Found only ", i, " out of ", n_stations, " stations from: ", trim (file_name)
      error = 1
    end if
  end if
  close (11)

  print *, "Done with read_station_list"

end subroutine read_station_list_prev
