subroutine generate_date_series (start_date, stop_date, date_series)

  IMPLICIT  NONE
  integer, intent(in):: start_date, stop_date ! start/stop date (yyyymmdd)
  integer, intent(out), allocatable :: date_series(:, :) ! (month numbers, year/month/start_day/stop_day)
  integer :: monthnum, start_y, stop_y, start_m, stop_m, start_d, stop_d, yy, mm, i, ntimes

  start_y = start_date / 10000
  start_m = mod(start_date, 10000) / 100
  start_d = mod(start_date, 100)
  stop_y = stop_date / 10000
  stop_m = mod(stop_date, 10000) / 100
  stop_d = mod(stop_date, 100)
  monthnum = (stop_y - start_y + 1) * 12 - (start_m - 1) - (12 - stop_m)
  
  
  allocate(date_series(monthnum, 4))
  
  if (monthnum .eq. 1) then
  	  date_series(1, 1) = start_y
  	  date_series(1, 2) = start_m
  	  date_series(1, 3) = start_d
  	  date_series(1, 4) = stop_d
  else
	  do i = 1, monthnum
		if (i .eq. 1) then
		    mm = start_m
		    yy = start_y
		else
		    mm = date_series(i-1, 2) + 1
		    if (mm .gt. 12) then
		        mm = 1
		        yy = date_series(i-1, 1) + 1
		    else
		    	yy = date_series(i-1, 1)  
		    end if
		end if
		! days in the month
		if ((mm .eq. 1) .or. (mm .eq. 3) .or. (mm .eq. 5) .or. (mm .eq. 7) .or. (mm .eq. 8) .or. (mm .eq. 10) .or. (mm .eq. 12)) then
			ntimes = 31
		end if
		if ((mm .eq. 4) .or. (mm .eq. 6) .or. (mm .eq. 9) .or. (mm .eq. 11)) then 
			ntimes = 30
		end if
		if (mm .eq. 2) then
			if (((mod(yy,4) .eq. 0) .and. (mod(yy,100) .ne. 0)) .or. (mod(yy,400) .eq. 0)) then
				ntimes = 29
			else
				ntimes = 28
			end if 
		end if 
		date_series(i, 1) = yy
		date_series(i, 2) = mm
		if (i .eq. 1) then
			date_series(i, 3) = start_d
		else
			date_series(i, 3) = 1
		end if
		if (i .eq. monthnum) then
			date_series(i, 4) = stop_d
		else
			date_series(i, 4) = ntimes	
		end if
	  end do
  end if

end subroutine generate_date_series
