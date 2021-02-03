program simple_test
  ! gfortran -L/usr/lib -lSystem simple_test.f90 -o simple_test.exe
   IMPLICIT  NONE
  
  interface
    subroutine generate_date_series (start_date, stop_date, date_series)
     IMPLICIT  NONE
     integer, intent(in):: start_date, stop_date ! start/stop date (yyyymmdd)
     integer, intent(out), allocatable :: date_series(:, :) ! (month numbers, year/month/start_day/stop_day)
    end subroutine generate_date_series
  end interface
  
  
  integer :: date1 = 19910504
  integer :: date2 = 19920913
  integer :: num(2)
  integer :: i
  integer, allocatable :: date_series(:, :)
  
  call generate_date_series(date1, date2, date_series)
  
  num = shape(date_series)
  
  do i = 1, num(1)
  	print *, date_series(i, :)
  end do
  

end program simple_test


