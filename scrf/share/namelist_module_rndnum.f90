module namelist_module_rndnum
  use nrtype
  implicit none
 
  ! integer (i4b) :: nens !number of ensemble members to generate
  integer (i4b) :: start_ens, stop_ens ! start and stop numbers of ens. members to make
  integer (i4b) :: start_date, stop_date ! start/stop date: yyyymmdd
  integer (i4b) :: cross_cc_flag ! whether cross cc is used (>0: use. <0: don't use)

  real (dp) :: weight_judge ! if spcc weight > weight_judge, uniform parameter will be used instead of distributed parameter

  character(len=1024)   :: exp2p_file ! netcdf file of parameters for Exponential 2-parameter
  character (len=1024) :: grid_name !name of grid file
  character (len=1024) :: out_spcorr_prefix !base output name for spatial correlation structure
  character (len=1024) :: out_rndnum_prefix !base output name for spatial temporal correlated random numbers
  character (len=1024) :: cc_file !netcdf file of CC_lag1 or CC_cross
  character (len=1024) :: cross_file_prefix !file prefix of random numbers that have cross correlation with the target variable

  ! define namelist required variables
  ! namelist / params / nens, ntimes, grid_name, out_name_base, qpe_nc_name, clen, start_time
  namelist / params / start_date, stop_date, start_ens, stop_ens, cross_cc_flag, weight_judge, exp2p_file, &
    & cross_file_prefix, cc_file, grid_name, out_spcorr_prefix, out_rndnum_prefix
 
  save
contains
 
    ! AWW-16 - updated to process namelist file given as argument
  subroutine read_namelist_rndnum (namelist_filename)
    implicit none
 
      !input variables
    character (len=200), intent (in) :: namelist_filename !AWW-2016
 
      !local variables
    integer :: ierr
 
    open (unit=30, file=namelist_filename, form="FORMATTED")
 
    read (unit=30, nml=params, iostat=ierr)
    if (ierr /= 0) then
      write (*, '(/," ***** ERROR: Problem reading namelist PARAMS",/)')
      rewind (unit=30)
      read (unit=30, nml=params)
      stop " ***** ERROR: Problem reading namelist PARAMS"
    end if
 
    return
  end subroutine
 
end module namelist_module_rndnum
