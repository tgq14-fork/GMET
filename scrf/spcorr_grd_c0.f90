subroutine spcorr_grd_c0 (nspl1, nspl2, c0, grid)
! ----------------------------------------------------------------------------------------
! Creator:
!   Martyn Clark, 2006, 2008
!
! Modified:
!    Andy Newman, Aug 2013
!    Andy Wood, Apr 2016
! ----------------------------------------------------------------------------------------
! Purpose:
!   Used to define weights needed to compute a spatially correlated grid of random numbers
!    (use the precipitation grid as a template)
!
! ----------------------------------------------------------------------------------------
! Method:
!
!   Based on the matrix-decomposition method described by Clark and Slater (JHM 2006).
!
!   Here the precipitation grid is used as a template.  Later, random fields generated
!   on the precipitation grid are used to perturb precipitation (on the precipitation grid,
!   before it is disaggregated and assigned to individual sub-basins).
!
!   Random fields used to perturb model states are generated by assigning random numbers
!   from the precipitation grid to individual sub-basins.
!
! ----------------------------------------------------------------------------------------
! I/O:
!
!   Input(s):
!   ---------
!    NSPL1: number of location points (1st spatial dimension)
!    NSPL2: number of location points (2nd spatial dimension)
!
! ----------------------------------------------------------------------------------------
! Structures Populated:
!
!   SPCORR in MODULE gridweight.f90
!    * includes the weights and variance assigned to each grid
!
! ----------------------------------------------------------------------------------------
! Future revisions:
!
!   This routine takes alot of time (and memory) for large grids.  Since many gridpoints
!   use the same combination of previously generated points (i.e., w.r.t. their relative
!   location compared to the target point), it should be possible to identify different
!   categories of gridpoints, and make computations for each of those categories.
!
! ----------------------------------------------------------------------------------------
  use nrtype ! variable types (DP, I4B, etc.)
  use nr, only: ludcmp, lubksb ! Num. Recipies
  use nrutil, only: arth ! Num. Recipies utilities
  use trig_degrees, only: sind, cosd ! added for the gfortran compiler
  use linkstruct ! linkage structures
  use gridweight ! grid correlation structure
  use namelist_module, only: clen ! correlation length
!
  implicit none
! input
  integer (i4b), intent (in) :: nspl1 ! # points (1st spatial dimension)
  integer (i4b), intent (in) :: nspl2 ! # points (2nd spatial dimension)
  real (dp), intent (in) :: c0
  type (coords), intent (in) :: grid ! input coordniate structure containing grid information
! define hyper parameters
  integer (i4b) :: nnst ! number of nests
  integer (i4b) :: nloc ! number of local points used
!REAL(DP)                                     :: CLEN        ! correlation length
! define looping variables
  integer (i4b) :: ires ! loop through resolution increments
  integer (i4b) :: incr ! spacing between grid cells
  integer (i4b) :: isp1, isp2 ! loop through lat-lon grid
  integer (i4b) :: jsp1, jsp2 ! loop through local grid
  integer (i4b) :: iprc ! counter for # points generated
! identify previously generated points
  integer (i4b) :: k ! counter for prev generated points
  integer (i4b) :: maxp ! max # of prev generated points
  integer (i4b) :: iprev, jprev ! loop through previously gen points
  integer (i4b) :: ierr ! error code for allocate statements
  logical (lgt), dimension (:, :), allocatable :: gmsk ! x-y mask
  integer (i4b) :: npts ! # previously generated points
  integer (i4b), dimension (:), allocatable :: ipos, jpos ! (i,j) index of prev generated points
! compute correlation between points
  real (dp) :: lat1, lon1 ! lat-lon of 1st prev generated point
  real (dp) :: lat2, lon2 ! lat-lon of 2nd prev generated point
  real (dp) :: dist ! distance between points (km)
  real (sp), dimension (:, :), allocatable :: corr ! correlation among prev generated points
  real (sp), dimension (:), allocatable :: gvec ! corr btw prev gen pts & current point
! compute weights and variance
  real (sp), dimension (:), allocatable :: twgt ! copy of GVEC
  integer (i4b), dimension (:), allocatable :: indx ! row permutation (ludcmp)
  real (sp) :: tmp ! identifier (ludcmp)
  real (dp), dimension (:), allocatable :: wght ! weight applied to each previous point
  real (dp) :: sdev ! standard deviation of estimate
  integer (i4b) :: n ! # of prev generated points
  logical (lgt), save :: init = .true. ! true for first call to generate_rndnum_c0
!
! output (none) -- structures updated in gridweight.f90
! ----------------------------------------------------------------------------------------
! (0) CHECK THAT SPCORR IS NOT POPULATED ALREADY
! ----------------------------------------------------------------------------------------
  if ( .not. init) then
    if (associated(spcorr)) return
  end if
! ----------------------------------------------------------------------------------------
! (1) DEFINE HYPER-PARAMETERS
! ----------------------------------------------------------------------------------------
  nnst = 10 ! number of nests
  nloc = 3 ! number of local points to include in the estimation
!
! ----------------------------------------------------------------------------------------
! (2) ALLOCATE SPACE FOR OUTPUT ARRAYS
! ----------------------------------------------------------------------------------------
! define the maximum number of previously generated points
  maxp = (nloc*2+1) ** 2
! allocate space for the mask
  if (allocated(gmsk)) then
    deallocate (gmsk, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for mask ')
  end if
  allocate (gmsk(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, ' problem allocating space for mask ')
! nullify pointers
  if (init) nullify (spcorr, iorder, jorder)! in MODULE gridweight
! allocate space for the correlation structure (SPCORR has a pointer structure)
  if (associated(spcorr)) then
    deallocate (spcorr, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for the spatial correlation st&
   &ructure ')
  end if
  allocate (spcorr(nspl1, nspl2), stat=ierr)
  if (ierr .ne. 0) call exit_scrf (1, ' problem allocating spatial corr struct ')
  if (init) then
    do isp1 = 1, nspl1
      do isp2 = 1, nspl2
        nullify (spcorr(isp1, isp2)%ipos, spcorr(isp1, isp2)%jpos, spcorr(isp1, isp2)%wght)
      end do
    end do
  end if
! allocate space for the processing order
  if (associated(iorder)) then
    deallocate (iorder, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for the processing order ')
  end if
  if (associated(jorder)) then
    deallocate (jorder, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for the processing order ')
  end if
  allocate (iorder(nspl1*nspl2), jorder(nspl1*nspl2), stat=ierr) ! iorder is row number
  if (ierr .ne. 0) call exit_scrf (1, ' problem allocating space for the processing order ')
! ----------------------------------------------------------------------------------------
! (3) LOOP THROUGH THE DIFFERENT GRID RESOLUTIONS (PROCESS COARSE RESOLUTION FIRST)
! ----------------------------------------------------------------------------------------
  sdev = 0D0 ! Initialize SDEV (added to account for the first point)  EÖH
  iprc = 0 ! counter for the number of grid points processed
  gmsk = .false. ! initialize a logical array to identify which points have been processed
  do ires = nnst - 1, 0, - 1
    incr = 2 ** ires ! increment (2**4=16, 2**3=8, 2**2=4, 2**1=2, 2**0=1)
    print *, 'Working on Loop: ', ires
 ! ---------------------------------------------------------------------------------------
 ! (4) LOOP THROUGH THE LAT-LON OF THE GRID AT A GIVEN RESOLUTION
 ! ---------------------------------------------------------------------------------------
    do isp1 = 1, nspl1, incr
      do isp2 = 1, nspl2, incr
   ! check that "current" point has not been generated yet
        if ( .not. gmsk(isp1, isp2)) then
    ! allocate space to store the (i,j) position, and weights
          if (allocated(ipos)) then
            deallocate (ipos, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating ipos ')
          end if
          if (allocated(jpos)) then
            deallocate (jpos, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating jpos ')
          end if
          if (allocated(wght)) then
            deallocate (wght, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating wght ')
          end if
          allocate (ipos(maxp), jpos(maxp), stat=ierr)
          if (ierr .ne. 0) call exit_scrf (1, ' problem allocating (i,j) position ')
          allocate (wght(maxp), stat=ierr)
          if (ierr .ne. 0) call exit_scrf (1, ' problem allocating space for weights ')
    ! increment IPRC
          iprc = iprc + 1
    ! save the (i,j) position of iprc
          iorder (iprc) = isp1
          jorder (iprc) = isp2
    ! ------------------------------------------------------------------------------------
    ! (5) IDENTIFY PREVIOUSLY GENERATED POINTS
    ! ------------------------------------------------------------------------------------
          k = 0 ! initialize the number of previous points generated to zero
    ! loop through points in the local neighbourhood
          do jsp1 = max (1, isp1-(incr*nloc)), min (isp1+(incr*nloc), nspl1)
            do jsp2 = max (1, isp2-(incr*nloc)), min (isp2+(incr*nloc), nspl2)
      ! check to see if the "local" point has been generated previously
              if (gmsk(jsp1, jsp2)) then ! local point has been previously generated
                k = k + 1
       ! save the (i,j) position of the previously generated point
                ipos (k) = jsp1
                jpos (k) = jsp2
              end if ! if the point has been previously generated
            end do ! jsp1
          end do ! jsp2
    ! include the (i,j) of the current point
          k = k + 1
          ipos (k) = isp1
          jpos (k) = isp2
    ! ...and save the number of points
          npts = k
    ! check that there are at least two points
          if (k .ge. 2) then
     ! ------------------------------------------------------------------------------------
     ! (6) COMPUTE THE CORRELATION AMONG PREVIOUSLY GENERATED POINTS
     ! ------------------------------------------------------------------------------------
            if (allocated(corr)) then
              deallocate (corr, stat=ierr)
              if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating corr ')
            end if
            if (allocated(gvec)) then
              deallocate (gvec, stat=ierr)
              if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating gvec ')
            end if
            if (allocated(twgt)) then
              deallocate (twgt, stat=ierr)
              if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating twgt ')
            end if
            if (allocated(indx)) then
              deallocate (indx, stat=ierr)
              if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating indx ')
            end if
            allocate (corr(k-1, k-1), gvec(k-1), twgt(k-1), indx(k-1), stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem allocating space for corr, gvec, twgt, or&
           & indx ')
     ! Note that the vector of previously generated points includes the current point as its
     ! last element.  The correlation among all previously generated points are computed over
     ! elements (1...k-1) and saved in the matrix corr.  The correlation between previously
     ! generated points (1...k-1) and the current point (k) is saved in ther vector gvec.
            do iprev = 1, k
              do jprev = 1, iprev
                if (iprev .eq. jprev) then
                  if (iprev .le. k-1) corr (iprev, jprev) = 1.0D0 ! grid points are the same, correlation=1
                else
        ! identify lat-lon of previously generated points
                  lon1 = grid%lon (ipos(iprev), jpos(iprev))! NOTE, iprev, lon
                  lon2 = grid%lon (ipos(jprev), jpos(jprev))! NOTE, jprev, lon
                  lat1 = grid%lat (ipos(iprev), jpos(iprev))! NOTE, iprev, lat
                  lat2 = grid%lat (ipos(jprev), jpos(jprev))! NOTE, jprev, lat
        ! compute distance (km) - on the surface of a sphere
                  dist = 6378.0D0 * acos &
                 & (sind(lat1)*sind(lat2)+cosd(lat1)*cosd(lat2)*cosd(lon1-lon2))
        ! compute correlation
                  if (iprev .le. k-1) then
         ! correlation among all previously generated points (1...k-1,1...k-1) -- corr
                    corr (iprev, jprev) = c0 * exp (-(dist/clen))
                    corr (jprev, iprev) = corr (iprev, jprev)
                  else
         ! correlation between all previously generated points and the current point -- gvec
                    if (jprev .le. k-1) gvec (jprev) = c0 * exp (-(dist/clen))
                  end if ! switch between corr and gvec
                end if ! if the points are the same
              end do ! jprev
            end do ! iprev
     ! ------------------------------------------------------------------------------------
     ! (7) COMPUTE THE WEIGHTS
     ! ------------------------------------------------------------------------------------
     ! Note that the vector of previously generated points includes the current point as its
     ! last element.  The correlation among all previously generated points are computed over
     ! elements (1...k-1) and saved in the matrix corr.  The correlation between previously
     ! generated points (1...k-1) and the current point (k) is saved in ther vector gvec.
     ! special case of the bi-variate normal
            if (k .eq. 2) then
              wght (1) = gvec (1)
              sdev = sqrt (1.-gvec(1)**2.)
     ! all other points
            else
      ! temporary weight (GVEC is over-written)
              twgt (1:k-1) = gvec (1:k-1)
      ! estimate weights
              call ludcmp (corr, indx, tmp)
              call lubksb (corr, indx, twgt)
      ! save weights and variance
              wght (1:k-1) = twgt (1:k-1)
              sdev = sqrt (1.-dot_product(gvec(1:k-1), twgt(1:k-1)))
            end if ! ( if k gt 2 )
     ! deallocate correlation arrays
            if (allocated(corr)) then
              deallocate (corr, stat=ierr)
              if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating corr ')
            end if
            if (allocated(gvec)) then
              deallocate (gvec, stat=ierr)
              if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating gvec ')
            end if
            if (allocated(twgt)) then
              deallocate (twgt, stat=ierr)
              if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating twgt ')
            end if
            if (allocated(indx)) then
              deallocate (indx, stat=ierr)
              if (ierr .ne. 0) call exit_scrf (1, ' probelm deallocating indx ')
            end if
          end if ! ( if k ge 2 )
    ! flag the point as "generated"
          gmsk (isp1, isp2) = .true.
    ! -------------------------------------------------------------------------------------
    ! (8) SAVE WEIGHTS IN THE SPATIAL CORRELATION STRUCTURE
    ! -------------------------------------------------------------------------------------
    ! allocate space for the (i,j) position of previously generated points
          if (associated(spcorr(isp1, isp2)%ipos)) then
            deallocate (spcorr(isp1, isp2)%ipos, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for ipos, jpos ')
          end if
          if (associated(spcorr(isp1, isp2)%jpos)) then
            deallocate (spcorr(isp1, isp2)%jpos, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for ipos, jpos ')
          end if
          allocate (spcorr(isp1, isp2)%ipos(npts-1), spcorr(isp1, isp2)%jpos(npts-1), stat=ierr)
          if (ierr .ne. 0) call exit_scrf (1, ' problem allocating space for the (i,j) position of &
         &previously generated points ')
    ! allocate space for the weights assigned to previously generated points
          if (associated(spcorr(isp1, isp2)%wght)) then
            deallocate (spcorr(isp1, isp2)%wght, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for weights assigned t&
           &o previously generated points ')
          end if
          allocate (spcorr(isp1, isp2)%wght(npts-1), stat=ierr)
          if (ierr .ne. 0) call exit_scrf (1, ' problem allocating space for weights assigned to pr&
         &eviously generated points ')
    ! populate the structures (-1 excludes the current (i,j) point)
          spcorr(isp1, isp2)%ipos(1:npts-1) = ipos (1:npts-1)! i-position
          spcorr(isp1, isp2)%jpos(1:npts-1) = jpos (1:npts-1)! j-position
          spcorr(isp1, isp2)%wght(1:npts-1) = wght (1:npts-1)! weight assigned to previously generated points
          spcorr(isp1, isp2)%sdev = sdev ! standard deviation of estimate
!
    ! deallocate (i,j) position and weights
          if (allocated(ipos)) then
            deallocate (ipos, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating ipos ')
          end if
          if (allocated(jpos)) then
            deallocate (jpos, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating jpos ')
          end if
          if (allocated(wght)) then
            deallocate (wght, stat=ierr)
            if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating wght` ')
          end if
        end if ! if the point has not been generated yet
      end do ! ilat
    end do ! ilon
  end do ! ires
! deallocate mask
  if (allocated(gmsk)) then
    deallocate (gmsk, stat=ierr)
    if (ierr .ne. 0) call exit_scrf (1, ' problem deallocating space for mask ')
  end if
!IF (INIT) INIT=.FALSE.
!
! ----------------------------------------------------------------------------------------
  return
end subroutine spcorr_grd_c0
