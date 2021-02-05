subroutine field_rand_nopointer (nspl1, nspl2, sp_wght, sp_sdev, sp_ipos, sp_jpos, sp_num, iorder1d, jorder1d, cfield)
! no longer use 

! ----------------------------------------------------------------------------------------
! Purpose:
!
!   Used to estimate a correlated field of random numbers for the basin gridpoints
!
! ----------------------------------------------------------------------------------------
!
! Modified:
!	Andy Newman: Aug 2013
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
! Structures Used:
!
!   (1) INPUTDAT2D
!   --------------
!   Correspondence between basins and gridpoints
!
!   (2) GRIDWEIGHT
!   --------------
!   Weights applied to each grid point to generate spatially correlated random fields
!
! ----------------------------------------------------------------------------------------
!
  use nrtype ! variable types (DP, I4B, etc.)
  use nr, only: gasdev ! Num. Recipies
  use inputdat2d ! use to relate basins to gridpoints
!
  use gridweight ! correlation structure  ! deleted by TGQ
  implicit none
!
! input
  integer (i4b), intent (in) :: nspl1 ! # points (1st spatial dimension)
  integer (i4b), intent (in) :: nspl2 ! # points (2nd spatial dimension)
  
  ! Add by TGQ
  real (dp), intent (in) :: sp_wght(:, :, :)
  real (dp), intent (in) :: sp_sdev(:, :)
  integer (i4b), intent (in) :: sp_ipos(:,:,:), sp_jpos(:,:,:)
  integer (i4b), intent (in) ::sp_num(:,:)
  integer (i4b), intent (in) :: iorder1d(:), jorder1d(:)
  
! Add by TGQ
  
!
  ! Output
  real (dp), dimension (nspl1, nspl2), intent (out) :: cfield ! correlated random field
!
!
  ! Internal variables
  integer (i4b) :: nlon ! number of x points in the 2-d grid
  integer (i4b) :: nlat ! number of y points in the 2-d grid
  integer (i4b) :: ierr ! error code for alolocate statement
  integer (i4b) :: igrd ! loop through the gridpoints
  integer (i4b) :: ilon, ilat ! (i,j) position of the igrd-th point
  integer (i4b) :: jlon, jlat ! (i,j) position of prev generated point
  integer (i4b) :: iprev ! loop thru previously generated points
  integer (i4b) :: nprv ! number of previously generated points
  integer (i4b) :: ilnk ! index to the nearest grid point
  real (dp), dimension (:), allocatable :: vprv ! vector of previously generated points
  real (dp) :: xbar ! conditional mean at each grid
  real (sp) :: aran ! a single random number
  real (dp), dimension (:, :), allocatable :: cran ! grid of correlated random numbers
  integer (i4b) :: irch ! loop through the stream segments

  ! ----------------------------------------------------------------------------------------
  ! (1) GET THE NUMBER OF X AND Y POINTS AND ALLOCATE SPACE FOR THE RANDOM GRID
  ! ---------------------------------------------------------------------------------------- 
  nlon = size (sp_num, 1)
  nlat = size (sp_num, 2)
  allocate (cran(nlon, nlat), stat=ierr)
  cran = -9999.0
  
  if (ierr .ne. 0) call exit_scrf (1, 'CRAN: PROBLEM ALLOCATING SPACE')
  ! ----------------------------------------------------------------------------------------
  ! (1) GENERATE GRID OF RANDOM NUMBERS
  ! ----------------------------------------------------------------------------------------
  ! loop through the grid points
  do igrd = 1, nlon * nlat
  ! identify the (i,j) position of the igrd-th point
    ilon = iorder1d (igrd)
    ilat = jorder1d (igrd)
  ! assign a random number to the first grid-point
    if ((ilon .gt. 0) .and. (ilat .gt. 0)) then ! whether the grid should be run
    
		if ((ilon .eq. 1) .and. (ilat .eq. 1)) then
		  call gasdev (aran)
		  cran (ilon, ilat) = aran
	  ! process gridpoints 2,...,n
		else
		! get the number of "previously generated points"
		  nprv = sp_num(ilon, ilat)
		  allocate (vprv(nprv), stat=ierr)
		  if (ierr .ne. 0) call exit_scrf (1, 'problem allocating array [field_rand.f90]')
		! build a vector of previously generated points
		  do iprev = 1, nprv
			jlon = sp_ipos(ilon, ilat, iprev)! i-position of previously generated point
			jlat = sp_jpos(ilon, ilat, iprev)! j-position of previously generated point
			vprv (iprev) = cran (jlon, jlat)! (previously generated point)
		  end do ! iprev
		! and generate the "current" point
		  call gasdev (aran)
! 		  sp_wght(ilon, ilat, 1:nprv) = sp_wght(ilon, ilat, 1:nprv) / sum(sp_wght(ilon, ilat, 1:nprv))
		  xbar = dot_product (vprv(1:nprv), sp_wght(ilon, ilat, 1:nprv))
		  cran (ilon, ilat) = xbar + sp_sdev(ilon, ilat) * aran  
	  
		! free up VPRV so that we can use it again for the next grid
		deallocate (vprv, stat=ierr)
		  if (ierr .ne. 0) call exit_scrf (1, 'problem deallocating array [field_rand.f90]')
		end if ! (if not the first point)
    
    end if 
  end do ! igrd
  
  cfield = cran
!
  return
!
end subroutine field_rand_nopointer
