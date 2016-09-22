module GenParam_types_flt

  implicit none

  type Genpar_flt
     character(len=3)  :: obstag
     character(len=3)  :: modtag

     integer :: ndim     ! Dimensions problem

     integer :: niter    ! Number of iterations to estimate the filters
     real    :: eps      ! Epsilong for regularization
     real    :: thresh_d ! Threshold Hyperbolic norm data space
     real    :: thresh_m ! Threshold Hyperbolic norm modelq space
  end type Genpar_flt

end module GenParam_types_flt
