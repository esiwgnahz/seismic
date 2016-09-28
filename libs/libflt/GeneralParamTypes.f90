module GenParam_types_flt

  implicit none

  type Genpar_flt
     character(len=3)  :: obstag
     character(len=3)  :: modtag
     character(len=4)  :: fmodtag
     character(len=4)  :: filttag
     character(len=7)  :: filtpchtag

     integer :: ndim     ! Dimensions problem

     integer :: niter    ! Number of iterations to estimate the filters
     real    :: eps      ! Epsilong for regularization
     real    :: thresh_d ! Threshold Hyperbolic norm data space
     real    :: thresh_m ! Threshold Hyperbolic norm modelq space

     integer :: nlaplac  ! Number of filter coefficients for Laplace filter in 1D
                         ! Final number is 3xnlaplac for 3D and 2xnlaplac in 2D

     logical :: prec     ! Whether or not to use preconditioning 
  end type Genpar_flt

end module GenParam_types_flt
