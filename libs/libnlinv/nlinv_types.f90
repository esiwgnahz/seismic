module nlinv_types_mod

  use nlinv_io_mod

  implicit none

contains

  type nlinvparam

     ! Input file
     type(sepfile_type) :: mod ! Model vector
     type(sepfile_type) :: fct ! Objective function
     type(sepfile_type) :: gdt ! Gradient vector
     type(sepfile_type) :: modi! Starting model

     ! Optional weight file, set to 1 if not present
     type(sepfile_type) :: wgh ! Gradient weight

     ! Input file for LBFGS
     character(len=*) :: diagFilename
     character(len=*) :: iflagFilename
     character(len=*) :: myInfoFilename
     character(len=*) :: lbfgsDatFilename
     character(len=*) :: mcsrchDatFilename
     character(len=*) :: workingArrayFilename
     
     ! Output file for LBFGS
     character(len=*) :: diagFilenameOut
     character(len=*) :: iflagFilenameOut
     character(len=*) :: myInfoFilenameOut
     character(len=*) :: lbfgsDatFilenameOut
     character(len=*) :: mcsrchDatFilenameOut
     character(len=*) :: workingArrayFilenameOut

     ! Clipping parameter influencing behavior
     integer          :: clip_type
     ! Parameter controling how the gradient is scaled in first iteration
     integer          :: stp1_opt
     ! min and max values for each parameter
     double precision, dimension(:), allocatable :: xmin,xmax
     ! Number of parameters to invert for
     integer          :: nparams   
     ! Size of model space for 1 parameter
     integer          :: NDIM
     ! Number of gradient steps to keep in memory for limited-memory bfgs
     integer          :: MSAVE
     ! Size working array in double precision integer
     integer(kind=8)  :: NWORK
     ! BFGS constants for line search
     double precision :: XTOL
     double precision :: EPS       
     ! BFGS control parameters
     integer          :: iflag
     integer          :: myinfo

     integer,dimension(2) :: iprint

  end type nlinvparam

  type nlinvarray
     
     double precision, dimension(:), allocatable :: gd   ! Gradient double precision
     double precision, dimension(:), allocatable :: xd   ! Model double precision
     double precision, dimension(:), allocatable :: wd   ! Working array double precision
     double precision, dimension(:), allocatable :: diagd! Diagonal of inverse Hessian
     double precision                            :: fd   ! Objective function double precision

  end type nlinvarray

  subroutine nlinv_cleanarray(param)
    type(nlinvarray) ::       param

    if(allocated(param%gd))    deallocate(param%gd)
    if(allocated(param%xd))    deallocate(param%xd)
    if(allocated(param%wd))    deallocate(param%wd)
    if(allocated(param%diagd)) deallocate(param%diagd)

  end subroutine nlinv_cleanarray


end module nlinv_types_mod
