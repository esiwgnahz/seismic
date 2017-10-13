module nlinv_types_mod

  use nlinv_io_mod

  implicit none

  type nlinvsepfile
     
     !*****************************************************************
     ! These arrays are read using seplib commands sreed and from_param
     !*****************************************************************
     ! Input file
     type(sepfile_type) :: mod ! Model vector
     type(sepfile_type) :: fct ! Objective function
     type(sepfile_type) :: gdt ! Gradient vector
     type(sepfile_type) :: modi! Starting model

     ! Optional weight file, set to 1 if not present
     type(sepfile_type) :: wgh ! Gradient weight
    
     ! min and max values for each parameter
     real, dimension(:), allocatable :: xmin, xmax

     !************************************************************************
     ! These parameters are read from the command line using seplib from_param
     !************************************************************************
     ! Clipping parameter influencing behavior
     integer          :: clip_type
     ! Parameter controling how the gradient is scaled in first iteration
     integer          :: stp1_opt
     ! Parameter controling the lbfgs version, with (1) without (0) clipping
     ! and masking. With clipping is for FWI for instance. Without clipping 
     ! is for LSRTM for instance.
     !************************************************************************
     integer          :: lbfgs_type

  end type nlinvsepfile

  type nlinvparam

     ! Input file for LBFGS
     character(len=1024) :: diagFilename
     character(len=1024) :: iflagFilename
     character(len=1024) :: myInfoFilename
     character(len=1024) :: lbfgsDatFilename
     character(len=1024) :: mcsrchDatFilename
     character(len=1024) :: workingArrayFilename
     
     ! Output file for LBFGS
     character(len=1024) :: diagFilenameOut
     character(len=1024) :: iflagFilenameOut
     character(len=1024) :: myInfoFilenameOut
     character(len=1024) :: lbfgsDatFilenameOut
     character(len=1024) :: mcsrchDatFilenameOut
     character(len=1024) :: workingArrayFilenameOut

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

contains

  subroutine nlinv_cleanarray(param)
    type(nlinvarray) ::       param

    if(allocated(param%gd))    deallocate(param%gd)
    if(allocated(param%xd))    deallocate(param%xd)
    if(allocated(param%wd))    deallocate(param%wd)
    if(allocated(param%diagd)) deallocate(param%diagd)

  end subroutine nlinv_cleanarray

  subroutine nlinvsepfile_clean(param)
    type(nlinvsepfile)   ::     param

    call deallocate_sepfile(param%mod)
    call deallocate_sepfile(param%fct)
    call deallocate_sepfile(param%gdt)
    call deallocate_sepfile(param%modi)
    call deallocate_sepfile(param%wgh)

    if(allocated(param%xmin)) deallocate(param%xmin)
    if(allocated(param%xmax)) deallocate(param%xmax)
    
  end subroutine nlinvsepfile_clean

end module nlinv_types_mod
