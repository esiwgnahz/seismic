module nlinv_read_mod

  use sep
  use nlinv_types_mod
  use nlinv_io_mod

  implicit none

contains

  function lbfgs_setup_sepfile(nlinv_sepfile)
    implicit none
    type(nlinvsepfile) ::      nlinv_sepfile 
    real, dimensio(:), allocatable :: tmp

    lbfgs_setup_sepfile = .true.

    ! We set the tag for the sep files
    nlinv_sepfile%mod%tag='Model'
    nlinv_sepfile%fct%tag='Function'
    nlinv_sepfile%gdt%tag='Gradient'
    nlinv_sepfile%modi%tag='Initial_Model'
    nlinv_sepfile%wgh%tag='Gradient_Weight'

    ! We first read the sep files for model, gradient and function and 
    ! make sure they exist and check dimension consitencies.
    if (.not.exist_file(nlinv_sepfile%mod%tag)) then
       write(0,*) 'Model file is missing, exit now'
       lbfgs_setup_sepfile=.false.
       return
    else
       call read_sepfile(nlinv_sepfile%mod)
       allocate(nlinv_array%xd(product(nlinv_sepfile%mod%n)))
       nlinv_array%xd=dble(nlinv_sepfile%mod%array)
       deallocate(nlinv_sepfile%mod%array)
    end if

    if (.not.exist_file(nlinv_sepfile%fct%tag)) then
       write(0,*) 'Function file is missing, exit now'
       lbfgs_setup_sepfile=.false.
       return
    else
       call read_sepfile(nlinv_sepfile%fct)
       nlinv_array%fd=dble(nlinv_sepfile%fct%array(1))
    end if

    if (.not.exist_file(nlinv_sepfile%gdt%tag)) then
       write(0,*) 'Gradient file is missing, exit now'
       lbfgs_setup_sepfile=.false.
       return
    else
       call read_sepfile(nlinv_sepfile%gdt)
       allocate(nlinv_array%gd(product(nlinv_sepfile%gdt%n)))
       nlinv_array%gd=dble(nlinv_sepfile%gdt%array)
       deallocate(nlinv_sepfile%gdt%array)
    end if

    if (.not.exist_file(nlinv_sepfile%modi%tag)) then
       write(0,*) 'Model init file is missing, exit now'
       lbfgs_setup_sepfile=.false.
       return
    else
       call read_sepfile(nlinv_sepfile%modi)
    end if

    if (.not.exist_file(nlinv_sepfile%wgh%tag)) then
       ! We set the weight to one if the file doesn't exist
       allocate(nlinv_sepfile%wgh%array(product(nlinv_sepfile%gdt%n)))
       nlinv_sepfile%wgh%array=1.
    else
       call read_sepfile(nlinv_sepfile%wgh)
    end if

    ! Do mod and grad have the same size?
    if (nlinv_sepfile%gdt%n.ne.nlinv_sepfile%mod%n) then
       write(0,*) 'Gradient and Model have different sizes, exit now'
       lbfgs_setup_sepfile=.false.
       return
    end if

    ! Do gdt and wgh have the same size?
    if (nlinv_sepfile%gdt%n.ne.nlinv_sepfile%wgh%n) then
       write(0,*) 'Gradient and Gradient weight have different sizes, exit now'
       lbfgs_setup_sepfile=.false.
       return
    end if

    ! Do mod and modi have the same size?
    if (nlinv_sepfile%mod%n.ne.nlinv_sepfile%modi%n) then
       write(0,*) 'Model and initial model have different sizes, exit now'
       lbfgs_setup_sepfile=.false.
       return
    end if

    ! Set the size for xmin and xmax: last axis is usually the nparam axis
    allocate(nlinv_sepfile%xmin(size(nlinv_sepfile%mod%n))
    allocate(nlinv_sepfile%xmax(size(nlinv_sepfile%mod%n))
    
    allocate(tmp(size(nlinv_sepfile%xmin))); tmp=0.
    call from_param('xmin',nlinv_sepfile%xmin,tmp)
    call from_param('xmax',nlinv_sepfile%xmax,tmp)
    deallocate(tmp)

  end function lbfgs_setup_sepfile

  function lbfgs_setup(nlinv_param,nlinv_array,nlinv_sepfile)
    implicit none
    type(nlinvparam) ::nlinv_param
    type(nlinvarray) ::            nlinv_array
    type(nlinvsepfile) ::                      nlinv_sepfile 
    logical :: lbfgs_setup
    logical :: iflagExists
    logical :: myInfoExists

    lbfgs_setup=.true.

    ! We set the filenames for the lbfgs files
    nlinv_param%diagFilename='Diag.dat'
    nlinv_param%iflagFilename='Iflag.dat'
    nlinv_param%myInfoFilename='Info.dat'
    nlinv_param%lbfgsDatFilename='LbfgsDat.dat'
    nlinv_param%mcsrchDatFilename='McsrchDat.dat'
    nlinv_param%workingArrayFilename='WorkArray.dat'

    nlinv_param%diagFilenameOut='DiagOut.dat'
    nlinv_param%iflagFilenameOut='IflagOut.dat'
    nlinv_param%myInfoFilenameOut='InfoOut.dat'
    nlinv_param%lbfgsDatFilenameOut='LbfgsDatOut.dat'
    nlinv_param%mcsrchDatFilenameOut='McsrchDatOut.dat'
    nlinv_param%workingArrayFilenameOut='WorkArrayOut.dat'

    ! Now setup iflag and myinfo, which control the behavior of lbfgs
    inquire( file=nlinv_param%iflagFilename, exist=iflagExists)
    if (iflagExists) then
       open(10,file=nlinv_param%iflagFilename)
       read(10,*) nlinv_param%iflag
       close(10)
       write(0,*) 'Reading iflag=', nlinv_param%iflag
    else
       nlinv_param%iflag = 0	
       write(0,*) 'iflag control file not present. Setting iflag=', nlinv_param%iflag
    end if
       
    inquire( file=nlinv_param%myInfoFilename, exist=myInfoExists)
    if (myInfoExists) then
       open(40,file=nlinv_param%myInfoFilename)
       read(40,*) nlinv_param%myinfo
       close(40)
       write(0,*) 'Reading myinfo= ', nlinv_param%myinfo
    else
       nlinv_param%myinfo = 0	
       write(0,*) 'myinfo control file not present. Setting myinfo=', nlinv_param%myinfo
    end if

    ! Set some parameters

    ! Dimension of model space
    nlinv_param%NDIM   =product(nlinv_sepfile%gdt%n)
    ! The last axis is usually the dimension of the number of parameters to invert for
    nlinv_param%nparam=nlinv_sepfile%gdt%n(size(nlinv_sepfile%gdt%n))
    ! We keep the history of the last 5 iterations for the inverse Hessian
    nlinv_param%MSAVE  =5 
    ! Set size of working array
    nlinv_param%NWORK  =nlinv_param%NDIM*(2*nlinv_param%MSAVE+1)+2*nlinv_param%MSAVE
    ! Set some constants for line search
    nlinv_param%XTOL = epsilon(XTOL)
    nlinv_param%EPS  = 1.0D-20

    ! Printing parameters
    nlinv_param%iprint(1)= 1
    nlinv_param%iprint(2)= 0
    
    allocate(nlinv_array%wd(NWORK))
    allocate(nlinv_array%diagd(NDIM))

    if (nlinv_param%iflag.ne.0) then
       write(0,*) 'Reading diag'
       open(88,file=nlinv_param%diagFilename, FORM='UNFORMATTED')
       rewind(88)
       read(88) nlinv_array%diagd
       close(88)
       
       write(0,*) 'Reading working array'
       open(89,file=nlinv_param%workingArrayFilename, FORM='UNFORMATTED')
       rewind(89)
       read(89) nlinv_array%wd
       close(89)
    end if

    !*************************************************
    ! Reading parameter values using seplib from param
    !*************************************************
    call from_param('lbfgs_type',nlinv_param%lbfgs_type,1)
    call from_param('bound_type',nlinv_param%clip_type,1)
    call from_param('stp1_opt',nlinv_nparam%stp1_opt,1)
    !*************************************************

  end function lbfgs_setup

end module nlinv_read_mod
