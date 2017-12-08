program lbgs_prog

  use sep
  use nlinv_io_mod
  use nlinv_types_mod
  use nlinv_read_mod

  use LBFGS_doc_mod

  implicit none

  type(nlinvsepfile) :: invsep
  type(nlinvparam)   :: invpar
  type(nlinvarray)   :: invarr

  logical            :: sep_setup_works
  logical            :: lbfgs_setup_works

  character(len=1024):: logFilename,iterlogFilename
  logical            :: existLog

  ! LBFGS variables
  DOUBLE PRECISION   :: GTOL,STPMIN,STPMAX
  INTEGER            :: MP,LP
  external :: LB4
  COMMON /LB5/MP,LP,GTOL,STPMIN,STPMAX

  call sep_init()
  call LBFGS_doc()

  logFilename='log.dat'
  iterlogFilename='bfgsiter.dat'

  inquire(file=logFilename,exist=existLog)
  if(existLog) then
     open(1010,File=logFilename,status='old',position='append',action='write')
  else
     open(1010,File=logFilename,status='new',action='write')
  end if

  inquire(file=iterlogFilename,exist=existLog)
  if(existLog) then
     open(1020,File=iterlogFilename,status='old',position='append',action='write')
  else 
     open(1020,File=iterlogFilename,status='new',action='write')
  end if

  sep_setup_works=lbfgs_setup_sepfile(invsep)
  if(.not.sep_setup_works) stop 'SEP file setup failed, exit now'
  
  lbfgs_setup_works=lbfgs_setup(invpar,invarr,invsep)
  if(.not.lbfgs_setup_works) stop 'LBFGS setup failed, exit now'

  if(invsep%lbfgs_type.eq.0) then
     call LBFGSOLD(invpar%NDIM,invpar%MSAVE,invarr%xd,invarr%fd,invarr%gd,&
          &        .False.,invarr%diagd,invpar%iprint,invpar%EPS,&
          &        invpar%XTOL,invsep%stp_init,invarr%wd,invpar%iflag,invpar%myinfo,&
          &        invpar%lbfgsDatFilename,invpar%lbfgsDatFilenameOut,&
          &        invpar%mcsrchDatFilename,invpar%mcsrchDatFilenameOut)
  else
!     call LBFGS(invpar%NDIM,invpar%MSAVE,invarr%xd,invarr%fd,invarr%gd,&
!          &     .False.,invarr%diagd,invpar%iprint,invpar%EPS,&
!          &     invpar%XTOL,invsep%stp_init,invarr%wd,invpar%iflag,invpar%myinfo,&
!          &     invpar%lbfgsDatFilename, invpar%lbfgsDatFilenameOut,&
!          &     invpar%mcsrchDatFilename, invpar%mcsrchDatFilenameOut,&
!          &     invsep%xmin,invsep%xmax,invsep%stp1_opt,invsep%clip_type)
  end if

  write(1010,*) "Writing working array"
  open(89,file=invpar%workingArrayFilenameOut, status='replace', FORM='UNFORMATTED', BUFFERED='YES')
  rewind(89)
  write(89) invarr%wd
  close(89)

  write(1010,*) "Writing diag"
  open(88,file=invpar%diagFilenameOut,status='replace', FORM='UNFORMATTED', BUFFERED='YES')
  rewind(88)
  write(88) invarr%diagd
  close(88)
 
  write(1010,*) 'Writing iflag=', invpar%iflag
  open(30,file=invpar%iflagFilenameOut,status='replace')
  rewind(30)
  write(30,"(I2)") invpar%iflag
  close(30)

  write(1010,*) 'Writing myinfo=', invpar%myinfo
  open(20,file=invpar%myInfoFilenameOut,status='replace')
  rewind(20)
  write(20,"(I2)") invpar%myinfo
  close(20) 

  close(1010)
  close(1020)

  allocate(invsep%mod%array(invpar%NDIM))
  invsep%mod%array=sngl(invarr%xd)
  call write_sepfile(invsep%mod)

  call nlinvsepfile_clean(invsep)
  call nlinv_cleanarray(invarr)

end program lbgs_prog
