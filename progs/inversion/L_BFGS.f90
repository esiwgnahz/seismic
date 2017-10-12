program lbgs_prog

  use sep
  use nlinv_io_mod
  use nlinv_types_mod
  use nlinv_read_mod

  implicit none
  
  type(nlinvsepfile) :: invsep
  type(nlinvparam)   :: invpar
  type(nlinvarray)   :: invarr

  logical            :: sep_setup_works
  logical            :: lbfgs_setup_works

  ! LBFGS variables
  DOUBLE PRECISION   :: GTOL,STPMIN,STPMAX
  INTEGER            :: MP,LP
  external :: LB2
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX

  call sep_init()

  sep_setup_works=lbfgs_setup_sepfile(invsep)
  if(.not.sep_setup_works) stop 'SEP file setup failed, exit now'
  
  lbfgs_setup_works=lbfgs_setup(invpar,invarr,invsep)
  if(.not.lbfgs_setup_works) stop 'LBFGS setup failed, exit now'

  if(invpar%lbfgs_type.eq.0) then
     call LBFGSOLD(invpar%NDIM,invpar%MSAVE,invarr%xd,invarr%fd,invarr%gd,&
          &        .False.,invarr%diagd,invpar%iprint,invpar%EPS,&
          &        invpar%XTOL,invarr%wd,invpar%iflag,invpar%myinfo,&
          &        invpar%lbfgsDatFilename,invpar%lbfgsDatFilenameOut,&
          &        invpar%mcsrchDatFilename,invpar%mcsrchDatFilenameOut)
  else
!     call LBFGS(invpar%NDIM,invpar%MSAVE,invarr%xd,invarr%fd,invarr%gd,&
!          &     .False.,invarr%diagd,invpar%iprint,invpar%EPS,&
!          &     invpar%lbfgsDatFilename, invpar%lbfgsDatFilenameOut,&
!          &     invpar%mcsrchDatFilename, invpar%mcsrchDatFilenameOut,&
!          &     invsep%xmin,invsep%xmax,invpar%stp1_opt,invpar%clip_type)
  end if

  write(0,*) "Writing working array"
  open(89,file=invpar%workingArrayFilenameOut, status='replace', FORM='UNFORMATTED', BUFFERED='YES')
  rewind(89)
  write(89) invarr%wd
  close(89)

  write(0,*) "Writing diag"
  open(88,file=invpar%diagFilenameOut,status='replace', FORM='UNFORMATTED', BUFFERED='YES')
  rewind(88)
  write(88) invarr%diagd
  close(88)
 
  write(0,*) 'Writing iflag=', invpar%iflag
  open(30,file=invpar%iflagFilenameOut,status='replace')
  rewind(30)
  write(30,*) invpar%iflag
  close(30)

  write(0,*) 'Writing myinfo=', invpar%myinfo
  open(20,file=invpar%myInfoFilenameOut,status='replace')
  rewind(20)
  write(20,*) invpar%myinfo
  close(20) 

  allocate(invsep%mod%array(invpar%NDIM))
  invsep%mod%array=sngl(invarr%xd)
  call write_sepfile(invsep%mod)

  call nlinvsepfile_clean(invsep)
  call nlinv_cleanarray(invarr)

end program lbgs_prog
