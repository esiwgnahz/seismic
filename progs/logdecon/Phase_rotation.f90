! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Phase_rotation

  use Source_mod
  use sep

  implicit none

  real,    dimension(:,:), allocatable   :: data

  integer, dimension(:), pointer :: n                  ! data-size, dimension(ndim)
  real,    dimension(:), pointer :: d, o

  character(len=1024)    :: label
  integer :: ndim,npanels,i
  real :: phase

  call sep_init(SOURCE)

  ndim=sep_dimension()

  if (ndim.le.2) then
     allocate(n(3), o(3), d(3))
     do i=1,ndim
        call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
     end do
     n(ndim+1:)=1
     o(ndim+1:)=0.
     d(ndim+1:)=1
  else
     allocate(n(ndim), o(ndim), d(ndim))
     do i=1,ndim
        call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
     end do
  end if

  npanels=n(3)

  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:       Dimensions Input Data      '
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:'

  do i=1,ndim
     write(0,*) 'INFO: - - - - - - - '
     write(0,*) 'INFO: n(',i,')=',n(i)
     write(0,*) 'INFO: o(',i,')=',o(i)
     write(0,*) 'INFO: d(',i,')=',d(i)
  end do

  call from_param('phase',phase,90.) 
  phase=phase*3.14159265/180
     
  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO: Starting processing of ',npanels
  write(0,*) 'INFO: ------------------------------------'
  write(0,*) 'INFO:'

  allocate(data(n(1),n(2)))

  PANELS:do i=1,npanels

     call sreed('in',data,4*n(1)*n(2))
     call SourceRandom(data,n(1),d(1),n(2),phase)
     call srite('out',data,4*n(1)*n(2))
     
  end do PANELS
  
call sep_close()

end program PHASE_ROTATION
