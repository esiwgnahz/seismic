! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 

program Add_long

  use sep

  implicit none

  double precision,    dimension(:), allocatable:: file1,file2,file3

  integer, dimension(:), allocatable:: n1,n2
  real,    dimension(:), allocatable:: o1,o2
  real,    dimension(:), allocatable:: d1,d2
  character(len=1024)               :: label

  integer :: ndim,i,j

  call sep_init(SOURCE)
  ndim=sep_dimension()
  if (ndim.eq.2) ndim=ndim+1
  allocate(n1(ndim), o1(ndim), d1(ndim))
  allocate(n2(ndim), o2(ndim), d2(ndim))
  n1=0.
  n2=0.

  do i=1,ndim
     call sep_get_data_axis_par("in",i,n1(i),o1(i),d1(i),label)
  end do
  do i=1,ndim
     call sep_get_data_axis_par("file2",i,n2(i),o2(i),d2(i),label)
  end do

  n1(3)=max(n1(3),1)
  n2(3)=max(n2(3),1)

  write(0,*) n1
  write(0,*) n2
  if (product(n1).ne.product(n2)) call erexit('ERROR: files are not the same dimension, exit now')
  allocate(file1(product(n1)),file2(product(n1)),file3(product(n1)))

  write(0,*) 'INFO: Adding stdin and file2'
  write(0,*) 'INFO: n=',n1

  file3=0.d0;file2=0.d0;file1=0.d0

  call sreed('in',   file1,8*product(n1))
  call sreed('file2',file2,8*product(n1))

  file3=file2+file1
  
  write(0,*) sngl(file3)
  call srite('file3',file3,8*product(n1))
        
  call to_history('n1',n1(1),'file3')
  call to_history('n2',n1(2),'file3')
  call to_history('n3',n1(3),'file3')
  call to_history('esize',8,'file3')

end program Add_long
