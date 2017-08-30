! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program MakeSegyHeader

  use sep 
  implicit none

  double precision, dimension(:,:), allocatable:: header

  integer, dimension(:),   allocatable:: n
  real,    dimension(:),   allocatable:: o
  real,    dimension(:),   allocatable:: d
  character(len=1024),    dimension(:), allocatable  :: label

  integer :: n23,ndim,i,j,k,l,index,indexarray
  real    :: scaling

  call sep_init(SOURCE)
  ndim=sep_dimension()

  allocate(n(3), o(ndim), d(ndim), label(ndim))
  n=1
  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label(i))
  end do

  call from_param('scaling',scaling,100000.)
  
  n23=product(n(2:))

  allocate(header(n23,ndim-1))

  write(0,*) 'INFO: There are ',n23,' traces in this file'

  do i=1,ndim-1
     write(0,*) 'INFO: Writing header information for',i,'/',ndim-1,' label=',label(i)
     index=0
     do k=1,n(3)
        do l=1,n(2)
           if (i==1) indexarray=l
           if (i==2) indexarray=k
           index=index+1
           header(index,i)=(dble(scaling)*dble((o(i+1)+(indexarray-1)*d(i+1))))
           !write(0,*) index,i,header(index,i)
        end do
     end do
  end do
  
  call to_history('n1',n23,'header')
  call to_history('n2',ndim-1,'header')
  call to_history('n3',1,'header')
  call srite('header',sngl(header),4*n23*(ndim-1))

end program MakeSegyHeader

