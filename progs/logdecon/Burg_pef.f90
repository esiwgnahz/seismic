! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
Program Burg_pef
  
  !use tcai1
  use sep
  
  implicit none
  
  integer :: n1,n2,n3,n123,na,ny,stat
  real, dimension(:), allocatable :: d,dpef,pef
  real, dimension(:), pointer     :: bb
  real, dimension(:), allocatable :: yy

  call sep_init()

  call from_param   ("na",na,10)
  call from_history ('n1', n1)
  call from_history ('n2', n2,1)
  call from_history ('n3', n3,1)

  n123=n1*n2*n3

  write (0,*) 'INFO: Number of data points =',n123
  write (0,*) 'INFO: Size filter           =',na

  ny=n123+na-1

  call to_history('n1',na,'filter')
  call to_history('n2',1,'filter')
  call to_history('n3',1,'filter')

  allocate (d(n123),dpef(n123),pef(na),yy(ny),bb(na))

  call sep_read(d)
  call burgc(d,dpef,pef,n123,na)
  write(0,*) 'INFO: Done with pef estimation'
  call srite('filter',pef,4*na)
  !bb=pef
  !call tcai1_init(bb)
  !stat=tcai1_lop(.false.,.false.,d,yy)
  !write(0,*) 'INFO: Done with convolution'
  !call srite('out',yy,4*n123)
  call srite('out',dpef,4*n123)
  deallocate(d,dpef,pef,yy)
  call sep_close()

end program Burg_pef

! Burg pef estimation from Jon's book
subroutine burgc(input,output,filter,n123,na)
  implicit none
  real,intent(in):: input(n123)
  real,intent(out):: filter(na),output(n123)
  integer ::lx,lc,i,j,n123,na
  real,allocatable::x(:),ep(:),em(:),c(:),a(:),s(:)
  real ::bot,top,epi

  lx=n123
  lc=na

  allocate(x(lx),ep(lx),em(lx),c(lc),a(lc),s(lx))

  s=0.0

  c=0.0
  a=0.0
  a(1)=1.0
  x=input
  em=x
  ep=x
  do j=2,lc
    top=0.0
    bot=0.0
    do i=j,lx
      bot=bot+ep(i)*ep(i)+em(i-j+1)*em(i-j+1)
      top=top+ep(i)*em(i-j+1)
    end do
    
    c(j)=2.0D+0*top/(bot+epsilon(bot))
    do i=j,lx
      epi=ep(i)
      ep(i)=ep(i)-c(j)*em(i-j+1)
      em(i-j+1)=em(i-j+1)-c(j)*epi
    end do
    a(j)=0.0
    do i=1,j
      s(i)=a(i)-c(j)*a(j-i+1)
    end do
    do i=1,j
      a(i)=s(i)
    end do
  end do

  output=ep
  filter=a

  deallocate(x,ep,em,c,a,s)

end subroutine burgc
