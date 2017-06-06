! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Taper

  implicit none
  integer			                           	:: n1,n2,n3,tap,n2t,i,j,left,right
  real			                          		:: o1,o2,o3,d1,d2,d3,thet,pi,it,n2b
  real,    dimension (:,:),   allocatable :: in,ou

  pi = 3.14159265

  call initpar()

  !################################
  ! get parameters from .H file
  !################################

  call hetch('n1','i',n1)
  call hetch('n2','i',n2)
  call hetch('n3','i',n3)
  call hetch('o1','r',o1)
  call hetch('o2','r',o2)
  call hetch('o3','r',o3)
  call hetch('d1','r',d1)
  call hetch('d2','r',d2)
  call hetch('d3','r',d3)

  !################################
  ! get parameter from par
  !################################

  left=1
  right=1
  call getch("tap","i",tap);write(0,*)'from par tap',tap
  call getch("left","i",left);write(0,*) 'from par left',left
  call getch("right","i",right);write(0,*) 'from par right',right
  if (tap >= 50) call erexit ('Incorrect tap value')

  !################################
  ! memory allocation/reading input
  !################################

  allocate (in (n1,n2))
  allocate (ou (n1,n2))
  call sreed('in',in,size(in)*4)

  !################################
  ! put input dimensions in output
  !################################

  call putch ("n1","i",n1)
  call putch ("n2","i",n2)
  call putch ("n3","i",n3)
  call putch ("o1","r",o1)
  call putch ("o2","r",o2)
  call putch ("o3","r",o3)
  call putch ("d1","r",d1)
  call putch ("d2","r",d2)
  call putch ("d3","r",d3)


  !################################
  ! taper the data with cos^2
  !################################

  n2t=int(n2*tap/100)
  if (2*n2t >= n2) call erexit ('Incorrect tap value')
  n2b = n2t

  ou = in

  if (left.eq.1) then
     do j=1,n1
        do i=1,n2t
           it      = i
           thet    = ((it-1)/(n2b-1))*(pi/2)
           ou(j,i) = in(j,i)*sin(thet)*sin(thet)
        end do
     end do
  end if
  
  if (right.eq.1) then
     do j=1,n1
        do i=n2-n2t,n2
           it      = i
           thet    = ((n2-it)/(n2b))*(pi/2)
           ou(j,i) = in(j,i)*sin(thet)*sin(thet)
        end do
     end do
  end if

  !################################
  ! End taper
  !################################

  call srite('out',ou,size(ou)*4)
  call exit (0)

end program Taper
