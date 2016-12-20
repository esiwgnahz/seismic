program TaperEdges

  implicit none
  integer			                           	:: n1,n2,n3,ntap,i,j,k
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

  call getch("ntap","i",ntap);write(0,*)'INFO: from par ntap=',ntap
  if (ntap >= int(n2/2)) call erexit ('Incorrect ntap value, too big')

  !################################
  ! memory allocation/reading input
  !################################

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

  allocate (in (n1,n2))
  allocate (ou (n1,n2))
  do k=1,n3
     call sreed('in',in,size(in)*4)

     !################################
     ! taper the data with cos^2
     !################################

     ou = in

     do j=1,n1
        do i=1,ntap
           it      = i
           thet    = ((it-1)/(ntap-1))*(pi/2)
           ou(j,i) = in(j,i)*sin(thet)*sin(thet)
        end do
     end do
     !     
     do j=1,n1
        do i=n2-ntap,n2
           it      = i
           thet    = ((n2-it)/(ntap))*(pi/2)
           ou(j,i) = in(j,i)*sin(thet)*sin(thet)
        end do
     end do

     !################################
     ! End taper
     !################################

     call srite('out',ou,size(ou)*4)
  end do
  call exit (0)

end program TaperEdges
