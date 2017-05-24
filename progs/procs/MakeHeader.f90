program MakeHeaders

  use sep

  implicit none

  integer :: nsx,nsy,nhx,nhy,ngx,ngy
  real    :: dsx,dsy,dhx,dhy,dgx,dgy
  real    :: osx,osy,ohx,ohy,ogx,ogy

  real    :: gx,gy,sx,sy,hx,hy

  integer :: counter,i,j,k,l

  real, dimension(:,:), allocatable:: ou
  character(len=3) :: stype

  call sep_init(SOURCE)

  call from_param('survey_type',stype,'sea')

  call from_param('nsx',nsx,1)
  call from_param('nsy',nsy,1)
  call from_param('osx',osx,0.)
  call from_param('osy',osy,0.)
  call from_param('dsx',dsx,10.)
  call from_param('dsy',dsy,10.)

  if (stype.eq.'lan') then
     write(0,*) 'LAND ACQUISITION'

     call from_param('ngx',ngx,1)
     call from_param('ngy',ngy,1)
     call from_param('ogx',ogx,0.)
     call from_param('ogy',ogy,0.)
     call from_param('dgx',dgx,10.)
     call from_param('dgy',dgy,10.)

     allocate(ou(4,nsx*nsy*ngx*ngy))

     counter=0
     do i=1,nsy     
        sy=(i-1)*dsy+osy
        do j=1,nsx
           sx=(j-1)*dsx+osx
           do k=1,ngy                  
              gy=(k-1)*dgy+ogy        
              do l=1,ngx        
                 counter=counter+1              
                 gx=(l-1)*dgx+ogx
                 ou(1,counter)=sx
                 ou(2,counter)=gx
                 ou(3,counter)=sy
                 ou(4,counter)=gy
              end do
           end do
        end do
     end do

  else if (stype.eq.'sea') then
     write(0,*) 'STREAMER ACQUISITION'

     call from_param('nhx',nhx,1)
     call from_param('nhy',nhy,1)
     call from_param('ohx',ohx,0.)
     call from_param('ohy',ohy,0.)
     call from_param('dhx',dhx,10.)
     call from_param('dhy',dhy,10.)
     
     allocate(ou(4,nsx*nsy*nhx*nhy))
     
     counter=0
     do i=1,nsy     
        sy=(i-1)*dsy+osy
        do j=1,nsx
           sx=(j-1)*dsx+osx
           do k=1,nhy          
              hy=(k-1)*dhy+ohy         
              gy=hy+sy          
              do l=1,nhx              
                 counter=counter+1
                 hx=(l-1)*dhx+ohx
                 
                 gx=hx+sx
                 ou(1,counter)=sx
                 ou(2,counter)=gx
                 ou(3,counter)=sy
                 ou(4,counter)=gy
                 
              end do
           end do
        end do
     end do
     
  end if
  call srite('header',ou,4*4*counter)
  call to_history('n1',4,'header')
  call to_history('n2',counter,'header')
  deallocate(ou)

end program MakeHeaders
