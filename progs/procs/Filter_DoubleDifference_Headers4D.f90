! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program FilterHeadersDD4D

  use sep
  use DataSpace_types

  implicit none

  type(TraceSpace), dimension(:), allocatable :: datahead1 ! First data header
  type(TraceSpace), dimension(:), allocatable :: datahead2 ! Second data header
  real, dimension(:), allocatable             :: mask1     ! header masks for header 1. Use sfheaderwindow to take the traces out
  real, dimension(:), allocatable             :: mask2     ! header masks for header 2. Use sfheaderwindow to take the traces out

  type(TraceSpace) :: datasourcecoord1
  type(TraceSpace) :: datasourcecoord2

  integer :: index_gx,index_gy,index_sx,index_sy
  integer :: ntraces1,ntraces2,nkeys1,nkeys2
  real, allocatable, dimension(:) :: tracekeys1
  real, allocatable, dimension(:) :: tracekeys2

  integer :: ntsnap,i,j,ni,nj
  real    :: distances,distancer,maxd

  call sep_init()
  call FilterHeadersDD4D_doc()

  write(0,*) 'INFO:'
  write(0,*) 'INFO: -- Filter Double Difference Headers Starting -- '
  write(0,*) 'INFO:'

  write(0,*) 'INFO:------------------------'
  write(0,*) 'INFO: Starting reading traces'

  call from_aux('header1','n2',ntraces1);allocate(datahead1(ntraces1))
  call from_aux('header2','n2',ntraces2);allocate(datahead2(ntraces2))
  call from_aux('header1','n1',nkeys1);  allocate(tracekeys1(nkeys1))
  call from_aux('header2','n1',nkeys2);  allocate(tracekeys2(nkeys2))

  write(0,*) 'INFO: Header 1 has ',ntraces1,' traces with ',nkeys1,' headers'
  write(0,*) 'INFO: Header 2 has ',ntraces2,' traces with ',nkeys2,' headers'

  call from_param('keygx',index_gx)
  call from_param('keysx',index_sx)
  call from_param('keygy',index_gy)
  call from_param('keysy',index_sy)
  call from_param('max_distance',maxd)

  write(0,*) 'INFO: Reading header information gx, gy header1'
  do i=1,ntraces1
     tracekeys1=0.
     call sreed('header1',tracekeys1,4*nkeys1)   
     datahead1(i)%coord(2)=tracekeys1(index_gx)    
     datahead1(i)%coord(3)=tracekeys1(index_gy)
  end do

  write(0,*) 'INFO: Reading header information gx, gy header2'
  do i=1,ntraces2
     tracekeys2=0.
     call sreed('header2',tracekeys2,4*nkeys2)   
     datahead2(i)%coord(2)=tracekeys2(index_gx)    
     datahead2(i)%coord(3)=tracekeys2(index_gy)
  end do

  allocate(mask1(ntraces1)); mask1=0.
  allocate(mask2(ntraces2)); mask2=0.

  datasourcecoord1%coord(1)=tracekeys1(index_sx) 
  datasourcecoord1%coord(2)=tracekeys1(index_sy)
  datasourcecoord2%coord(1)=tracekeys2(index_sx) 
  datasourcecoord2%coord(2)=tracekeys2(index_sy)

  write(0,*) 'INFO: Source header 1 coordinates sx ',datasourcecoord1%coord(1),' sy ',datasourcecoord1%coord(2)
  write(0,*) 'INFO: Source header 2 coordinates sx ',datasourcecoord2%coord(1),' sy ',datasourcecoord2%coord(2)
   
  distances=sqrt((datasourcecoord1%coord(1)-datasourcecoord2%coord(1))**2+(datasourcecoord1%coord(2)-datasourcecoord2%coord(2))**2)
  
  write(0,*) 'INFO:      - Distance between two sources = ',distances

  do i=1,ntraces2
     do j=1,ntraces1
        distancer=sqrt((datahead1(j)%coord(2)-datahead2(i)%coord(2))**2+(datahead1(j)%coord(3)-datahead2(i)%coord(3))**2)
        if (distancer.le.maxd) then
           mask2(i)=1.
           exit
        end if
     end do
  end do

  do i=1,ntraces1
     do j=1,ntraces2
        distancer=sqrt((datahead1(i)%coord(2)-datahead2(j)%coord(2))**2+(datahead1(i)%coord(3)-datahead2(j)%coord(3))**2)
        if (distancer.le.maxd) then
           mask1(i)=1.
        end if
     end do
  end do

  write(0,*) 'INFO: checking number of traces left in header 1 ',sum(mask1)
  write(0,*) 'INFO: checking number of traces left in header 2 ',sum(mask2)

  call srite('mask1',mask1,4*ntraces1)
  call srite('mask2',mask2,4*ntraces2)
  call to_history('n1',ntraces1,'mask1')
  call to_history('n1',ntraces2,'mask2')

end program FILTERHEADERSDD4D

subroutine FILTERHEADERSDD4D_doc()
  
    call sep_add_doc_line("NAME")
    call sep_add_doc_line("    Filter_DoubleDifference_Headers4D.x - Creates two mask from two header files to keep interesection of trace coordinates ")
    call sep_add_doc_line("                                          that can be used for double differencing in 4D analysis")
    call sep_add_doc_line("                                          Assumes that the headers comes from the same shot position, more or less")
    call sep_add_doc_line("SYNOPSIS")
    call sep_add_doc_line("    Filter_DoubleDifference_Headers4D.x header1= header2= mask1= mask2= max_distance= keygx,y= keysx,y= > /dev/null ")
    call sep_add_doc_line("")
    call sep_add_doc_line("INPUT PARAMETERS")
    call sep_add_doc_line("")
    call sep_add_doc_line("    max_distance - float")
    call sep_add_doc_line("               Maximum distance allowed between two receivers from two headers")
    call sep_add_doc_line("")
    call sep_add_doc_line("    keysx - integer")
    call sep_add_doc_line("                Shot X key index in header file")
    call sep_add_doc_line("")
    call sep_add_doc_line("    keysy - integer")
    call sep_add_doc_line("                Shot Y key index in header file")
    call sep_add_doc_line("")
    call sep_add_doc_line("    keygx - integer")
    call sep_add_doc_line("                Receiver X key index in header file")
    call sep_add_doc_line("")
    call sep_add_doc_line("    keygy - integer")
    call sep_add_doc_line("                Receiver Y key index in header file")
    call sep_add_doc_line("")
    call sep_add_doc_line("")
    call sep_add_doc_line("INPUT FILES")
    call sep_add_doc_line("    header1 - sepfile")
    call sep_add_doc_line("             (nh1,ntraces1) = input")
    call sep_add_doc_line("             First header file with nh1 keys and ntraces1 traces")
    call sep_add_doc_line("")
    call sep_add_doc_line("    header2 - sepfile")
    call sep_add_doc_line("             (nh1,ntraces1) = input")
    call sep_add_doc_line("             First header file with nh2 keys and ntraces2 traces")
    call sep_add_doc_line("")
    call sep_add_doc_line("OUTPUT FILES")
    call sep_add_doc_line("    mask1 - sepfile")
    call sep_add_doc_line("            (ntraces1) = output")
    call sep_add_doc_line("            mask file with 0 if we don't keep and 1 if we keep trace")
    call sep_add_doc_line("")
    call sep_add_doc_line("    mask2 - sepfile")
    call sep_add_doc_line("            (ntraces2) = output")
    call sep_add_doc_line("            mask file with 0 if we don't keep and 1 if we keep trace")
    call sep_add_doc_line("")
    call sep_add_doc_line("DESCRIPTION")
    call sep_add_doc_line("    Creates two mask files that can be used to filter out traces that are not commong to both shots")
    call sep_add_doc_line("    This is useful for 4D analysis")
    call sep_add_doc_line("")
    call sep_add_doc_line("COMMENTS")
    call sep_add_doc_line("     (1) This code creates mask that can be used by sfheaderwindow to then remove unwanted traces.")
    call sep_add_doc_line("     (2) A conversion to int will be needed in that case, using sfdd < maskin.H type=int > maskou.H ")
    call sep_add_doc_line("")
    call doc('SOURCE')

end subroutine FILTERHEADERSDD4D_doc
