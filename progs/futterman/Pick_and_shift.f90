! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Pick_and_shift

  use sep
  use picking
  use statistics

  use omp_lib
  implicit none

  real, dimension(:), allocatable :: trace,ER,RM,ou,pck,pcko,weight,tracein,traceou


  integer, dimension(:), allocatable :: n
  real,    dimension(:), allocatable :: d,o

  integer, dimension(:), allocatable :: nt_pick
  integer :: ndim,nl,nw,ntimes,i,j,k,pad2
  logical :: pick_wb,shift_forward,pad2_trace
  character(len=1024) :: label
  real :: pick_t0,tmax,tpick_max,dt_shift
  integer :: ntmax,n2dim

  call sep_init(SOURCE)
  ndim=sep_dimension()

  allocate(n(ndim),o(ndim),d(ndim))

  do i=1,ndim
     call sep_get_data_axis_par('in',i,n(i),o(i),d(i),label)
  end do

  if (ndim.eq.1) then
     n2dim=1
  else
     n2dim=product(n(2:))
  end if
  call from_param('pick_wb',pick_wb,.false.)

  allocate(pcko(n2dim))
  allocate(nt_pick(n2dim))

  call from_param('shift_forward',shift_forward,.true.)

  if (shift_forward) then

     allocate(trace(n(1)))
     write(0,*) 'INFO: Picking first breaks'
     call from_param('n_leading',nl,20)
     call from_param('n_smooth',nw,int(1.5*nl))
     call from_param('miss_iter',ntimes,1)

     allocate(pck(n2dim))

     allocate(ER(n(1)))
     allocate(RM(n(1)))
     allocate(ou(n(1)))

     do i=1,n2dim
        call sreed('in',trace,4*n(1))
        if (minval(trace).eq.maxval(trace)) then
           ou=0.
           pck(i)=1e19
           cycle
        end if
        trace=trace/maxval(abs(trace))
        call coppens(trace,ER,nl)
        call running_median(ER,RM,nw)
        call derivative(RM,ou)
        call find_maxindex(ou,nt_pick(i))
        pck(i)=o(1)+d(1)*nt_pick(i)
     end do

     pcko=pck
     ! Takes care of zero traces
     do i=2,n2dim
        if (pck(i).eq.1e19) then
           pcko(i)=pck(i-1)
        end if
     end do

     do j=1,ntimes
        do i=n2dim-1,1,-1
           if (abs(pcko(i)-pcko(i+1)).gt.0.1*pcko(i+1)) then
              pcko(i)=pcko(i+1)
           end if
        end do
     end do
     write(0,*) 'INFO: Done with picking first breaks'
 
     call sseek('in',0,0)
     call srite('pick',pcko,4*n2dim)
     call sep_put_data_axis_par('pick',1,n2dim,o(2),d(2),label)
     deallocate(trace)
     tmax=(n(1)-1)*d(1)+o(1)
     tpick_max=minval(pcko)
     ntmax=max(1,int((tmax-tpick_max)/d(1)+1.5)) 

     call from_param('pad2',pad2_trace,.true.)
     if (pad2_trace) ntmax=pad2(ntmax) ! make sure ntmax is a power of 2 for FFTs
     call auxclose('pick')
  else

     call from_aux('shifted','n1',ntmax)
     
  end if

  allocate(trace(ntmax))
  allocate(tracein(n(1)))
  allocate(traceou(ntmax))

  call from_param('dt_shift',dt_shift,0.15)

  if (shift_forward) then

     do i=1,n2dim
       
        tracein=0.

        if (modulo(i,max(1,int(10*n2dim/100))).eq.0) write(0,*) 'INFO: Shifting forward trace',i,'/',n2dim
        nt_pick(i)=max(1,int((pcko(i)-o(1)-dt_shift)/d(1)+1.5))
        call sreed('in',tracein,4*n(1))
        trace=0
        trace(1:n(1)-nt_pick(i)+1)=tracein(nt_pick(i):n(1))
        call srite('shifted',trace,4*ntmax)
        
     end do

     call sep_put_data_axis_par('shifted',1,ntmax,tpick_max,d(1),label)
     do i=2,ndim
        call sep_put_data_axis_par('shifted',i,n(i),o(i),d(i),label)
     end do
     write(0,*) 'INFO: Done with forward shifting'

  else 

     call sreed('pick',pcko,4*n2dim)

     do i=1,n2dim

        tracein=0.
        trace=0.
        traceou=0.

        if (modulo(i,max(1,int(10*n2dim/100))).eq.0) write(0,*) 'INFO: Shifting backward trace',i,'/',n2dim

        nt_pick(i)=max(1,int((pcko(i)-o(1)-dt_shift)/d(1)+1.5))
        call sreed('shifted',trace,4*ntmax)
        tracein(nt_pick(i):)=trace(1:n(1)-nt_pick(i)+1)
         
        call srite('out',tracein,4*n(1))   

     end do
     call sep_put_data_axis_par('out',1,n(1),o(1),d(1),label)

     do i=2,ndim
        call sep_put_data_axis_par('out',i,n(i),o(i),d(i),label)
     end do
     write(0,*) 'INFO: Done with backward shifting'
  end if
  
end program Pick_and_shift

integer function pad2(n)
  integer n
  pad2=1 
  do while (pad2.lt.n)
     pad2=pad2*2
  end do
  return
end function pad2
