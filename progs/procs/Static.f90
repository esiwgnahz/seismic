! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Static3D

  use sep
  
  use Readsouvelrho_mod
  use GeneralParam_types
  use DataSpace_types

  type(GeneralParam) :: genpar
  type(DataSpace)    :: dat

  type(TraceSpace), dimension(:), allocatable :: datavec
  type(TraceSpace), dimension(:), allocatable :: sourcevec
  real, dimension(:), allocatable :: trace

  real :: t0,dt,ot,t
  real :: tstat
  integer :: i,nt,it,ntstat,nsample
  logical :: verb

  call sep_init()

  allocate(sourcevec(1))
  call from_aux('traces','n1',sourcevec(1)%dimt%nt)
  call from_aux('traces','d1',sourcevec(1)%dimt%dt)
  call readtraces(datavec,sourcevec,genpar)
  call readtstat(datavec,sourcevec,genpar)

  nt=datavec(1)%dimt%nt
  ot=datavec(1)%dimt%ot
  dt=datavec(1)%dimt%dt

  call from_param('verb',verb,.true.)
  call from_param('num_threads',genpar%nthreads,4)
  
  call omp_set_num_threads(genpar%nthreads)

  if (verb) then
     write(0,*) 'INFO:'
     write(0,*) 'INFO: Simple static correction of 3D shot gathers using tstat'
     write(0,*) 'INFO:'
     write(0,*) 'INFO: Nthreads=',genpar%nthreads     
     write(0,*) 'INFO: Number of traces',size(datavec)
  end if

  !$OMP PARALLEL DO PRIVATE(i,tstat,ntstat,trace)
  do i=1,size(datavec)
     tstat=datavec(i)%tstat
     ntstat=tstat*0.001/dt
     !write(0,*) tstat,ntstat,nt
     allocate(trace(nt))
     trace=datavec(i)%trace(:,1)
     datavec(i)%trace(:,1)=0
     datavec(i)%trace(ntstat+1:,1)=trace(1:nt-ntstat)
     deallocate(trace)
  end do
  !$OMP END PARALLEL DO

  do i=1,size(datavec)
     call srite('out',datavec(i)%trace(:,1),4*nt)
  end do
  call to_history('n1',nt,'out')
  call to_history('n2',size(datavec),'out')
  call to_history('d1',dt,'out')
  call to_history('d2',1.,'out')
  call to_history('o1',ot,'out')
  call to_history('o2',0.,'out')

  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do

  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec,sourcevec)

  write(0,*) 'INFO:'
  write(0,*) 'INFO: -- Static correction Done -- '
  write(0,*) 'INFO:'
end program Static3D
