! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program THREEDBORN

  use sep
  use Readsouvelrho_mod

  use to_disk_to_memory_BORN_mod

  use ReadParam_mod
  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types

  implicit none

  type(GeneralParam) :: genpar
  type(ModelSpace)   :: mod
  type(DataSpace)    :: dat
  type(FDbounds)     :: bounds
  type(ModelSpace_elevation) :: elev

  type(TraceSpace), dimension(:), allocatable :: datavec
  type(TraceSpace), dimension(:), allocatable :: sourcevec

  integer :: i,j,k 
  integer :: ntsnap
  double precision :: memory_needed

  call sep_init()

  write(0,*) 'INFO:'
  write(0,*) 'INFO: -- Born Modeling Starting -- '
  write(0,*) 'INFO:'

  mod%veltag='vel'
  mod%rhotag='rho'
  mod%reftag='ref'
  mod%waFtag='wave_fwd'

  if (.not.exist_file(mod%veltag)) call erexit('ERROR: need velocity file')
  if (.not.exist_file(mod%reftag)) call erexit('ERROR: need reflectivity file')
  if (.not.exist_file('data')) call erexit('ERROR: need dat file for output')

  call read_3D_params(genpar)
  genpar%Born=.true.
  call readsou(sourcevec,genpar)
  call readtraces(datavec,sourcevec,genpar)
  call readcoords(datavec,sourcevec,genpar)
  call extract_coord_source_receiver_patch(datavec,sourcevec(1),mod,genpar)
  call read_window_vel(mod,genpar,bounds)
  call read_window_ref(mod,genpar,bounds)

  genpar%ntsnap=int(genpar%nt/genpar%snapi)
  
  write(0,*) 'INFO: ------- Size model space for propagation --------'
  write(0,*) 'INFO:'
  write(0,*) 'INFO: bounds%nmin1',bounds%nmin1,'bounds%nmax1',bounds%nmax1
  write(0,*) 'INFO: bounds%nmin2',bounds%nmin2,'bounds%nmax2',bounds%nmax2
  write(0,*) 'INFO: bounds%nmin3',bounds%nmin3,'bounds%nmax3',bounds%nmax3
  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------------------'

  allocate(elev%elev(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  elev%elev=0.

  memory_needed=dble(mod%nxw)*dble(mod%nz)*dble(mod%nyw)*1e-9*4*dble(genpar%ntsnap)
  write(0,*) 'INFO:'
  write(0,*) 'INFO: Memory needed to write wavefield =',memory_needed,'Gb'

  call omp_set_num_threads(genpar%nthreads)

  if (dble(genpar%max_memory).gt.memory_needed) then
     write(0,*) 'INFO: writing wavefield in memory'
     write(0,*) 'INFO:'
     call BORN_to_memory(mod,genpar,dat,bounds,elev,datavec,sourcevec)
  else
     write(0,*) 'INFO: writing wavefield on disk file',mod%waFtag
     write(0,*) 'INFO:'
     call auxinout(mod%waFtag)
     call BORN_to_disk(mod,genpar,dat,bounds,elev,datavec,sourcevec)
     call auxclose(mod%waFtag)
  end if

  do i=1,size(datavec)
     call srite('data',datavec(i)%trace(:,1),4*sourcevec(1)%dimt%nt)
  end do

  call to_history('n1',sourcevec(1)%dimt%nt,'data')
  call to_history('n2',size(datavec),'data')
  call to_history('d1',sourcevec(1)%dimt%dt,'data')
  call to_history('d2',1.,'data')
  call to_history('o1',0.,'data')
  call to_history('o2',0.,'data')

  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do

  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec,sourcevec)
  call deallocateModelSpace(mod)

  write(0,*) 'INFO:'
  write(0,*) 'INFO: -- Born Modeling End -- '
  write(0,*) 'INFO:'

end program THREEDBORN
