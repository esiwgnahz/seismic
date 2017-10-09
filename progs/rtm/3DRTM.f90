! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program THREEDRTM

  use sep
  use Readsouvelrho_mod

  use to_disk_to_memory_RTM_mod

  use ReadParam_mod
  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types

  use THREED_Doc

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
  call RTM_doc()

  write(0,*) 'INFO:'
  write(0,*) 'INFO: -- RTM Starting -- '
  write(0,*) 'INFO:'

  mod%veltag='vel'
  mod%vel2tag='vel2'
  mod%rhotag='rho'
  mod%waFtag='wave_fwd'
  genpar%Born=.false.
  mod%exist_vel2=.false.

  if (.not.exist_file(mod%veltag)) call erexit('ERROR: need vel file')
  if (exist_file(mod%vel2tag)) then
     mod%exist_vel2=.true.
     write(0,*) 'INFO: ---- A second velocity file is used for bacward propagation ----'
     write(0,*) 'INFO: ---- This usually means that PS or SP waves are imaged      ----'
  end if

  call read_3D_params(genpar)
  call readsou(sourcevec,genpar)

  call readtraces(datavec,sourcevec,genpar)
  call readcoords(datavec,sourcevec,genpar)
  call extract_coord_source_receiver_patch(datavec,sourcevec(1),mod,genpar)
  call read_window_vel(mod,genpar,bounds)

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

  memory_needed=dble(mod%nxw)*dble(mod%nz)*dble(mod%nyw)*9.3132257e-10*4*dble(genpar%ntsnap)
  write(0,*) 'INFO:'
  write(0,*) 'INFO: Memory needed to write wavefield =',memory_needed,'Gb'

  call omp_set_num_threads(genpar%nthreads)

  if (dble(genpar%max_memory).gt.memory_needed) then
     write(0,*) 'INFO: writing wavefield in memory'
     write(0,*) 'INFO:'
     if (genpar%LSRTM) then
        if (genpar%CHEAPLSRTM) then
           call CHEAP_LSRTM_to_memory(mod,genpar,dat,bounds,elev,datavec,sourcevec)
        else
           call LSRTM_to_memory(mod,genpar,dat,bounds,elev,datavec,sourcevec)
        end if
     else
        call RTM_to_memory(mod,genpar,dat,bounds,elev,datavec,sourcevec)
     end if
  else
     write(0,*) 'INFO: writing wavefield on disk file',mod%waFtag
     write(0,*) 'INFO:'
     call auxinout(mod%waFtag)
     if (genpar%LSRTM) then
        if (genpar%CHEAPLSRTM) then
           call CHEAP_LSRTM_to_disk(mod,genpar,dat,bounds,elev,datavec,sourcevec)
        else
           call LSRTM_to_disk(mod,genpar,dat,bounds,elev,datavec,sourcevec)
        end if
     else
        call RTM_to_disk(mod,genpar,dat,bounds,elev,datavec,sourcevec)
     end if
     call auxclose(mod%waFtag)
  end if

  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do

  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec)
  deallocate(sourcevec)

  call mod_copy_image_to_disk(mod)

  call to_history('n1',mod%nz,'image')
  call to_history('n2',mod%nx,'image')
  call to_history('d1',mod%dz,'image')
  call to_history('d2',mod%dx,'image')
  call to_history('o1',mod%oz,'image')
  call to_history('o2',mod%ox,'image')
  if (genpar%twoD) then
     call to_history('n3',2,'image')
  else
     call to_history('n3',mod%ny,'image')
     call to_history('d3',mod%dy,'image')
     call to_history('o3',mod%oy,'image')
     call to_history('n4',2,'image')
  end if

  call deallocateModelSpace(mod)

  write(0,*) 'INFO:'
  write(0,*) 'INFO: -- RTM End -- '
  write(0,*) 'INFO:'

end program THREEDRTM
