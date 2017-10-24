! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program THREEDAFWI

  use sep
  use Readsouvelrho_mod

  use to_disk_to_memory_AFWI_mod

  use Mute_gather

  use ReadParam_mod
  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types
  use Inversion_types

  use THREED_Doc

  implicit none

  type(GeneralParam) :: genpar
  type(ModelSpace)   :: mod
  type(DataSpace)    :: dat
  type(FDbounds)     :: bounds
  type(ModelSpace_elevation) :: elev

  type(TraceSpace), dimension(:), allocatable :: datavec
  type(TraceSpace), dimension(:), allocatable :: sourcevec

  type(InversionParam)                        :: invparam
  type(MuteParam)                             :: mutepar

  integer :: i,j,k
  integer :: ntsnap
  double precision :: memory_needed

  call sep_init()
  call RTM_doc()

  write(0,*) 'INFO:'
  write(0,*) 'INFO: -- 3DAFWI for one shot Starting -- '
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
     write(0,*) 'INFO: ---- A second velocity file is used for backward propagation ----'
     write(0,*) 'INFO: ---- This usually means that PS or SP waves are imaged      ----'
  end if

  call read_3D_params(genpar)
  call readsou(sourcevec,genpar)

  call readtraces(datavec,sourcevec,genpar)
  call readcoords(datavec,sourcevec,genpar)
  call extract_coord_source_receiver_patch(datavec,sourcevec(1),mod,genpar)
  call read_window_vel(mod,genpar,bounds)

  invparam%nparam=1
  if (genpar%withRho) invparam%nparam=2 
  call from_param('data_nrm_type',invparam%dat_nrm_type_char,'L2norm')
  call from_param('data_threshold',invparam%dat_thresh,0.)
  call from_param('vprho_param',invparam%vprho_param,0)

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

  call Init_MuteParam(mutepar)
  call MuteParam_compute_mask_1shot(mutepar,datavec,sourcevec)

  if (dble(genpar%max_memory).gt.memory_needed) then
     write(0,*) 'INFO: writing wavefield in memory'
     write(0,*) 'INFO:'
     call FWI_to_memory(mod,invparam,mutepar,genpar,dat,bounds,elev,datavec,sourcevec)
  else
     write(0,*) 'INFO: writing wavefield on disk file',mod%waFtag
     write(0,*) 'INFO:'
     call auxinout(mod%waFtag)
     call FWI_to_disk(mod,invparam,mutepar,genpar,dat,bounds,elev,datavec,sourcevec)
     call auxclose(mod%waFtag)
  end if

  call MuteParam_deallocate(mutepar)
  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do

  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec)
  deallocate(sourcevec)

  call mod_copy_gradient_to_disk(mod,invparam%nparam,invparam%vprho_param)

  call to_history('n1',mod%nz,'gradient')
  call to_history('n2',mod%nx,'gradient')
  call to_history('d1',mod%dz,'gradient')
  call to_history('d2',mod%dx,'gradient')
  call to_history('o1',mod%oz,'gradient')
  call to_history('o2',mod%ox,'gradient')
  if (genpar%twoD) then
     call to_history('n3',invparam%nparam+1,'gradient')
  else
     call to_history('n3',mod%ny,'gradient')
     call to_history('d3',mod%dy,'gradient')
     call to_history('o3',mod%oy,'gradient')
     call to_history('n4',invparam%nparam+1,'gradient')
  end if

  call deallocateModelSpace(mod)

  write(0,*) 'INFO:'
  write(0,*) 'INFO: -- 3DAFWI one shot End -- '
  write(0,*) 'INFO:'

end program THREEDAFWI
