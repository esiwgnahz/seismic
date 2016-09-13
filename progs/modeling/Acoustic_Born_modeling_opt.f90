program Acoustic_Born_modeling_opt

  use sep
  use Readsouvelrho_mod
  use ExtractPadModel_mod
  use FDcoefs_assign
  use Propagator_mod
  use Interpolate_mod

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

  integer :: i,j,k ,counting(4),count_rate,count_max
  integer :: ntsnap,totcount(3)

  call sep_init()

  totcount=0.

  call from_param('fmax',genpar%fmax,30.)
  call from_param('ntaper',genpar%ntaper,20)
  call from_param('snapi',genpar%snapi,4)
  call from_param('aperture_x',genpar%aperture(1))
  call from_param('aperture_y',genpar%aperture(2))

  call from_param('num_threads',genpar%nthreads,4)
  call from_param('lsinc',genpar%lsinc,7)

  write(0,*) 'INFO: -------- Parameters ---------'
  write(0,*) 'INFO:'
  write(0,*) 'INFO: fmax  =',genpar%fmax
  write(0,*) 'INFO: ntaper=',genpar%ntaper
  write(0,*) 'INFO: snapi =',genpar%snapi
  write(0,*) 'INFO: aperture_x=',genpar%aperture(1)
  write(0,*) 'INFO: aperture_y=',genpar%aperture(2)
  write(0,*) 'INFO:'
  write(0,*) 'INFO: num_threads=',genpar%nthreads
  write(0,*) 'INFO: lsinc      =',genpar%lsinc
  write(0,*) 'INFO: ----------------------------'

  call omp_set_num_threads(genpar%nthreads)

  genpar%Born=.true.
  mod%veltag='vel'
  mod%rhotag='rho'
  mod%reftag='ref'
  mod%waFtag='wave_fwd'

  call auxinout(mod%waFtag)
  call from_param('twoD',genpar%twoD,.false.)

  if (.not.genpar%twoD) then     
     genpar%nbound=4
  else
     genpar%nbound=0
  end if

  call from_param('rec_type',genpar%rec_type,0)
  call from_param('shot_type',genpar%shot_type,0)
  call from_param('surf_type',genpar%surf_type,0)

  write(0,*) 'INFO: -------- Source/Receiver parameters ------'
  write(0,*) 'INFO:'
  write(0,*) 'INFO: rec_type =',genpar%rec_type
  write(0,*) 'INFO: shot_type=',genpar%shot_type
  write(0,*) 'INFO: surf_type=',genpar%surf_type
  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------------'

  call readsou(sourcevec,genpar)
 
  if (genpar%withRho) then
     genpar%coefpower=1
  else
     genpar%coefpower=2
  end if

  call readtraces(datavec,sourcevec,genpar)
  call readcoords(datavec,sourcevec,genpar)
  call extract_coord_source_receiver_patch(datavec,sourcevec,mod,genpar)
  call read_window_vel(mod,genpar,bounds)
  call read_window_ref(mod,genpar,bounds)

  genpar%ntsnap=int(genpar%nt/genpar%snapi)

  allocate(wfld_fwd%wave(mod%nz,mod%nxw,mod%nyw,genpar%ntsnap,1))
  
  write(0,*) 'INFO: ------- Size model space for propagation --------'
  write(0,*) 'INFO:'
  write(0,*) 'INFO: bounds%nmin1',bounds%nmin1,'bounds%nmax1',bounds%nmax1
  write(0,*) 'INFO: bounds%nmin2',bounds%nmin2,'bounds%nmax2',bounds%nmax2
  write(0,*) 'INFO: bounds%nmin3',bounds%nmin3,'bounds%nmax3',bounds%nmax3
  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------------------'

  allocate(elev%elev(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  elev%elev=0.
  genpar%tmin=1
  genpar%tmax=sourcevec(1)%dimt%nt
  genpar%tstep=1

  call to_history('n1',mod%nz,'wave_fwd')
  call to_history('n2',mod%nxw,'wave_fwd')
  call to_history('d1',mod%dz,'wave_fwd')
  call to_history('d2',mod%dx,'wave_fwd')
  call to_history('o1',genpar%omodel(1),'wave_fwd')
  call to_history('o2',genpar%omodel(2),'wave_fwd')
  if (genpar%twoD) then
     call to_history('n3',genpar%ntsnap,'wave_fwd')
  else
     call to_history('n3',mod%nyw,'wave_fwd')
     call to_history('d3',mod%dy,'wave_fwd')
     call to_history('o3',genpar%omodel(3),'wave_fwd')
     call to_history('n4',genpar%ntsnap,'wave_fwd')
  end if

  call system_clock(counting(1),count_rate,count_max)

  write(0,*) 'INFO: Starting forward modeling'
  if (genpar%twoD) then
     call propagator_acoustic(                        &
     & FD_acoustic_init_coefs,                        &
     & FD_2nd_2D_derivatives_scalar_forward_grid,     &
     & Injection_sinc,                                &
     & FD_2nd_time_derivative_grid,                   &
     & FDswaptime_pointer,                            &
     & bounds,mod,elev,genpar,                        &
     & sou=sourcevec,wfld=wfld_fwd,ExtractWave=Extraction_wavefield_copy_to_disk)
  else
     call propagator_acoustic(                        &
     & FD_acoustic_init_coefs,                        &
     & FD_2nd_3D_derivatives_scalar_forward_grid,     &
     & Injection_sinc,                                &
     & FD_2nd_time_derivative_grid,                   &
     & FDswaptime_pointer,                            &
     & bounds,mod,elev,genpar,                        &
     & sou=sourcevec,wfld=wfld_fwd,ExtractWave=Extraction_wavefield_copy_to_disk)
  end if

  call system_clock(counting(2),count_rate,count_max)

  write(0,*) 'INFO: Done with forward modeling'
  
  mod%wvfld=>wfld_fwd
  mod%wvfld%counter=0

!  do i=1,genpar%ntsnap
!     call srite('wave_fwd',mod%wvfld%wave(1:mod%nz,1:mod%nxw,1:mod%nyw,i,1),4*mod%nxw*mod%nyw*mod%nz)
!  end do

  ! Set datavec to zero: modeling!
  do i=1,size(datavec)
     datavec(i)%trace(:,1)=0.
  end do

  call sseek(mod%waFtag,0,0)

  write(0,*) 'INFO: Starting backward propagation'
  if (genpar%twoD) then
     call propagator_acoustic(                          &
     & FD_acoustic_init_coefs,                          &
     & FD_2nd_2D_derivatives_scalar_forward_grid,       &
     & Injection_Born_from_disk,                        &
     & FD_2nd_time_derivative_grid,                     &
     & FDswaptime_pointer,                              &
     & bounds,mod,elev,genpar,                          &
     & datavec=datavec,ExtractData=Extraction_array_sinc)
  else
     call propagator_acoustic(                          &
     & FD_acoustic_init_coefs,                          &
     & FD_2nd_3D_derivatives_scalar_forward_grid,       &
     & Injection_Born_from_disk,                        &
     & FD_2nd_time_derivative_grid,                     &
     & FDswaptime_pointer,                              &
     & bounds,mod,elev,genpar,                          &
     & datavec=datavec,ExtractData=Extraction_array_sinc)    
  end if
  write(0,*) 'INFO: Done with backward propagation'

  call system_clock(counting(3),count_rate,count_max)

  do i=1,size(datavec)
     call srite('data',datavec(i)%trace(:,1),4*sourcevec(1)%dimt%nt)
  end do

  call to_history('n1',sourcevec(1)%dimt%nt,'data')
  call to_history('n2',size(datavec),'data')
  call to_history('d1',sourcevec(1)%dimt%dt,'data')
  call to_history('d2',1.,'data')
  call to_history('o1',0.,'data')
  call to_history('o2',0.,'data')

  call to_history('n1',mod%nz,'wave_fwd')
  call to_history('n2',mod%nxw,'wave_fwd')
  call to_history('d1',mod%dz,'wave_fwd')
  call to_history('d2',mod%dx,'wave_fwd')
  call to_history('o1',genpar%omodel(1),'wave_fwd')
  call to_history('o2',genpar%omodel(2),'wave_fwd')
  if (genpar%twoD) then
     call to_history('n3',genpar%ntsnap,'wave_fwd')
  else
     call to_history('n3',mod%nyw,'wave_fwd')
     call to_history('d3',mod%dy,'wave_fwd')
     call to_history('o3',genpar%omodel(3),'wave_fwd')
     call to_history('n4',genpar%ntsnap,'wave_fwd')
  end if

  call auxclose(mod%waFtag)
  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do
  call deallocateWaveSpace(wfld_fwd)
  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec,sourcevec)

  call system_clock(counting(4),count_rate,count_max)

  do i=1,3
     totcount(i)=totcount(i)+float(counting(i+1)-counting(i))/float(count_rate)
  end do

  write(0,*) 'INFO ------------------------------'
  write(0,*) 'INFO Total Wall-clock time       = ',sum(totcount)
  write(0,*) 'INFO ------------------------------'
  write(0,*) 'INFO  * Source propagation       = ',100*totcount(1)/sum(totcount),'%',totcount(1)
  write(0,*) 'INFO  * Receiver propagation     = ',100*totcount(2)/sum(totcount),'%',totcount(2)
  write(0,*) 'INFO  * Copy to disk             = ',100*totcount(3)/sum(totcount),'%',totcount(3)
  write(0,*) 'INFO ------------------------------'

end program Acoustic_Born_modeling_opt
