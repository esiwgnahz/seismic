program Acoustic_rtm_opt_memory

  use sep
  use Readsouvelrho_mod
  use ExtractPadModel_mod
  use FDcoefs_assign
  use Propagator_mod
  use Interpolate_mod
  use Imaging_mod
  use Taper_mod

  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types

  implicit none

  type(GeneralParam) :: genpar
  type(ModelSpace)   :: mod
  type(DataSpace)    :: dat
  type(FDbounds)     :: bounds
  type(ModelSpace_elevation) :: elev

  type(WaveSpace), target                     :: wfld_fwd
  type(TraceSpace), dimension(:), allocatable :: datavec
  type(TraceSpace), dimension(:), allocatable :: sourcevec

  integer :: i,j,k,counting(4),count_rate,count_max
  integer :: ntsnap
  real :: d1,totcount(3)

  call sep_init()
  
  totcount=0.

  genpar%lsinc=7

  call from_param('fmax',genpar%fmax,30.)
  call from_param('ntaper',genpar%ntaper,20)
  call from_param('snapi',genpar%snapi,4)
  call from_param('aperture_x',genpar%aperture(1))
  call from_param('aperture_y',genpar%aperture(2))
  call from_param('num_threads',genpar%nthreads,4)

  write(0,*) 'INFO: -------- Parameters ---------'
  write(0,*) 'INFO:'
  write(0,*) 'INFO: fmax  =',genpar%fmax
  write(0,*) 'INFO: ntaper=',genpar%ntaper
  write(0,*) 'INFO: snapi =',genpar%snapi
  write(0,*) 'INFO: aperture_x=',genpar%aperture(1)
  write(0,*) 'INFO: aperture_y=',genpar%aperture(2)
  write(0,*) 'INFO:'
  write(0,*) 'INFO: num_threads=',genpar%nthreads
  write(0,*) 'INFO: ----------------------------'

  call omp_set_num_threads(genpar%nthreads)

  mod%veltag='vel'
  mod%rhotag='rho'
  genpar%Born=.false.

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
  genpar%tmin=1
  genpar%tmax=sourcevec(1)%dimt%nt
  genpar%tstep=1

  call system_clock(counting(1),count_rate,count_max)

  allocate(wfld_fwd%wave(mod%nz,mod%nxw,mod%nyw,genpar%ntsnap,1))
  
  write(0,*) 'INFO: Starting forward modeling'
  if (genpar%twoD) then
     if (.not.genpar%withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_2D_derivatives_scalar_forward_grid,     &
        & Injection_sinc,                                &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &    
        & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd) 
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_2D_derivatives_acoustic_forward_grid,       &
        & Injection_rho_sinc,                            &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &    
        & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd)
     end if
  else
     if (.not.genpar%withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_3D_derivatives_scalar_forward_grid,     &
        & Injection_sinc,                                &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &    
        & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd) 
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_3D_derivatives_acoustic_forward_grid,       &
        & Injection_rho_sinc,                            &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &    
        & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd)
     end if
  end if

  call system_clock(counting(2),count_rate,count_max)

  write(0,*) 'INFO: Done with forward modeling'
  
  allocate(mod%imagesmall(mod%nz,mod%nxw,mod%nyw))
  allocate(mod%illumsmall(mod%nz,mod%nxw,mod%nyw))
  mod%imagesmall=0.
  mod%illumsmall=0.

  genpar%tmax=1
  genpar%tmin=sourcevec(1)%dimt%nt
  genpar%tstep=-1

  mod%counter=0
  mod%wvfld=>wfld_fwd

  call compute_taper(mod)

  write(0,*) 'INFO: Starting backward propagation'
  if (genpar%twoD) then
     if (.not.genpar%withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_2D_derivatives_scalar_adjoint_grid,     &
        & Injection_sinc,                                &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &
        & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_memory)
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_2D_derivatives_acoustic_forward_grid,       &
        & Injection_rho_sinc,                            &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &
        & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_memory)
     end if
  else
     if (.not.genpar%withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_3D_derivatives_scalar_forward_grid,     &
        & Injection_sinc,                                &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &
        & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_memory)
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_3D_derivatives_acoustic_forward_grid,       &
        & Injection_rho_sinc,                            &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &
        & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_memory)
     end if
  end if
  write(0,*) 'INFO: Done with backward propagation'

  call to_history('n1',sourcevec(1)%dimt%nt,'modeled_data')
  call to_history('n2',size(datavec),'modeled_data')
  call to_history('d1',sourcevec(1)%dimt%dt,'modeled_data')
  call to_history('d2',1.,'modeled_data')
  call to_history('o1',0.,'modeled_data')
  call to_history('o2',0.,'modeled_data')

  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do
  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec)
  deallocate(sourcevec)

  call deallocateWaveSpace(wfld_fwd)
  
  call system_clock(counting(3),count_rate,count_max)

  call mod_copy_image_to_disk(mod)

  call system_clock(counting(4),count_rate,count_max)

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

  do i=1,3
     totcount(i)=totcount(i)+float(counting(i+1)-counting(i))/float(count_rate)
  end do

  write(0,*) 'INFO ------------------------------'
  write(0,*) 'INFO Total Wall-clock time       = ',sum(totcount)
  write(0,*) 'INFO ------------------------------'
  write(0,*) 'INFO  * Source propagation       = ',100*totcount(1)/sum(totcount),'%',totcount(1)
  write(0,*) 'INFO  * Receiver propagation     = ',100*totcount(2)/sum(totcount),'%',totcount(2)
  write(0,*) 'INFO  * Copy image/illum to disk = ',100*totcount(3)/sum(totcount),'%',totcount(3)
  write(0,*) 'INFO ------------------------------'

end program Acoustic_rtm_opt_memory
