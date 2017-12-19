! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Acoustic_modeling_sep

  use sep
  use Readsouvelrho_mod
  use ExtractPadModel_mod
  use FDcoefs_assign
  use Propagator_mod
  use Interpolate_mod


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

  type(WaveSpace)                             :: wfld_fwd
  type(TraceSpace), dimension(:), allocatable :: datavec
  type(TraceSpace), dimension(:), allocatable :: sourcevec

  integer :: i,j,k 
  integer :: ntsnap
  logical :: i_want_wavefield

  call sep_init()

  mod%veltag='vel'
  mod%rhotag='rho'
  genpar%Born=.false.

  call from_param('i_want_wavefield',i_want_wavefield,.false.)
  call read_3D_params(genpar)
  call readsou(sourcevec,genpar)

  call readtraces(datavec,sourcevec,genpar)
  call readcoords(datavec,sourcevec,genpar)
  call extract_coord_source_receiver_patch(datavec,sourcevec(1),mod,genpar)
  call read_window_vel(mod,genpar,bounds)

  genpar%ntsnap=int(genpar%nt/genpar%snapi)

  allocate(wfld_fwd%wave(mod%nz,mod%nxw,mod%nyw,genpar%ntsnap,1))
  
  write(0,*) 'INFO: ------- Size model space for propagation --------'
  write(0,*) 'INFO:'
  write(0,*) 'INFO: bounds%nmin1',bounds%nmin1,'bounds%nmax1',bounds%nmax1
  write(0,*) 'INFO: bounds%nmin2',bounds%nmin2,'bounds%nmax2',bounds%nmax2
  write(0,*) 'INFO: bounds%nmin3',bounds%nmin3,'bounds%nmax3',bounds%nmax3
  write(0,*) 'INFO:'
  write(0,*) 'INFO: ------------------------------------------------'

  call omp_set_num_threads(genpar%nthreads)

  allocate(elev%elev(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  elev%elev=0.
  genpar%tmin=1
  genpar%tmax=sourcevec(1)%dimt%nt
  genpar%tstep=1

  ! Set datavec to zero: modeling!
  do i=1,size(datavec)
     datavec(i)%trace(:,1)=0.
  end do

  write(0,*) 'before wave propagator'
  if (genpar%twoD) then
     if (.not.genpar%withRho) then
        call propagator_acoustic(                          &
        & FD_acoustic_init_coefs,                          &
        & FD_2nd_2D_derivatives_scalar_forward_grid,       &
        & Injection_sinc,                                  &
        & FD_2nd_time_derivative_grid,                     & 
        & FDswaptime_pointer,                              &
        & bounds,mod,elev,genpar,                          &
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,     &
        & ExtractData=Extraction_array_sinc,ExtractWave=Extraction_wavefield)
     else
        call propagator_acoustic(                          &
        & FD_acoustic_rho_init_coefs,                      &
        & FD_2D_derivatives_acoustic_forward_grid,         &
        & Injection_sinc,                              &
        & FD_2nd_time_derivative_grid,                     &
        & FDswaptime_pointer,                              &
        & bounds,mod,elev,genpar,                          &
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,     &
        & ExtractData=Extraction_array_sinc,ExtractWave=Extraction_wavefield)
     end if
  else
     if (.not.genpar%withRho) then
        call propagator_acoustic(                          &
        & FD_acoustic_init_coefs,                          &
        & FD_2nd_3D_derivatives_scalar_forward_grid,       &
        & Injection_sinc,                                  &
        & FD_2nd_time_derivative_grid,                     &
        & FDswaptime_pointer,                              &
        & bounds,mod,elev,genpar,                          &
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,     &
        & ExtractData=Extraction_array_sinc,ExtractWave=Extraction_wavefield) 
     else
        call propagator_acoustic(                          &
        & FD_acoustic_rho_init_coefs,                      &
        & FD_3D_derivatives_acoustic_forward_grid,         &
        & Injection_sinc,                              &
        & FD_2nd_time_derivative_grid,                     &
        & FDswaptime_pointer,                              &
        & bounds,mod,elev,genpar,                          &
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,     &
        & ExtractData=Extraction_array_sinc,ExtractWave=Extraction_wavefield)
     end if
  end if

  write(0,*) 'after wave propagator'
  
  do i=1,size(datavec)
     call srite('data',datavec(i)%trace(:,1),4*sourcevec(1)%dimt%nt)
  end do


  call to_history('n1',sourcevec(1)%dimt%nt,'data')
  call to_history('n2',size(datavec),'data')
  call to_history('d1',sourcevec(1)%dimt%dt,'data')
  call to_history('d2',1.,'data')
  call to_history('o1',0.,'data')
  call to_history('o2',0.,'data')

  if (i_want_wavefield) then

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

     do i=1,genpar%ntsnap
        call srite('wave_fwd',wfld_fwd%wave(1:mod%nz,1:mod%nxw,1:mod%nyw,i,1),4*mod%nxw*mod%nyw*mod%nz)
     end do
     
  end if

  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do
  call deallocateWaveSpace(wfld_fwd)
  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec)
  deallocate(sourcevec)

end program Acoustic_modeling_sep
