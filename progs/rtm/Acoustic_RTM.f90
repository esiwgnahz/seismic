! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Acoustic_rtm

  use sep
  use Readsouvelrho_mod
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

  type(WaveSpace)                             :: wfld_fwd,wfld_bwd
  type(TraceSpace), dimension(:), allocatable :: datavec
  type(TraceSpace), dimension(:), allocatable :: sourcevec

  integer :: i,j,k 
  integer :: ntsnap

  call sep_init()

  call from_param('fmax',genpar%fmax,30.)
  call from_param('ntaper',genpar%ntaper,20)
  call from_param('snapi',genpar%snapi,4)
  call from_param('lsinc',genpar%lsinc,7)

  call from_param('num_threads',genpar%nthreads,4)

  call omp_set_num_threads(genpar%nthreads)

  mod%veltag='vel'
  mod%rhotag='rho'

  call from_param('twoD',genpar%twoD,.false.)

  if (.not.genpar%twoD) then     
     genpar%nbound=4
  else
     genpar%nbound=0
  end if

  call from_param('rec_type',genpar%rec_type,0)
  call from_param('shot_type',genpar%shot_type,0)
  call from_param('surf_type',genpar%surf_type,0)

  call readsou(sourcevec,genpar)
 
  if (genpar%withRho) then
     genpar%coefpower=1
  else
     genpar%coefpower=2
  end if

  call readvel(mod,genpar,bounds)
  call readsoucoord(sourcevec,mod) 

  if (genpar%twoD) then
     allocate(datavec(mod%nx))
  else
     allocate(datavec(mod%nx+mod%ny))
  end if

  ! Acquisition on 2 perpendicular lines
  do i=1,mod%nx
     allocate(datavec(i)%trace(sourcevec(1)%dimt%nt,1)) ! 1 component trace    
     datavec(i)%trace=0.
     datavec(i)%coord(1)=genpar%omodel(1)  ! Z
     datavec(i)%coord(2)=genpar%omodel(2)+(i-1)*mod%dx ! X
     datavec(i)%coord(3)=(mod%ny-1)*mod%dy/2+genpar%omodel(3)
     datavec(i)%dimt%nt=sourcevec(1)%dimt%nt
  end do
  if (.not.genpar%twoD) then
     do i=1,mod%ny
        allocate(datavec(i+mod%nx)%trace(sourcevec(1)%dimt%nt,1)) ! 1 component trace    
        datavec(i+mod%nx)%trace=0.
        datavec(i+mod%nx)%coord(1)=genpar%omodel(1)  ! Z
        datavec(i+mod%nx)%coord(3)=genpar%omodel(3)+(i-1)*mod%dy ! X
        datavec(i+mod%nx)%coord(2)=(mod%nx-1)*mod%dx/2+genpar%omodel(2)
        datavec(i)%dimt%nt=sourcevec(1)%dimt%nt
     end do
  end if

  genpar%ntsnap=int(genpar%nt/genpar%snapi)

  allocate(wfld_fwd%wave(mod%nz,mod%nx,mod%ny,genpar%ntsnap,1))
  
  write(0,*) 'bounds%nmin1',bounds%nmin1,'bounds%nmax1',bounds%nmax1
  write(0,*) 'bounds%nmin2',bounds%nmin2,'bounds%nmax2',bounds%nmax2
  write(0,*) 'bounds%nmin3',bounds%nmin3,'bounds%nmax3',bounds%nmax3

  allocate(elev%elev(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  elev%elev=0.
  write(0,*) 'before wave propagator'
  genpar%tmin=1
  genpar%tmax=sourcevec(1)%dimt%nt
  genpar%tstep=1
  if (genpar%twoD) then
     if (.not.genpar%withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_2D_derivatives_scalar_forward_grid,     &
        & Injection_sinc,                                &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &    
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,   &
        & ExtractData=Extraction_array_sinc              ) 
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_2D_derivatives_acoustic_forward_grid,       &
        & Injection_sinc,                            &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &    
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,   &
        & ExtractData=Extraction_array_sinc              )
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
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,   &
        & ExtractData=Extraction_array_sinc              ) 
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_3D_derivatives_acoustic_forward_grid,       &
        & Injection_sinc,                            &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &    
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,   &
        & ExtractData=Extraction_array_sinc              )
     end if
  end if

  write(0,*) 'after wave propagator'
  
  do i=1,size(datavec)
     call srite('data',datavec(i)%trace(:,1),4*sourcevec(1)%dimt%nt)
  end do

  do i=1,genpar%ntsnap
     call srite('wave_fwd',wfld_fwd%wave(1:mod%nz,1:mod%nx,1:mod%ny,i,1),4*mod%nx*mod%ny*mod%nz)
  end do

  call deallocateWaveSpace(wfld_fwd)
  allocate(wfld_bwd%wave(mod%nz,mod%nx,mod%ny,genpar%ntsnap,1))

  genpar%tmax=1
  genpar%tmin=sourcevec(1)%dimt%nt
  genpar%tstep=-1
  if (genpar%twoD) then
     if (.not.genpar%withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_2D_derivatives_scalar_adjoint_grid,     &
        & Injection_sinc,                                &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &
        & sou=datavec,wfld=wfld_bwd)
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_2D_derivatives_acoustic_forward_grid,       &
        & Injection_sinc,                            &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &
        & sou=datavec,wfld=wfld_bwd)
     end if
  else
     if (.not.genpar%withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_3D_derivatives_scalar_adjoint_grid,     &
        & Injection_sinc,                                &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &
        & sou=datavec,wfld=wfld_bwd)
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_3D_derivatives_acoustic_forward_grid,       &
        & Injection_sinc,                            &
        & FD_2nd_time_derivative_grid,                   &
        & FDswaptime_pointer,                            &
        & bounds,mod,elev,genpar,                        &
        & sou=datavec,wfld=wfld_bwd)
     end if
  end if
  write(0,*) 'after wave propagator'

  do i=genpar%ntsnap,1,-1
     call srite('wave_bwd',wfld_bwd%wave(1:mod%nz,1:mod%nx,1:mod%ny,i,1),4*mod%nx*mod%ny*mod%nz)
  end do

  call to_history('n1',sourcevec(1)%dimt%nt,'data')
  call to_history('n2',size(datavec),'data')
  call to_history('d1',sourcevec(1)%dimt%dt,'data')
  call to_history('d2',1.,'data')
  call to_history('o1',0.,'data')
  call to_history('o2',0.,'data')

  call to_history('n1',mod%nz,'wave_bwd')
  call to_history('n2',mod%nx,'wave_bwd')
  call to_history('d1',mod%dz,'wave_bwd')
  call to_history('d2',mod%dx,'wave_bwd')
  call to_history('o1',genpar%omodel(1),'wave_bwd')
  call to_history('o2',genpar%omodel(2),'wave_bwd')
  if (genpar%twoD) then
     call to_history('n3',genpar%ntsnap,'wave_bwd')
  else
     call to_history('n3',mod%ny,'wave_bwd')
     call to_history('d3',mod%dy,'wave_bwd')
     call to_history('o3',genpar%omodel(3),'wave_bwd')
     call to_history('n4',genpar%ntsnap,'wave_bwd')
  end if

  call to_history('n1',mod%nz,'wave_fwd')
  call to_history('n2',mod%nx,'wave_fwd')
  call to_history('d1',mod%dz,'wave_fwd')
  call to_history('d2',mod%dx,'wave_fwd')
  call to_history('o1',genpar%omodel(1),'wave_fwd')
  call to_history('o2',genpar%omodel(2),'wave_fwd')
  if (genpar%twoD) then
     call to_history('n3',genpar%ntsnap,'wave_fwd')
  else
     call to_history('n3',mod%ny,'wave_fwd')
     call to_history('d3',mod%dy,'wave_fwd')
     call to_history('o3',genpar%omodel(3),'wave_fwd')
     call to_history('n4',genpar%ntsnap,'wave_fwd')
  end if


  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do
  call deallocateWaveSpace(wfld_bwd)
  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec)

end program Acoustic_rtm
