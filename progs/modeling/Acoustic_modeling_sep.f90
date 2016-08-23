program Acoustic_modeling_sep

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

  type(WaveSpace)                             :: wfld_fwd
  type(TraceSpace), dimension(:), allocatable :: datavec
  type(TraceSpace), dimension(:), allocatable :: sourcevec

  integer :: i,j,k 
  integer :: ntsnap

  call sep_init()
  
  genpar%lsinc=7

  call from_param('fmax',genpar%fmax,30.)
  call from_param('ntaper',genpar%ntaper,20)

  genpar%snapi=4

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

  call readsou(sourcevec,genpar)
 
  if (genpar%withRho) then
     genpar%coefpower=1
  else
     genpar%coefpower=2
  end if

  call readvel(mod,genpar,bounds)
  call readtraces(datavec,sourcevec,genpar)
  call readcoords(datavec,sourcevec,genpar)
  call are_traces_within_model(datavec,mod)
  call are_traces_within_model(sourcevec,mod)

  genpar%ntsnap=int(genpar%nt/genpar%snapi)

  allocate(wfld_fwd%wave(mod%nz,mod%nx,mod%ny,genpar%ntsnap,1))
  
  write(0,*) 'bounds%nmin1',bounds%nmin1,'bounds%nmax1',bounds%nmax1
  write(0,*) 'bounds%nmin2',bounds%nmin2,'bounds%nmax2',bounds%nmax2
  write(0,*) 'bounds%nmin3',bounds%nmin3,'bounds%nmax3',bounds%nmax3

  allocate(elev%elev(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  elev%elev=0.
  genpar%tmin=1
  genpar%tmax=sourcevec(1)%dimt%nt
  genpar%tstep=1

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
        & ExtractData=Extraction_array_sinc)
     else
        call propagator_acoustic(                          &
        & FD_acoustic_rho_init_coefs,                      &
        & FD_2D_derivatives_acoustic_forward_grid,         &
        & Injection_rho_sinc,                              &
        & FD_2nd_time_derivative_grid,                     &
        & FDswaptime_pointer,                              &
        & bounds,mod,elev,genpar,                          &
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,     &
        & ExtractData=Extraction_array_sinc)
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
        & ExtractData=Extraction_array_sinc) 
     else
        call propagator_acoustic(                          &
        & FD_acoustic_rho_init_coefs,                      &
        & FD_3D_derivatives_acoustic_forward_grid,         &
        & Injection_rho_sinc,                              &
        & FD_2nd_time_derivative_grid,                     &
        & FDswaptime_pointer,                              &
        & bounds,mod,elev,genpar,                          &
        & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,     &
        & ExtractData=Extraction_array_sinc)
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

  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do
  call deallocateWaveSpace(wfld_fwd)
  call deallocateTraceSpace(sourcevec(1))
  deallocate(datavec)
  deallocate(sourcevec)

end program Acoustic_modeling_sep
