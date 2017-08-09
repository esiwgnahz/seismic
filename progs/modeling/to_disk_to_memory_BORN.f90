! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module to_disk_to_memory_BORN_mod

  use sep

  use Readsouvelrho_mod
  use ExtractPadModel_mod
  use FDcoefs_assign
  use Propagator_mod
  use Interpolate_mod
  use Injection_mod
  use Imaging_mod
  use Taper_mod

  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types

  implicit none

contains

  subroutine BORN_to_disk(model,genpar,dat,bounds,elev,datavec,sourcevec)

    type(GeneralParam) :: genpar
    type(ModelSpace)   :: model
    type(DataSpace)    :: dat
    type(FDbounds)     :: bounds
    type(ModelSpace_elevation) :: elev

    type(TraceSpace), dimension(:), allocatable :: datavec
    type(TraceSpace), dimension(:), allocatable :: sourcevec

    integer :: i,j,k 
    integer :: ntsnap

    call to_history('n1',model%nz,'wave_fwd')
    call to_history('n2',model%nxw,'wave_fwd')
    call to_history('d1',model%dz,'wave_fwd')
    call to_history('d2',model%dx,'wave_fwd')
    call to_history('o1',genpar%omodel(1),'wave_fwd')
    call to_history('o2',genpar%omodel(2),'wave_fwd')
    if (genpar%twoD) then
       call to_history('n3',genpar%ntsnap,'wave_fwd')
    else
       call to_history('n3',model%nyw,'wave_fwd')
       call to_history('d3',model%dy,'wave_fwd')
       call to_history('o3',genpar%omodel(3),'wave_fwd')
       call to_history('n4',genpar%ntsnap,'wave_fwd')
    end if

    genpar%tmin=1
    genpar%tmax=sourcevec(1)%dimt%nt
    genpar%tstep=1

    write(0,*) 'INFO: Starting forward modeling'
    if (genpar%twoD) then
       call propagator_acoustic(                        &
       & FD_acoustic_init_coefs,                        &
       & FD_2nd_2D_derivatives_scalar_forward_grid,     &
       & Injection_sinc,                                &
       & FD_2nd_time_derivative_grid,                   &
       & FDswaptime_pointer,                            &
       & bounds,model,elev,genpar,                      &
       & sou=sourcevec,ExtractWave=Extraction_wavefield_copy_to_disk)
    else
       call propagator_acoustic(                        &
       & FD_acoustic_init_coefs,                        &
       & FD_2nd_3D_derivatives_scalar_forward_grid,     &
       & Injection_sinc,                                &
       & FD_2nd_time_derivative_grid,                   &
       & FDswaptime_pointer,                            &
       & bounds,model,elev,genpar,                      &
       & sou=sourcevec,ExtractWave=Extraction_wavefield_copy_to_disk)
    end if

    write(0,*) 'INFO: Done with forward modeling'

    model%counter=0
    call compute_taper(model)

    ! Set datavec to zero: modeling!
    do i=1,size(datavec)
       datavec(i)%trace(:,1)=0.
    end do

    call sseek(model%waFtag,0,0)

    write(0,*) 'INFO: Starting backward propagation'
    if (genpar%twoD) then
       call propagator_acoustic(                          &
       & FD_acoustic_init_coefs,                          &
       & FD_2nd_2D_derivatives_scalar_forward_grid,       &
       & Injection_LSRTM_sinc,                            &
       & FD_2nd_time_derivative_grid,                     &
       & FDswaptime_pointer,                              &
       & bounds,model,elev,genpar,                        &
       & ImagingCondition=Injection_Born_from_disk,       &
       & datavec=datavec,ExtractData=Extraction_array_LSRTM_sinc)
    else
       call propagator_acoustic(                          &
       & FD_acoustic_init_coefs,                          &
       & FD_2nd_3D_derivatives_scalar_forward_grid,       &
       & Injection_LSRTM_sinc,                            &
       & FD_2nd_time_derivative_grid,                     &
       & FDswaptime_pointer,                              &
       & bounds,model,elev,genpar,                        &
       & ImagingCondition=Injection_Born_from_disk,       &
       & datavec=datavec,ExtractData=Extraction_array_LSRTM_sinc)    
    end if
    write(0,*) 'INFO: Done with backward propagation'

  end subroutine BORN_to_disk

  subroutine BORN_to_memory(model,genpar,dat,bounds,elev,datavec,sourcevec)

    type(GeneralParam) :: genpar
    type(ModelSpace)   :: model
    type(DataSpace)    :: dat
    type(FDbounds)     :: bounds
    type(ModelSpace_elevation) :: elev

    type(TraceSpace), dimension(:), allocatable :: datavec
    type(TraceSpace), dimension(:), allocatable :: sourcevec
    type(WaveSpace), target                     :: wfld_fwd

    integer :: i,j,k 
    integer :: ntsnap

    genpar%tmin=1
    genpar%tmax=sourcevec(1)%dimt%nt
    genpar%tstep=1

    allocate(wfld_fwd%wave(model%nz,model%nxw,model%nyw,genpar%ntsnap,1))
    wfld_fwd%wave=0.
    write(0,*) 'INFO: Starting forward modeling'
    if (genpar%twoD) then
       call propagator_acoustic(                        &
       & FD_acoustic_init_coefs,                        &
       & FD_2nd_2D_derivatives_scalar_forward_grid,     &
       & Injection_sinc,                                &
       & FD_2nd_time_derivative_grid,                   &
       & FDswaptime_pointer,                            &
       & bounds,model,elev,genpar,                      &
       & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd)
    else
       call propagator_acoustic(                        &
       & FD_acoustic_init_coefs,                        &
       & FD_2nd_3D_derivatives_scalar_forward_grid,     &
       & Injection_sinc,                                &
       & FD_2nd_time_derivative_grid,                   &
       & FDswaptime_pointer,                            &
       & bounds,model,elev,genpar,                      &
       & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd)
    end if
    write(0,*) 'INFO: Done with forward modeling'

    model%wvfld=>wfld_fwd
    model%counter=0

    ! Set datavec to zero: modeling!
    do i=1,size(datavec)
       datavec(i)%trace(:,1)=0.
    end do

    write(0,*) 'INFO:'
    if (genpar%WriteFwdWvfld) then
       if (genpar%twoD) then
          do i=1,genpar%ntsnap
             if (mod(i,20).eq.0) write(0,*) 'INFO: Writing forward wavefield',i,'/',genpar%ntsnap         
             call srite('wave_fwd',model%wvfld%wave(1:model%nz,1:model%nxw,1:model%nyw,i,1),4*model%nxw*model%nyw*model%nz)
          end do
       else
          do i=1,int(genpar%ntsnap/8),genpar%ntsnap
             write(0,*) 'INFO: Writing forward wavefield',i,'/',genpar%ntsnap         
             call srite('wave_fwd',model%wvfld%wave(1:model%nz,1:model%nxw,1:model%nyw,i,1),4*model%nxw*model%nyw*model%nz)
          end do
       end if
    end if
    write(0,*) 'INFO:'

    call compute_taper(model)

    write(0,*) 'INFO: Starting 2nd forward propagation'
    if (genpar%twoD) then
       call propagator_acoustic(                          &
       & FD_acoustic_init_coefs,                          &
       & FD_2nd_2D_derivatives_scalar_forward_grid,       &
       & Injection_LSRTM_sinc,                            &
       & FD_2nd_time_derivative_grid,                     &
       & FDswaptime_pointer,                              &
       & bounds,model,elev,genpar,                        &
       & ImagingCondition=Injection_Born,                 &
       & datavec=datavec,ExtractData=Extraction_array_LSRTM_sinc)
    else
       call propagator_acoustic(                          &
       & FD_acoustic_init_coefs,                          &
       & FD_2nd_3D_derivatives_scalar_forward_grid,       &
       & Injection_LSRTM_sinc,                            &
       & FD_2nd_time_derivative_grid,                     &
       & FDswaptime_pointer,                              &
       & bounds,model,elev,genpar,                        &
       & ImagingCondition=Injection_Born,                 &
       & datavec=datavec,ExtractData=Extraction_array_LSRTM_sinc)    
    end if
    write(0,*) 'INFO: Done with 2nd forward propagation'
    call deallocateWaveSpace(wfld_fwd)

  end subroutine BORN_to_memory

end module to_disk_to_memory_BORN_mod
