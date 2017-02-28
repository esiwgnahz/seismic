module to_disk_to_memory_AMOD_mod

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

  subroutine AMOD_to_memory(model,genpar,dat,bounds,elev,datavec,sourcevec,wfld_fwd,ishot)

    type(GeneralParam) :: genpar
    type(ModelSpace)   :: model
    type(DataSpace)    :: dat
    type(FDbounds)     :: bounds
    type(ModelSpace_elevation) :: elev

    type(TraceSpace), dimension(:) :: datavec
    type(TraceSpace), dimension(:) :: sourcevec
    type(TraceSpace), dimension(:), allocatable :: source
    type(WaveSpace), target        :: wfld_fwd

    integer :: i,j,k,ishot
    integer :: ntsnap

    genpar%tmin=1
    genpar%tmax=sourcevec(ishot)%dimt%nt
    genpar%tstep=1

    allocate(source(1))
    allocate(source(1)%trace(sourcevec(ishot)%dimt%nt,1))
    source(1)%trace=0.
    source(1)=sourcevec(ishot)

    allocate(wfld_fwd%wave(model%nz,model%nxw,model%nyw,genpar%ntsnap,1))
    wfld_fwd%wave=0.

    if (genpar%verbose) write(0,*) 'INFO: Starting forward modeling'
    if (genpar%twoD) then
       if (genpar%withRho) then
          call propagator_acoustic(                          &
          & FD_acoustic_rho_init_coefs,                      &
          & FD_2D_derivatives_acoustic_forward_grid,         &
          & Injection_sinc_noomp,                            &
          & FD_2nd_time_derivative_grid_noomp,               &
          & FDswaptime_pointer,                              &
          & bounds,model,elev,genpar,                        &
          & sou=source,wfld=wfld_fwd,datavec=datavec,     &
          & ExtractData=Extraction_array_simple_afwi_noomp,ExtractWave=Extraction_wavefield)
       else
          call propagator_acoustic(                          &
          & FD_acoustic_init_coefs,                          &
          & FD_2nd_2D_derivatives_scalar_forward_grid_noomp, &
          & Injection_sinc_noomp,                            &
          & FD_2nd_time_derivative_grid_noomp,               &
          & FDswaptime_pointer,                              &
          & bounds,model,elev,genpar,                        &
          & sou=source,wfld=wfld_fwd,datavec=datavec,     &
          & ExtractData=Extraction_array_simple_afwi_noomp,ExtractWave=Extraction_wavefield)
       end if
    else
       if (genpar%withRho) then
          call propagator_acoustic(                          &
          & FD_acoustic_init_coefs,                          &
          & FD_3D_derivatives_acoustic_forward_grid_noomp,   &
          & Injection_sinc_noomp,                            &
          & FD_2nd_time_derivative_grid_noomp,               &
          & FDswaptime_pointer,                              &
          & bounds,model,elev,genpar,                        &
          & sou=source,wfld=wfld_fwd,datavec=datavec,     &
          & ExtractData=Extraction_array_simple_afwi_noomp,ExtractWave=Extraction_wavefield)
       else
          call propagator_acoustic(                          &
          & FD_acoustic_init_coefs,                          &
          & FD_2nd_3D_derivatives_scalar_forward_grid_noomp, &
          & Injection_sinc_noomp,                            &
          & FD_2nd_time_derivative_grid_noomp,               &
          & FDswaptime_pointer,                              &
          & bounds,model,elev,genpar,                        &
          & sou=source,wfld=wfld_fwd,datavec=datavec,     &
          & ExtractData=Extraction_array_simple_afwi_noomp,ExtractWave=Extraction_wavefield)
       end if
    end if
    if (genpar%verbose) write(0,*) 'INFO: Done with forward modeling'

    deallocate(source(1)%trace)
    deallocate(source)

  end subroutine AMOD_to_memory

end module to_disk_to_memory_AMOD_mod
