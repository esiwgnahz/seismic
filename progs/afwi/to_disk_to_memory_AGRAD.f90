module to_disk_to_memory_AGRAD_mod

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

  subroutine AGRAD_to_memory(model,genpar,dat,bounds,elev,datavec,wfld_fwd,nparam)

    type(GeneralParam) :: genpar
    type(ModelSpace)   :: model
    type(DataSpace)    :: dat
    type(FDbounds)     :: bounds
    type(ModelSpace_elevation) :: elev

    type(TraceSpace), dimension(:) :: datavec
    type(WaveSpace), target        :: wfld_fwd

    integer :: i,j,k,ishot,nparam
    integer :: ntsnap

    genpar%tmax=1
    genpar%tmin=datavec(1)%dimt%nt
    genpar%tstep=-1

    allocate(model%imagesmall_nparam(model%nz,model%nxw,model%nyw,nparam))
    allocate(model%illumsmall(model%nz,model%nxw,model%nyw))
    model%imagesmall_nparam=0.
    model%illumsmall=0.

    model%wvfld=>wfld_fwd
    model%counter=0
   
    call compute_taper(model)
    if (genpar%verbose) write(0,*) 'INFO: Starting 2nd forward propagation'
    if (genpar%twoD) then
       if (genpar%withRho) then
          call propagator_acoustic(                          &
          & FD_acoustic_rho_init_coefs,                      &
          & FD_2D_derivatives_acoustic_forward_grid,         &
          & Injection_simple_noomp,                          &
          & FD_2nd_time_derivative_grid_noomp,               &
          & FDswaptime_pointer,                              &
          & bounds,model,elev,genpar,                        &
          & sou=datavec,ImagingCondition=Imaging_condition_AFWI_RHOVP_from_memory_noomp2)
       else
          call propagator_acoustic(                          &
          & FD_acoustic_init_coefs,                          &
          & FD_2nd_2D_derivatives_scalar_forward_grid_noomp, &
          & Injection_simple_noomp,                          &
          & FD_2nd_time_derivative_grid_noomp,               &
          & FDswaptime_pointer,                              &
          & bounds,model,elev,genpar,                        &
          & sou=datavec,ImagingCondition=Imaging_condition_AFWI_RHOVP_from_memory_noomp2)
       end if
    else
       if (genpar%withRho) then
          call propagator_acoustic(                          &
          & FD_acoustic_rho_init_coefs,                      &
          & FD_3D_derivatives_acoustic_forward_grid_noomp,   &
          & Injection_simple_noomp,                          &
          & FD_2nd_time_derivative_grid_noomp,               &
          & FDswaptime_pointer,                              &
          & bounds,model,elev,genpar,                        &
          & sou=datavec,ImagingCondition=Imaging_condition_AFWI_RHOVP_from_memory_noomp2)
       else
          call propagator_acoustic(                          &
          & FD_acoustic_init_coefs,                          &
          & FD_2nd_3D_derivatives_scalar_forward_grid_noomp, &
          & Injection_simple_noomp,                           &
          & FD_2nd_time_derivative_grid_noomp,               &
          & FDswaptime_pointer,                              &
          & bounds,model,elev,genpar,                        &
          & sou=datavec,ImagingCondition=Imaging_condition_AFWI_RHOVP_from_memory_noomp2)
       end if
    end if
    if (genpar%verbose) write(0,*) 'INFO: Done with 2nd forward propagation'
    call deallocateWaveSpace(wfld_fwd)
    call compute_taper_close(model)

  end subroutine AGRAD_to_memory

end module to_disk_to_memory_AGRAD_mod
