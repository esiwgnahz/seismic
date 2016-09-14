module to_disk_to_memory_RTM_mod

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

    subroutine RTM_to_disk(mod,genpar,dat,bounds,elev,datavec,sourcevec)

      type(GeneralParam) :: genpar
      type(ModelSpace)   :: mod
      type(DataSpace)    :: dat
      type(FDbounds)     :: bounds
      type(ModelSpace_elevation) :: elev

      type(TraceSpace), dimension(:), allocatable :: datavec
      type(TraceSpace), dimension(:), allocatable :: sourcevec
      
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

      genpar%tmin=1
      genpar%tmax=sourcevec(1)%dimt%nt
      genpar%tstep=1

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
            & sou=sourcevec,ExtractWave=Extraction_wavefield_copy_to_disk) 
         else
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_2D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                            &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &    
            & sou=sourcevec,ExtractWave=Extraction_wavefield_copy_to_disk)
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
            & sou=sourcevec,ExtractWave=Extraction_wavefield_copy_to_disk) 
         else
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_3D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                            &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &    
            & sou=sourcevec,ExtractWave=Extraction_wavefield_copy_to_disk)
         end if
      end if

      write(0,*) 'INFO: Done with forward modeling'

      allocate(mod%imagesmall(mod%nz,mod%nxw,mod%nyw))
      allocate(mod%illumsmall(mod%nz,mod%nxw,mod%nyw))
      mod%imagesmall=0.
      mod%illumsmall=0.

      genpar%tmax=1
      genpar%tmin=sourcevec(1)%dimt%nt
      genpar%tstep=-1

      mod%counter=0
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
            & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_disk)
         else
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_2D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                            &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &
            & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_disk)
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
            & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_disk)
         else
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_3D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                            &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &
            & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_disk)
         end if
      end if
      write(0,*) 'INFO: Done with backward propagation'

    end subroutine RTM_to_disk
    
    subroutine LSRTM_to_disk(mod,genpar,dat,bounds,elev,datavec,sourcevec)

      type(GeneralParam) :: genpar
      type(ModelSpace)   :: mod
      type(DataSpace)    :: dat
      type(FDbounds)     :: bounds
      type(ModelSpace_elevation) :: elev

      type(TraceSpace), dimension(:), allocatable :: datavec
      type(TraceSpace), dimension(:), allocatable :: sourcevec
      
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

      genpar%tmin=1
      genpar%tmax=sourcevec(1)%dimt%nt
      genpar%tstep=1

      write(0,*) 'INFO: Starting forward modeling'
      if (genpar%twoD) then
         call propagator_acoustic(                        &
         & FD_acoustic_init_coefs,                        &
         & FD_2nd_2D_derivatives_scalar_forward_grid,     &
         & Injection_LSRTM_sinc,                          &
         & FD_2nd_time_derivative_grid,                   &
         & FDswaptime_pointer,                            &
         & bounds,mod,elev,genpar,                        &    
         & sou=sourcevec,ExtractWave=Extraction_wavefield_copy_to_disk)         
      else
         call propagator_acoustic(                        &
         & FD_acoustic_init_coefs,                        &
         & FD_2nd_3D_derivatives_scalar_forward_grid,     &
         & Injection_LSRTM_sinc,                          &
         & FD_2nd_time_derivative_grid,                   &
         & FDswaptime_pointer,                            &
         & bounds,mod,elev,genpar,                        &    
         & sou=sourcevec,ExtractWave=Extraction_wavefield_copy_to_disk)         
      end if

      write(0,*) 'INFO: Done with forward modeling'

      allocate(mod%imagesmall(mod%nz,mod%nxw,mod%nyw))
      allocate(mod%illumsmall(mod%nz,mod%nxw,mod%nyw))
      mod%imagesmall=0.
      mod%illumsmall=0.

      genpar%tmax=1
      genpar%tmin=sourcevec(1)%dimt%nt
      genpar%tstep=-1

      mod%counter=0
      call compute_taper(mod)

      write(0,*) 'INFO: Starting backward propagation'
      if (genpar%twoD) then
         call propagator_acoustic(                        &
         & FD_acoustic_init_coefs,                        &
         & FD_2nd_2D_derivatives_scalar_adjoint_grid,     &
         & Injection_LSRTM_sinc,                          &
         & FD_2nd_time_derivative_grid,                   &
         & FDswaptime_pointer,                            &
         & bounds,mod,elev,genpar,                        &
         & sou=datavec,ImagingCondition=Imaging_condition_LSRTM_sourceonly_from_disk)        
      else
         call propagator_acoustic(                        &
         & FD_acoustic_init_coefs,                        &
         & FD_2nd_3D_derivatives_scalar_adjoint_grid,     &
         & Injection_LSRTM_sinc,                          &
         & FD_2nd_time_derivative_grid,                   &
         & FDswaptime_pointer,                            &
         & bounds,mod,elev,genpar,                        &
         & sou=datavec,ImagingCondition=Imaging_condition_LSRTM_sourceonly_from_disk)       
      end if
      write(0,*) 'INFO: Done with backward propagation'

    end subroutine LSRTM_to_disk
    
    subroutine RTM_to_memory(mod,genpar,dat,bounds,elev,datavec,sourcevec)

      type(GeneralParam) :: genpar
      type(ModelSpace)   :: mod
      type(DataSpace)    :: dat
      type(FDbounds)     :: bounds
      type(ModelSpace_elevation) :: elev

      type(TraceSpace), dimension(:), allocatable :: datavec
      type(TraceSpace), dimension(:), allocatable :: sourcevec
      type(WaveSpace), target                     :: wfld_fwd
      
      genpar%tmin=1
      genpar%tmax=sourcevec(1)%dimt%nt
      genpar%tstep=1

      allocate(wfld_fwd%wave(mod%nz,mod%nxw,mod%nyw,genpar%ntsnap,1))

      wfld_fwd%wave=0.
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
            & Injection_sinc,                            &
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
            & Injection_sinc,                            &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &    
            & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd)
         end if
      end if

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
            & Injection_sinc,                            &
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
            & Injection_sinc,                            &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &
            & sou=datavec,ImagingCondition=Imaging_condition_sourceonly_from_memory)
         end if
      end if
      write(0,*) 'INFO: Done with backward propagation'

      call deallocateWaveSpace(wfld_fwd)
    end subroutine RTM_to_memory
    
    subroutine LSRTM_to_memory(mod,genpar,dat,bounds,elev,datavec,sourcevec)

      type(GeneralParam) :: genpar
      type(ModelSpace)   :: mod
      type(DataSpace)    :: dat
      type(FDbounds)     :: bounds
      type(ModelSpace_elevation) :: elev

      type(TraceSpace), dimension(:), allocatable :: datavec
      type(TraceSpace), dimension(:), allocatable :: sourcevec
      type(WaveSpace), target                     :: wfld_fwd

      genpar%tmin=1
      genpar%tmax=sourcevec(1)%dimt%nt
      genpar%tstep=1

      allocate(wfld_fwd%wave(mod%nz,mod%nxw,mod%nyw,genpar%ntsnap,1))

      wfld_fwd%wave=0.
      write(0,*) 'INFO: Starting forward modeling'
      if (genpar%twoD) then
         call propagator_acoustic(                        &
         & FD_acoustic_init_coefs,                        &
         & FD_2nd_2D_derivatives_scalar_forward_grid,     &
         & Injection_LSRTM_sinc,                          &
         & FD_2nd_time_derivative_grid,                   &
         & FDswaptime_pointer,                            &
         & bounds,mod,elev,genpar,                        &    
         & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd) 

      else
         call propagator_acoustic(                        &
         & FD_acoustic_init_coefs,                        &
         & FD_2nd_3D_derivatives_scalar_forward_grid,     &
         & Injection_LSRTM_sinc,                          &
         & FD_2nd_time_derivative_grid,                   &
         & FDswaptime_pointer,                            &
         & bounds,mod,elev,genpar,                        &    
         & sou=sourcevec,ExtractWave=Extraction_wavefield,wfld=wfld_fwd) 

      end if

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
         call propagator_acoustic(                        &
         & FD_acoustic_init_coefs,                        &
         & FD_2nd_2D_derivatives_scalar_adjoint_grid,     &
         & Injection_LSRTM_sinc,                          &
         & FD_2nd_time_derivative_grid,                   &
         & FDswaptime_pointer,                            &
         & bounds,mod,elev,genpar,                        &
         & sou=datavec,ImagingCondition=Imaging_condition_LSRTM_sourceonly_from_memory)
      else
         call propagator_acoustic(                        &
         & FD_acoustic_init_coefs,                        &
         & FD_2nd_3D_derivatives_scalar_adjoint_grid,     &
         & Injection_LSRTM_sinc,                          &
         & FD_2nd_time_derivative_grid,                   &
         & FDswaptime_pointer,                            &
         & bounds,mod,elev,genpar,                        &
         & sou=datavec,ImagingCondition=Imaging_condition_LSRTM_sourceonly_from_memory)
      end if
      write(0,*) 'INFO: Done with backward propagation'

      call deallocateWaveSpace(wfld_fwd)
    end subroutine LSRTM_to_memory
    
  end module to_disk_to_memory_RTM_mod
