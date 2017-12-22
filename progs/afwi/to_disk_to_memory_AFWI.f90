! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module to_disk_to_memory_AFWI_mod

  use sep

  use Mute_gather
  use Readsouvelrho_mod
  use ExtractPadModel_mod
  use FDcoefs_assign
  use FD_derivatives
  use Propagator_mod
  use Interpolate_mod
  use Injection_mod
  use Imaging_mod
  use Taper_mod

  use Inversion_mod
  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types

  use OF_Res_AdjSrc_mod

  implicit none

  contains

    subroutine AFWI_to_disk(mod,invparam,mutepar,genpar,dat,bounds,elev,datavec,sourcevec)

      type(MuteParam)     :: mutepar
      type(InversionParam):: invparam
      type(GeneralParam)  :: genpar
      type(ModelSpace)    :: mod
      type(DataSpace)     :: dat
      type(FDbounds)      :: bounds
      type(ModelSpace_elevation) :: elev

      type(TraceSpace), dimension(:), allocatable :: datavec
      type(TraceSpace), dimension(:), allocatable :: sourcevec

      type(TraceSpace), dimension(:), allocatable :: dmodvec

      integer :: ntraces,i
      integer(kind=8)  :: ntotal
      double precision :: f
      
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

      ntotal=size(datavec)*sourcevec(1)%dimt%nt
      ntraces=size(datavec)
      allocate(dmodvec(ntraces))
      do i=1,ntraces
         allocate(dmodvec(i)%trace(datavec(i)%dimt%nt,1))
         dmodvec(i)=datavec(i)
         dmodvec(i)%trace=0.
      end do

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
            & sou=sourcevec,datavec=dmodvec,                 &
            & ExtractData=Extraction_array_sinc_afwi,        &
            & ExtractWave=Extraction_wavefield_copy_to_disk) 
         else
            call allocate_sxx2d_szz2d_delp2d(bounds,mod)
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_2D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &    
            & sou=sourcevec,datavec=dmodvec,                 &
            & ExtractData=Extraction_array_sinc_afwi,        &
            & ExtractWave=Extraction_wavefield_copy_to_disk)
            call deallocate_sxx2d_szz2d_delp2d()
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
            & sou=sourcevec,datavec=dmodvec,                 &
            & ExtractData=Extraction_array_sinc_afwi,        &
            & ExtractWave=Extraction_wavefield_copy_to_disk) 
         else
            call allocate_sxx3d_syy3d_szz3d_delp3d(bounds,mod)
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_3D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &    
            & sou=sourcevec,datavec=dmodvec,                 &
            & ExtractData=Extraction_array_sinc_afwi,        &
            & ExtractWave=Extraction_wavefield_copy_to_disk)
            call deallocate_sxx3d_syy3d_szz3d_delp3d()
         end if
      end if

      write(0,*) 'INFO: Done with forward modeling'

      if (exist_file('dmod')) then
         do i=1,ntraces
            call srite('dmod',mutepar%maskgath(1)%gathtrace(i)%trace*dmodvec(i)%trace,4*dmodvec(i)%dimt%nt)
         end do
         call to_history('n1',dmodvec(1)%dimt%nt,'dmod')
         call to_history('n2',ntraces,'dmod')
      end if

      if (exist_file('residual')) then
         do i=1,ntraces
            call srite('residual',mutepar%maskgath(1)%gathtrace(i)%trace*(datavec(i)%trace-dmodvec(i)%trace),4*datavec(i)%dimt%nt)
         end do
         call to_history('n1',datavec(1)%dimt%nt,'residual')
         call to_history('n2',ntraces,'residual')
      end if

      call Compute_OF_RES_3D(invparam,datavec,dmodvec,mutepar%maskgath(1)%gathtrace,f)
      call to_history('n1',1,'function')
      call srite('function',sngl(f),4)

      if (invparam%nparam.eq.1) then
         allocate(mod%imagesmall(mod%nz,mod%nxw,mod%nyw))
         mod%imagesmall=0.
      else
         allocate(mod%imagesmall_nparam(mod%nz,mod%nxw,mod%nyw,2))
         mod%imagesmall_nparam=0.
      end if
      allocate(mod%illumsmall(mod%nz,mod%nxw,mod%nyw))
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
            & sou=dmodvec,ImagingCondition=Imaging_condition_AFWI_from_disk)
         else
            call allocate_sxx2d_szz2d_delp2d(bounds,mod)
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_2D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &
            & sou=dmodvec,ImagingCondition=Imaging_condition_AFWI_RHOVP_from_disk)
            call deallocate_sxx2d_szz2d_delp2d()
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
            & sou=dmodvec,ImagingCondition=Imaging_condition_AFWI_from_disk)
         else
            call allocate_sxx3d_syy3d_szz3d_delp3d(bounds,mod)
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_3D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &
            & sou=dmodvec,ImagingCondition=Imaging_condition_AFWI_RHOVP_from_disk)
            call deallocate_sxx3d_syy3d_szz3d_delp3d()
         end if
      end if
      write(0,*) 'INFO: Done with backward propagation'

      do i=1,size(dmodvec)
         call deallocateTraceSpace(dmodvec(i))
      end do
      deallocate(dmodvec)

    end subroutine AFWI_to_disk

    subroutine AFWI_to_memory(mod,invparam,mutepar,genpar,dat,bounds,elev,datavec,sourcevec)

      type(MuteParam)     :: mutepar
      type(InversionParam):: invparam
      type(GeneralParam)  :: genpar
      type(ModelSpace)    :: mod
      type(DataSpace)     :: dat
      type(FDbounds)      :: bounds
      type(ModelSpace_elevation) :: elev

      type(TraceSpace), dimension(:), allocatable :: datavec
      type(TraceSpace), dimension(:), allocatable :: sourcevec
      type(WaveSpace), target                     :: wfld_fwd
      
      type(TraceSpace), dimension(:), allocatable :: dmodvec

      integer :: ntraces,i
      double precision :: f
      integer(kind=8) :: ntotal

      genpar%tmin=1
      genpar%tmax=sourcevec(1)%dimt%nt
      genpar%tstep=1

      ntotal=size(datavec)*sourcevec(1)%dimt%nt
      ntraces=size(datavec)
      allocate(dmodvec(ntraces))
      do i=1,ntraces
         allocate(dmodvec(i)%trace(datavec(i)%dimt%nt,1))
         dmodvec(i)=datavec(i)
         dmodvec(i)%trace=0.
      end do
      
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
            & sou=sourcevec,datavec=dmodvec,                 &
            & ExtractData=Extraction_array_sinc_afwi,        &
            & ExtractWave=Extraction_wavefield,wfld=wfld_fwd) 
         else
            call allocate_sxx2d_szz2d_delp2d(bounds,mod)
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_2D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &    
            & sou=sourcevec,datavec=dmodvec,                 &
            & ExtractData=Extraction_array_sinc_afwi,        &
            & ExtractWave=Extraction_wavefield,wfld=wfld_fwd)
            call deallocate_sxx2d_szz2d_delp2d()
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
            & sou=sourcevec,datavec=dmodvec,                 &
            & ExtractData=Extraction_array_sinc_afwi,        &
            & ExtractWave=Extraction_wavefield,wfld=wfld_fwd) 
         else
            call allocate_sxx3d_syy3d_szz3d_delp3d(bounds,mod)
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_3D_derivatives_acoustic_forward_grid,      &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &    
            & sou=sourcevec,datavec=dmodvec,                 &
            & ExtractData=Extraction_array_sinc_afwi,        &
            & ExtractWave=Extraction_wavefield,wfld=wfld_fwd)
            call deallocate_sxx3d_syy3d_szz3d_delp3d()
         end if
      end if

      write(0,*) 'INFO: Done with forward modeling'

      if (exist_file('dmod')) then
         do i=1,ntraces
            call srite('dmod',mutepar%maskgath(1)%gathtrace(i)%trace*dmodvec(i)%trace,4*dmodvec(i)%dimt%nt)
         end do
         call to_history('n1',dmodvec(1)%dimt%nt,'dmod')
         call to_history('n2',ntraces,'dmod')
      end if

      if (exist_file('residual')) then
         do i=1,ntraces
            call srite('residual',mutepar%maskgath(1)%gathtrace(i)%trace*(datavec(i)%trace-dmodvec(i)%trace),4*datavec(i)%dimt%nt)
         end do
         call to_history('n1',datavec(1)%dimt%nt,'residual')
         call to_history('n2',ntraces,'residual')
      end if

      call Compute_OF_RES_3D(invparam,datavec,dmodvec,mutepar%maskgath(1)%gathtrace,f)
      call to_history('n1',1,'function')
      call srite('function',sngl(f),4)

      if (invparam%nparam.eq.1) then
         allocate(mod%imagesmall(mod%nz,mod%nxw,mod%nyw))
         mod%imagesmall=0.
      else
         allocate(mod%imagesmall_nparam(mod%nz,mod%nxw,mod%nyw,2))
         mod%imagesmall_nparam=0.
      end if
      allocate(mod%illumsmall(mod%nz,mod%nxw,mod%nyw))
      mod%illumsmall=0.

      genpar%tmax=1
      genpar%tmin=sourcevec(1)%dimt%nt
      genpar%tstep=-1

      mod%counter=0
      mod%wvfld=>wfld_fwd

      call compute_taper(mod)

      write(0,*) 'INFO: Starting backward propagation',genpar%withRho
      if (genpar%twoD) then
         if (.not.genpar%withRho) then
            call propagator_acoustic(                        &
            & FD_acoustic_init_coefs,                        &
            & FD_2nd_2D_derivatives_scalar_adjoint_grid,     &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &
            & sou=dmodvec,ImagingCondition=Imaging_condition_AFWI_from_memory)
         else
            call allocate_sxx2d_szz2d_delp2d(bounds,mod)
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_2D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &
            & sou=dmodvec,ImagingCondition=Imaging_condition_AFWI_RHOVP_from_memory)
            call deallocate_sxx2d_szz2d_delp2d()
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
            & sou=dmodvec,ImagingCondition=Imaging_condition_AFWI_from_memory)
         else
            call allocate_sxx3d_syy3d_szz3d_delp3d(bounds,mod)
            call propagator_acoustic(                        &
            & FD_acoustic_rho_init_coefs,                    &
            & FD_3D_derivatives_acoustic_forward_grid,       &
            & Injection_sinc,                                &
            & FD_2nd_time_derivative_grid,                   &
            & FDswaptime_pointer,                            &
            & bounds,mod,elev,genpar,                        &
            & sou=dmodvec,ImagingCondition=Imaging_condition_AFWI_RHOVP_from_memory)
            call deallocate_sxx3d_syy3d_szz3d_delp3d()
         end if
      end if
      write(0,*) 'INFO: Done with backward propagation'

      do i=1,size(dmodvec)
         call deallocateTraceSpace(dmodvec(i))
      end do
      deallocate(dmodvec)

      call deallocateWaveSpace(wfld_fwd)
    end subroutine AFWI_to_memory
    
  end module to_disk_to_memory_AFWI_mod
