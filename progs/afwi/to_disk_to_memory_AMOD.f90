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

  subroutine AMOD_to_memory(model,genpar,dat,bounds,elev,datavec,sourcevec,wfld_fwd)

    type(GeneralParam) :: genpar
    type(ModelSpace)   :: model
    type(DataSpace)    :: dat
    type(FDbounds)     :: bounds
    type(ModelSpace_elevation) :: elev

    type(TraceSpace), dimension(:) :: datavec
    type(TraceSpace), dimension(1) :: sourcevec
    type(WaveSpace), target        :: wfld_fwd

    integer :: i,j,k,ishot
    integer :: ntsnap

    genpar%tmin=1
    genpar%tmax=sourcevec(1)%dimt%nt
    genpar%tstep=1

    allocate(wfld_fwd%wave(model%nz,model%nxw,model%nyw,genpar%ntsnap,1))
    wfld_fwd%wave=0.

    if (genpar%verbose) write(0,*) 'INFO: Starting forward modeling'
    if (genpar%twoD) then
       call propagator_acoustic(                          &
       & FD_acoustic_init_coefs,                          &
       & FD_2nd_2D_derivatives_scalar_forward_grid_noomp, &
       & Injection_sinc_noomp,                            &
       & FD_2nd_time_derivative_grid_noomp,               &
       & FDswaptime_pointer,                              &
       & bounds,model,elev,genpar,                        &
       & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,     &
       & ExtractData=Extraction_array_sinc_noomp,ExtractWave=Extraction_wavefield)
    else
       call propagator_acoustic(                          &
       & FD_acoustic_init_coefs,                          &
       & FD_2nd_3D_derivatives_scalar_forward_grid_noomp, &
       & Injection_sinc_noomp,                            &
       & FD_2nd_time_derivative_grid_noomp,               &
       & FDswaptime_pointer,                              &
       & bounds,model,elev,genpar,                        &
       & sou=sourcevec,wfld=wfld_fwd,datavec=datavec,     &
       & ExtractData=Extraction_array_sinc_noomp,ExtractWave=Extraction_wavefield)
    end if
    if (genpar%verbose) write(0,*) 'INFO: Done with forward modeling'

!    if (ishot.eq.5) then
!       do i=1,genpar%ntsnap
!          call srite('wfld',wfld_fwd%wave(1:model%nz,1:model%nxw,1:model%nyw,i,1),4*model%nz*model%nxw*model%nyw)
!       end do
!       do j=1,size(datavec)
!           call srite('shotout',datavec(j)%trace(:,1),4*sourcevec(1)%dimt%nt)
!        end do
!        
!       
!       call to_history('n1',sourcevec(1)%dimt%nt,'shotout')
!       call to_history('n2',size(datavec),'shotout')
!       call to_history('n1',model%nz,'wfld')
!       call to_history('n2',model%nxw,'wfld')
!       call to_history('n3',model%nyw,'wfld')
!       call to_history('n4',genpar%ntsnap,'wfld')
!    end if
!
!    deallocate(wfld_fwd%wave)

  end subroutine AMOD_to_memory

end module to_disk_to_memory_AMOD_mod
