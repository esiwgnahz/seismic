program Acoustic_modeling

  use sep
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

  type(WaveSpace)                             :: wfld
  type(TraceSpace), dimension(:), allocatable :: datavec
  type(TraceSpace)                            :: sourcevec

  integer :: i,j,k 
  integer :: ntsnap

  logical :: withRho
  call sep_init()
  
  call from_history('n1',sourcevec%dimt%nt)
  call from_history('d1',sourcevec%dimt%dt)
  allocate(sourcevec%trace(sourcevec%dimt%nt,1))
  call sreed('in',sourcevec%trace(:,1),4*sourcevec%dimt%nt)

  call from_param('withRho',withRho,.false.)
  if (withRho) then
     genpar%coefpower=1
  else
     genpar%coefpower=2
  end if
  call from_param('twoD',genpar%twoD,.false.)

  allocate(datavec(201))
  do i=1,size(datavec)
     allocate(datavec(i)%trace(sourcevec%dimt%nt,1)) ! 1 component trace
  end do

  if (.not.genpar%twoD) then     
     genpar%nbound=4
  else
     genpar%nbound=0
  end if

  genpar%dt=sourcevec%dimt%dt
  genpar%dx=10
  genpar%dy=10
  genpar%dz=10

  genpar%nt=sourcevec%dimt%nt
  genpar%lsinc=7
  genpar%ntaper=10
  genpar%snapi=4

  genpar%rec_type=0
  genpar%surf_type=0
  genpar%shot_type=0

  mod%nx=201
  mod%ny=201
  mod%nz=101

  sourcevec%coord(1)=130.4   ! Z
  sourcevec%coord(2)=546.7  ! X
  sourcevec%coord(3)=800.1  ! Y

  genpar%ntsnap=int(genpar%nt/genpar%snapi)

  if (genpar%twoD) then
     mod%ny=1
     genpar%dy=1.
  end if

  mod%dx=genpar%dx
  mod%dy=genpar%dy
  mod%dz=genpar%dz

  elev%omodel=0.
  elev%delta(1)=genpar%dz
  elev%delta(2)=genpar%dx
  elev%delta(3)=genpar%dy

  do i=1,size(datavec)     
     datavec(i)%coord(1)=0.  ! Z
     datavec(i)%coord(2)=elev%omodel(2)+(i-1)*mod%dx ! X
     datavec(i)%coord(3)=800.1 ! Y
  end do

  allocate(wfld%wave(mod%nz,mod%nx,mod%ny,genpar%ntsnap,1))
  
  if (genpar%shot_type.gt.0 .or. genpar%surf_type.gt.0) then
     if (genpar%shot_type.gt.0) bounds%nmin1 = -(genpar%ntaper+genpar%lsinc/2+2)+1
     if (genpar%surf_type.gt.0) bounds%nmin1 = -genpar%lsinc+1
  else
     bounds%nmin1 = -genpar%ntaper+1
  endif
  bounds%nmax1 =  mod%nz+genpar%ntaper
  bounds%nmin2 = -genpar%ntaper+1
  bounds%nmax2 =  mod%nx+genpar%ntaper
  if (.not. genpar%twoD) then
     bounds%nmin3 = -genpar%ntaper+1
     bounds%nmax3 =  mod%ny+genpar%ntaper
  else
     bounds%nmin3 = 1
     bounds%nmax3 = 1
  end if

  write(0,*) 'bounds%nmin1',bounds%nmin1,'bounds%nmax1',bounds%nmax1
  write(0,*) 'bounds%nmin2',bounds%nmin2,'bounds%nmax2',bounds%nmax2
  write(0,*) 'bounds%nmin3',bounds%nmin3,'bounds%nmax3',bounds%nmax3

  allocate(mod%vel(bounds%nmin1:bounds%nmax1, bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  allocate(mod%rho(bounds%nmin1:bounds%nmax1, bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  allocate(mod%rho2(bounds%nmin1:bounds%nmax1, bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  allocate(elev%elev(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))

  mod%vel(bounds%nmin1:40,:,:)=2400.
  mod%vel(40:80,:,:)=2400. 
  mod%vel(80:,:,:)=2400. 
  mod%rho(bounds%nmin1:40,:,:)=1.
  mod%rho(40:80,:,:)=2.
  mod%rho(80:,:,:)=4.
  call Interpolate(mod,bounds)

  write(0,*) 'before wave propagator'
  if (genpar%twoD) then
     if (.not.withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_2D_derivatives_scalar_forward,          &
        & Extraction_array_simple,                       &
        & Injection_source_sinc_xz,                      &
        & FD_2nd_time_derivative,                        &
        & FDswaptime,                                    &
        & bounds,mod,sourcevec,datavec,wfld,elev,genpar)
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_2D_derivatives_acoustic_forward,            &
        & Extraction_array_simple,                       &
        & Injection_source_rho_sinc_xz,                  &
        & FD_2nd_time_derivative,                        &
        & FDswaptime,                                    &
        & bounds,mod,sourcevec,datavec,wfld,elev,genpar)
     end if
  else
     if (.not.withRho) then
        call propagator_acoustic(                        &
        & FD_acoustic_init_coefs,                        &
        & FD_2nd_3D_derivatives_scalar_forward,          &
        & Extraction_array_simple,                       &
        & Injection_source_sinc_xyz,                     &
        & FD_2nd_time_derivative_omp,                    &
        & FDswaptime_omp,                                &
        & bounds,mod,sourcevec,datavec,wfld,elev,genpar)
     else
        call propagator_acoustic(                        &
        & FD_acoustic_rho_init_coefs,                    &
        & FD_3D_derivatives_acoustic_forward,            &
        & Extraction_array_simple,                       &
        & Injection_source_rho_sinc_xz,                  &
        & FD_2nd_time_derivative,                        &
        & FDswaptime,                                    &
        & bounds,mod,sourcevec,datavec,wfld,elev,genpar)
     end if
  end if
  write(0,*) 'afterwave propagator'
  
  do i=1,size(datavec)
     call srite('data',datavec(i)%trace(:,1),4*sourcevec%dimt%nt)
  end do

  do i=1,genpar%ntsnap
     call srite('wave',wfld%wave(1:mod%nz,1:mod%nx,1:mod%ny,i,1),4*mod%nx*mod%ny*mod%nz)
  end do

  call to_history('n1',sourcevec%dimt%nt,'data')
  call to_history('n2',size(datavec),'data')
  call to_history('d1',sourcevec%dimt%dt,'data')
  call to_history('d2',1.,'data')
  call to_history('o1',0.,'data')
  call to_history('o2',0.,'data')

  call to_history('n1',mod%nz,'wave')
  call to_history('n2',mod%nx,'wave')
  call to_history('n3',mod%ny,'wave')
  call to_history('d1',mod%dz,'wave')
  call to_history('d2',mod%dx,'wave')
  call to_history('d3',mod%dy,'wave')
  call to_history('o1',elev%omodel(1),'wave')
  call to_history('o2',elev%omodel(2),'wave')
  call to_history('o3',elev%omodel(3),'wave')
  call to_history('n4',genpar%ntsnap,'wave')

  do i=1,size(datavec)
     call deallocateTraceSpace(datavec(i))
  end do
  call deallocateWaveSpace(wfld)
  call deallocateTraceSpace(sourcevec)
  deallocate(datavec)

end program Acoustic_modeling
