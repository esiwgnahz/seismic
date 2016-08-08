program Acoustic_modeling

  use FDcoefs_assign
  use Propagator_mod

  use DataSpace_types
  use ModelSpace_types
  use GeneralParam_types

  implicit none

  type(GeneralParam) :: genpar
  type(ModelSpace)   :: mod
  type(DataSpace)    :: dat
  type(FDbounds)     :: bounds
  type(ModelSpace_elevation) :: elev

  integer :: i,j,k 

  genpar%twoD=.false.

!  genpar%twoD=.true.

  if (.not.genpar%twoD) then     
     genpar%nbound=4
  else
     genpar%nbound=0
  end if

  genpar%dt=0.004
  genpar%dx=10
  genpar%dy=10
  genpar%dz=10

  genpar%nt=100
  genpar%lsinc=7
  genpar%ntaper=10

  genpar%rec_type=0
  genpar%surf_type=0
  genpar%shot_type=0

  mod%nx=51
  mod%ny=51
  mod%nz=51

  if (genpar%twoD) then
     mod%ny=1
     genpar%dy=1.
  end if

  mod%dx=genpar%dx
  mod%dy=genpar%dy
  mod%dz=genpar%dz

  dat%nt=genpar%nt
  dat%dt=genpar%dt
  dat%nx=mod%nx
  dat%ny=mod%ny

  elev%dz=genpar%dz

  allocate(dat%source(dat%nt))
  dat%source=0.
  dat%source(int(dat%nt/4))=1.

  allocate(dat%data(dat%nt,dat%nx,dat%ny))
  allocate(dat%wave(mod%nz,mod%nx,mod%ny,dat%nt))
  
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
  allocate(elev%elev(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  allocate(elev%elev_rec(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))
  allocate(elev%elev_sou(bounds%nmin2:bounds%nmax2, bounds%nmin3:bounds%nmax3))

!  do k=bounds%nmin3,bounds%nmax3
!     do j=bounds%nmin2,bounds%nmax2
!        do i=bounds%nmin1,bounds%nmax1
!           mod%vel(i,j,k)=j
!        end do
!     end do
!  end do

  mod%vel=2500. 
  elev%elev=0.
  elev%elev_rec=0.
  elev%elev_sou=0.
  elev%rec_z=0.
  elev%shot_z=0.
  elev%o1model=0. 
  elev%ishot_x=int(mod%nx/2)
  elev%ishot_y=int(mod%ny/2)

  write(0,*) 'before scalar wave propagator'
  if (genpar%twoD) then
     call propagator_acoustic(FD_acoustic_init_coefs, &
     & FD_2nd_2D_derivatives_scalar_forward,          &
     & Injection_source,                              &
     & FD_2nd_time_derivative,                        &
     & FDswaptime,bounds,mod,dat,elev,genpar)
  else
     call propagator_acoustic(FD_acoustic_init_coefs, &
     & FD_2nd_3D_derivatives_scalar_forward,          &
     & Injection_source,                              &
     & FD_2nd_time_derivative_omp,                    &
     & FDswaptime_omp,bounds,mod,dat,elev,genpar)
  end if
  write(0,*) 'after scalar wave propagator'
  
  call deallocateModelSpace(mod)
  call deallocateDataSpace(dat)
end program Acoustic_modeling
