program Acoustic_modeling

  use sep
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
  integer :: ntsnap

  call sep_init()
  
  call from_history('n1',dat%nt)
  call from_history('d1',dat%dt)
  allocate(dat%source(dat%nt))
  call sreed('in',dat%source,4*dat%nt)

  genpar%twoD=.false.

!  genpar%twoD=.true.

  if (.not.genpar%twoD) then     
     genpar%nbound=4
  else
     genpar%nbound=0
  end if

  genpar%dt=dat%dt
  genpar%dx=10
  genpar%dy=10
  genpar%dz=10

  genpar%nt=dat%nt
  genpar%lsinc=7
  genpar%ntaper=10
  genpar%snapi=4

  genpar%rec_type=0
  genpar%surf_type=0
  genpar%shot_type=0

  mod%nx=201
  mod%ny=201
  mod%nz=101

  genpar%ntsnap=int(genpar%nt/genpar%snapi)

  if (genpar%twoD) then
     mod%ny=1
     genpar%dy=1.
  end if

  mod%dx=genpar%dx
  mod%dy=genpar%dy
  mod%dz=genpar%dz

  dat%nx=mod%nx
  dat%ny=mod%ny

  elev%dz=genpar%dz

  allocate(dat%data(dat%nt,dat%nx,dat%ny))
  allocate(dat%wave(mod%nz,mod%nx,mod%ny,genpar%ntsnap))
  
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
  
  do i=1,dat%ny
     call srite('data',dat%data(1:dat%nt,1:dat%nx,i),4*dat%nx*dat%nt)
  end do

  do i=1,genpar%ntsnap
     call srite('wave',dat%wave(1:mod%nz,1:mod%nx,1:mod%ny,i),4*mod%nx*mod%ny*mod%nz)
  end do

  call to_history('n1',dat%nt,'data')
  call to_history('n2',dat%nx,'data')
  call to_history('n3',dat%ny,'data')

  call to_history('n1',mod%nz,'wave')
  call to_history('n2',mod%nx,'wave')
  call to_history('n3',mod%ny,'wave')
  call to_history('n4',dat%nt,'wave')
  call deallocateModelSpace(mod)
  call deallocateDataSpace(dat)
end program Acoustic_modeling
