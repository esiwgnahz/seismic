program TWODAFWI

  use sep
  use omp_lib
  use Inversion_mod
  use ReadParam_mod
  use DataSpace_types
  use Readsouvelrho_mod
  use ExtractPadModel_mod

  use Mute_gather
  use Inversion_types
  use ModelSpace_types
  use GeneralParam_types

  use compute_function_gradient

  implicit none

  type(GeneralParam) :: genpar
  type(ModelSpace)   :: mod
  type(FDbounds)     :: bounds
  type(ModelSpace_elevation) :: elev
  type(MuteParam)    :: mutepar

  type(ModelSpace),  dimension(:), allocatable :: modgath
  type(GatherSpace), dimension(:), allocatable :: shotgath
  type(TraceSpace),  dimension(:), allocatable :: source
  type(TraceSpace),  dimension(:), allocatable :: sourcegath
  type(FDbounds),    dimension(:), allocatable :: boundsgath
  type(GeneralParam),dimension(:), allocatable :: genpargath

  type(InversionParam)                         :: invparam

  real, dimension(:),              allocatable :: grad
  integer :: i,ntotaltraces,stat
  double precision    :: f
  call sep_init()

  mod%veltag='vel'
  mod%waFtag='wave_fwd'

  if (.not.exist_file(mod%veltag)) call erexit('ERROR: need velocity file')


  call read_3D_params(genpar)
  call readsou(source,genpar)
  call readgathercoords(shotgath,sourcegath,genpar)
  if (genpar%task.ne.'MOD') call readgathertraces(shotgath,source)
  !call readgathertraces(shotgath,source)
  call copysou2sougath(source,sourcegath)
  
  allocate(modgath(size(shotgath)))   ! Each shot has one model
  allocate(boundsgath(size(shotgath)))! Each shot has different bounds
  allocate(genpargath(size(shotgath)))! Each shot has different parameters
  do i=1,size(shotgath)
     call extract_coord_source_receiver_patch(shotgath(i)%gathtrace,sourcegath(i),modgath(i),genpar)
  end do

  call read_vel(mod,genpar)

  call compute_fct_gdt_init(mod,modgath,genpar,genpargath,shotgath,sourcegath,bounds,boundsgath)
  if (genpar%task.eq.'INV') then

     call read_inv_params(invparam)
     call Init_Inversion_Array(mod,invparam)
     call Init_MuteParam(mutepar)
     call MuteParam_compute_mask(mutepar,shotgath,sourcegath)
     call compute_fct_gft_dmute_init(mutepar)
     call l_bfgs(invparam,compute_fct_gdt,mod%vel)
     call MuteParam_deallocate(mutepar)
     call to_history('n1',mod%nz,'grad')
     call to_history('n2',mod%nx,'grad')
     call to_history('n3',mod%ny,'grad')
     
     call srite('inv',mod%vel,4*mod%nz*mod%nx*mod%ny)
     call to_history('n1',mod%nz,'inv')
     call to_history('n2',mod%nx,'inv')
     call to_history('n3',mod%ny,'inv')
     deallocate(invparam%vpmask,invparam%vpinit)
  else if (genpar%task.eq.'MOD') then
     stat=compute_mod()
  end if
  call compute_fct_gdt_nullify()

  do i=1,size(shotgath)
     call deallocateGatherSpace(shotgath(i))
     call deallocateTraceSpace(sourcegath(i))
     call deallocateModelSpace(modgath(i))
  end do

  do i=1,size(sourcegath)
     call deallocateTraceSpace(sourcegath(i))
  end do
  call deallocateTraceSpace(source(1))
  deallocate(mod%vel)

end program TWODAFWI
