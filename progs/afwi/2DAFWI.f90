program TWODAFWI

  use sep
  use omp_lib
  use Inversion_mod
  use ReadParam_mod
  use DataSpace_types
  use Readsouvelrho_mod
  use ExtractPadModel_mod

  use Mute_gather
  use Bandpass_mod
  use Smoothing_mod
  use Inversion_types
  use ModelSpace_types
  use GeneralParam_types
  use Sparse_regularization_mod
  use compute_function_gradient

  implicit none

  type(GeneralParam) :: genpar
  type(ModelSpace)   :: mod
  type(FDbounds)     :: bounds
  type(ModelSpace_elevation) :: elev

  type(ModelSpace),  dimension(:), allocatable :: modgath
  type(GatherSpace), dimension(:), allocatable :: shotgath
  type(TraceSpace),  dimension(:), allocatable :: source
  type(TraceSpace),  dimension(:), allocatable :: sourcegath
  type(FDbounds),    dimension(:), allocatable :: boundsgath
  type(GeneralParam),dimension(:), allocatable :: genpargath

  type(SparseRegParam)                         :: sparseparam
  type(SmoothingParam)                         :: smoothpar
  type(InversionParam)                         :: invparam
  type(BandPassParam)                          :: bpparam
  type(MuteParam)                              :: mutepar

  real, dimension(:),              allocatable :: grad
  integer :: i,j,ntotaltraces,stat
  double precision    :: f
  call sep_init()


  mod%veltag='vel'
  mod%waFtag='wave_fwd'

  if (.not.exist_file(mod%veltag)) call erexit('ERROR: need velocity file')

  call from_aux('coordfile','n2',ntotaltraces)
  call read_3D_params(genpar)
  call readsou(source,genpar)
  call readgathercoords(shotgath,sourcegath,genpar)
  if (genpar%task.ne.'MOD') call readgathertraces(shotgath,source)
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
     call Init_BandPassParam(bpparam)
     call BandpassSouTraces(bpparam,shotgath,sourcegath)
     call Init_SmoothingParam(smoothpar,mod,bpparam)
     call read_inv_params(invparam)
     if (invparam%eps.ne.0.) call Init_SparseRegularization(sparseparam,mod)
     invparam%ntotaltraces=ntotaltraces
     invparam%n1=genpar%nt
     call Init_Inversion_Array(mod,invparam)
     call Init_MuteParam(mutepar)
     call MuteParam_compute_mask(mutepar,shotgath,sourcegath)
     call compute_fct_gdt_dmute_smooth_inv_init(mutepar,smoothpar,invparam)
     if (invparam%eps.ne.0.) then
        call compute_fct_gdt_sparsepar_init(sparseparam)
        call l_bfgs(invparam,compute_fct_gdt_reg,mod%vel)
     else
        call l_bfgs(invparam,compute_fct_gdt,mod%vel)
     end if
     call MuteParam_deallocate(mutepar)
     if (invparam%eps.ne.0.) call SparseRegularization_filter_close(sparseparam)
     call srite('inv',mod%vel,4*mod%nz*mod%nx*mod%ny)
     call create_header(mod,invparam,genpar,ntotaltraces)
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

subroutine create_header(mod,invparam,genpar,ntotaltraces)
  use sep
  use Inversion_types
  use ModelSpace_types
  use GeneralParam_types
  
  type(ModelSpace)   :: mod
  type(GeneralParam) :: genpar
  type(InversionParam)                         :: invparam
  integer:: ntotaltraces
  
     if (exist_file('gradient')) then
        call to_history('n1',mod%nz,'gradient')
        call to_history('n2',mod%nx,'gradient')
        call to_history('n3',mod%ny,'gradient')
        call to_history('d1',mod%dz,'gradient')
        call to_history('d2',mod%dx,'gradient')
        call to_history('d3',mod%dy,'gradient')
        call to_history('o1',mod%oz,'gradient')
        call to_history('o2',mod%ox,'gradient')
        call to_history('o3',mod%oy,'gradient')
        call to_history('n4',invparam%iter,'gradient')
     end if
     
     if (exist_file('model')) then
        call to_history('n1',mod%nz,'model')
        call to_history('n2',mod%nx,'model')
        call to_history('n3',mod%ny,'model')
        call to_history('d1',mod%dz,'model')
        call to_history('d2',mod%dx,'model')
        call to_history('d3',mod%dy,'model')
        call to_history('o1',mod%oz,'model')
        call to_history('o2',mod%ox,'model')
        call to_history('o3',mod%oy,'model')
        call to_history('n4',invparam%iter,'model')
     end if

     if (exist_file('function')) then
        call to_history('n1',invparam%iter+1,'function')
     end if
     
     call to_history('n1',mod%nz,'inv')
     call to_history('n2',mod%nx,'inv')
     call to_history('n3',mod%ny,'inv')
     call to_history('d1',mod%dz,'inv')
     call to_history('d2',mod%dx,'inv')
     call to_history('d3',mod%dy,'inv')
     call to_history('o1',mod%oz,'inv')
     call to_history('o2',mod%ox,'inv')
     call to_history('o3',mod%oy,'inv')

     if (exist_file('residual')) then
        call to_history('n1',genpar%nt,'residual')
        call to_history('n2',ntotaltraces,'residual')
        call to_history('n3',invparam%iter,'residual')
        call to_history('d1',genpar%dt,'residual')
        call to_history('d2',1.,'residual')
        call to_history('o1',0.,'residual')
        call to_history('o2',0.,'residual')
     end if

     if (exist_file('mute')) then
        call to_history('n1',genpar%nt,'mute')
        call to_history('n2',ntotaltraces,'mute')
        call to_history('d1',genpar%dt,'mute')
        call to_history('d2',1.,'mute')
        call to_history('o1',0.,'mute')
        call to_history('o2',0.,'mute')
     end if

end subroutine create_header
