! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
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
  type(ModelSpace), target   :: mod
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


  mod%rhotag='rho'
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
 
  if ((genpar%task.eq.'INV').or.(genpar%task.eq.'inv')) then
     call Init_BandPassParam(bpparam)
     call BandpassSouTraces(bpparam,shotgath,sourcegath)
     call Init_SmoothingParam(smoothpar,mod,bpparam)
     call read_inv_params(invparam)
     if (invparam%wantreg.or.invparam%wantlog) call Init_SparseRegularization(sparseparam,mod)
     invparam%ntotaltraces=ntotaltraces
     invparam%n1=genpar%nt
     call Init_Inversion_Array(mod,invparam)

     if (invparam%wantreg.or.invparam%wantlog) call Init_SparseRegularization_Array(sparseparam,mod)
     call Init_MuteParam(mutepar)
     call MuteParam_compute_mask(mutepar,shotgath,sourcegath)
     call compute_fct_gdt_dmute_smooth_inv_init(mutepar,smoothpar,invparam)
!     if (invparam%wantreg.or.invparam%wantlog) then
!        call compute_fct_gdt_sparsepar_init(sparseparam)
!        call l_bfgs(invparam,compute_fct_gdt_reg,x)
!     else
        call l_bfgs(invparam,compute_fct_gdt,mod)
!     end if

     ! Freeze velocity model again, and make sure we don't update where we freeze
     call Freeze(invparam%parmin(1),invparam%parmax(1),invparam%modmask(:,1),&
     &           invparam%freeze_soft,mod%vel,mod%nx*mod%ny*mod%nz,invparam%modinit(:,1))
     call ApplyMask(mod%nx*mod%ny*mod%nz,invparam%modmask(:,1),mod%vel,invparam%modinit(:,1))
     if (invparam%nparam.eq.2) then
        if (invparam%vprho_param.eq.0) then
           call Freeze(invparam%parmin(2),invparam%parmax(2),invparam%modmask(:,2),&
           &           invparam%freeze_soft,mod%rho,mod%nx*mod%ny*mod%nz,invparam%modinit(:,2))
           call ApplyMask(mod%nx*mod%ny*mod%nz,invparam%modmask(:,2),mod%rho,invparam%modinit(:,2))
        else
           call Freeze_imp(invparam%parmin(2),invparam%parmax(2),invparam%modmask(:,2),&
           &           invparam%freeze_soft,mod%imp,mod%vel,mod%nx*mod%ny*mod%nz,invparam%modinit(:,2))
           call ApplyMask(mod%nx*mod%ny*mod%nz,invparam%modmask(:,2),mod%imp,invparam%modinit(:,2))
        end if
     end if

     call MuteParam_deallocate(mutepar)
     if (invparam%wantreg.or.invparam%wantlog) call SparseRegularization_filter_close(sparseparam)
     call srite('inv',mod%vel,4*mod%nz*mod%nx*mod%ny)
     if (invparam%nparam.eq.2) then
        if (invparam%vprho_param.eq.0) then
           call srite('inv',mod%rho,4*mod%nz*mod%nx*mod%ny)
        else
           call srite('inv',mod%imp,4*mod%nz*mod%nx*mod%ny)
        end if
     end if
     call create_header(mod,invparam,genpar,ntotaltraces)
     deallocate(invparam%modmask,invparam%modinit)
     deallocate(invparam%parmin,invparam%parmax)
  else if ((genpar%task.eq.'MOD').or.(genpar%task.eq.'mod')) then
     stat=compute_mod()  
  else if ((genpar%task.eq.'MIG').or.(genpar%task.eq.'mig')) then
     stat=compute_mig()
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
  if (allocated(mod%rho)) deallocate(mod%rho)
  if (allocated(mod%imp)) deallocate(mod%imp)

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
        call to_history('n4',invparam%nparam,'gradient')
        call to_history('d1',mod%dz,'gradient')
        call to_history('d2',mod%dx,'gradient')
        call to_history('d3',mod%dy,'gradient')
        call to_history('o1',mod%oz,'gradient')
        call to_history('o2',mod%ox,'gradient')
        call to_history('o3',mod%oy,'gradient')
        call to_history('n5',invparam%iter,'gradient')
        if (invparam%nparam.eq.2) then
           if (invparam%vprho_param.eq.0) then
              call to_history('label4','vel:rho','gradient')
           else
              call to_history('label4','vel:imp','gradient')
           endif
        end if
     end if
     
     if (exist_file('model')) then
        call to_history('n1',mod%nz,'model')
        call to_history('n2',mod%nx,'model')
        call to_history('n3',mod%ny,'model')
        call to_history('n4',invparam%nparam,'model')
        call to_history('d1',mod%dz,'model')
        call to_history('d2',mod%dx,'model')
        call to_history('d3',mod%dy,'model')
        call to_history('o1',mod%oz,'model')
        call to_history('o2',mod%ox,'model')
        call to_history('o3',mod%oy,'model')
        call to_history('n5',invparam%iter,'model')
        if (invparam%nparam.eq.2) then
           if (invparam%vprho_param.eq.0) then
              call to_history('label4','vel:rho','model')
           else
              call to_history('label4','vel:imp','model')
           endif
        end if
     end if

     if (exist_file('function')) then
        call to_history('n1',invparam%iter+1,'function')
     end if
     
     call to_history('n1',mod%nz,'inv')
     call to_history('n2',mod%nx,'inv')
     call to_history('n3',mod%ny,'inv')
     call to_history('n4',invparam%nparam,'inv')
     call to_history('d1',mod%dz,'inv')
     call to_history('d2',mod%dx,'inv')
     call to_history('d3',mod%dy,'inv')
     call to_history('o1',mod%oz,'inv')
     call to_history('o2',mod%ox,'inv')
     call to_history('o3',mod%oy,'inv')
     if (invparam%nparam.eq.2) then
        if (invparam%vprho_param.eq.0) then
           call to_history('label4','vel:rho','inv')
        else
           call to_history('label4','vel:imp','inv')
        endif
     end if

     if (exist_file('residual')) then
        call to_history('n1',genpar%nt,'residual')
        call to_history('n2',ntotaltraces,'residual')
        call to_history('n3',invparam%iter,'residual')
        call to_history('d1',genpar%dt,'residual')
        call to_history('d2',1.,'residual')
        call to_history('o1',genpar%t0,'residual')
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

subroutine ApplyMask(m,mask,x,xinit)
  integer ::         m
  real,             dimension(m) :: mask,x
  double precision, dimension(m) :: xinit

  x=mask*x+(1-mask)*sngl(xinit)

end subroutine ApplyMask

subroutine Freeze(xmin,xmax,mask,type,x,m,xinit)
  real ::         xmin,xmax
  integer ::                            m,j,i
  logical ::                     type
  real, dimension(m) ::     mask,     x
  double precision, dimension(m) ::       xinit

  IF (type) THEN
     do i=1,M
        x(i)=max(xmin,x(i))
        x(i)=min(xmax,x(i))
     end do
  ELSE     
     do i=1,M
        if ((mask(i).ne.0.).and.(mask(i).ne.2.)) then
           x(i)=max(xmin,x(i))
           x(i)=min(xmax,x(i))
        end if
     end do
     do i=1,M
        if(mask(i).eq.0.) x(i)=sngl(xinit(i))
        if(mask(i).eq.2.) x(i)=sngl(xinit(i))
     end do
  ENDIF

end subroutine Freeze


subroutine Freeze_imp(xmin,xmax,mask,type,x,v,m,xinit)
  real ::             xmin,xmax
  integer ::                                m,j,i
  logical ::                         type
  real, dimension(m) ::         mask,     x,v
  double precision, dimension(m) ::           xinit

  !xmin is rhomin
  !xmax is rhomax
  IF (type) THEN
     do i=1,M
        x(i)=max(xmin*v(i),x(i))
        x(i)=min(xmax*v(i),x(i))
     end do
  ELSE     
     do i=1,M
        if ((mask(i).ne.0.).and.(mask(i).ne.2.)) then
           x(i)=max(xmin*v(i),x(i))
           x(i)=min(xmax*v(i),x(i))
        end if
     end do
     do i=1,M
        if(mask(i).eq.0.) x(i)=sngl(xinit(i))
        if(mask(i).eq.2.) x(i)=sngl(xinit(i))
     end do
  ENDIF

end subroutine Freeze_imp
