module compute_function_gradient

  use sep
  use omp_lib
  use ReadParam_mod
  use DataSpace_types
  use Readsouvelrho_mod
  use ExtractPadModel_mod

  use Mute_gather
  use Smoothing_mod
  use Inversion_types
  use ModelSpace_types
  use GeneralParam_types

  use Sparse_regularization_mod

  use OF_Res_AdjSrc_mod
  use to_disk_to_memory_AMOD_mod
  use to_disk_to_memory_AGRAD_mod

  implicit none

  type(SparseRegParam),pointer,private :: sparseparam
  type(InversionParam),pointer,private :: invparam
  type(SmoothingParam),pointer,private :: smoothpar
  type(MuteParam),    pointer, private :: mutepar
  type(GeneralParam), pointer, private :: genpar
  type(ModelSpace),   pointer, private :: mod
  type(FDbounds),     pointer, private :: bounds

  type(ModelSpace),  dimension(:), pointer, private :: modgath
  type(GatherSpace), dimension(:), pointer, private :: shotgath
  type(TraceSpace),  dimension(:), pointer, private :: sourcegath
  type(FDbounds),    dimension(:), pointer, private :: boundsgath
  type(GeneralParam),dimension(:), pointer, private :: genpargath

contains

  subroutine compute_fct_gdt_sparsepar_init(sparsepar_in)
    type(SparseRegParam), target ::         sparsepar_in
    sparseparam=>sparsepar_in
  end subroutine compute_fct_gdt_sparsepar_init

  subroutine compute_fct_gdt_dmute_smooth_inv_init(mutepar_in,smooth_in,inv_init)
    type(MuteParam),      target ::                mutepar_in
    type(SmoothingParam), target ::                           smooth_in
    type(InversionParam), target ::                                     inv_init

    mutepar=>mutepar_in
    smoothpar=>smooth_in
    invparam=>inv_init
  end subroutine compute_fct_gdt_dmute_smooth_inv_init

  subroutine compute_fct_gdt_init(mod_in,modgath_in,genpar_in,genpargath_in,shotgath_in,sourcegath_in,bounds_in,boundsgath_in)

    type(GeneralParam), target :: genpar_in
    type(ModelSpace),   target :: mod_in
    type(FDbounds),     target :: bounds_in

    type(ModelSpace),  dimension(:), target :: modgath_in
    type(GatherSpace), dimension(:), target :: shotgath_in
    type(TraceSpace),  dimension(:), target :: sourcegath_in
    type(FDbounds),    dimension(:), target :: boundsgath_in
    type(GeneralParam),dimension(:), target :: genpargath_in

    genpar    =>genpar_in
    mod       =>mod_in
    bounds    =>bounds_in
    modgath   =>modgath_in
    shotgath  =>shotgath_in
    sourcegath=>sourcegath_in
    boundsgath=>boundsgath_in
    genpargath=>genpargath_in

  end subroutine compute_fct_gdt_init

  subroutine compute_fct_gdt_nullify()
    if(associated(sparseparam))nullify(sparseparam)
    if(associated(smoothpar))  nullify(smoothpar)
    if(associated(mutepar))    nullify(mutepar)
    if(associated(invparam))   nullify(invparam)
    if(associated(genpar))     nullify(genpar)
    if(associated(modgath))    nullify(modgath)
    if(associated(mod))        nullify(mod)
    if(associated(bounds))     nullify(bounds)
    if(associated(shotgath))   nullify(shotgath)
    if(associated(sourcegath)) nullify(sourcegath)
    if(associated(boundsgath)) nullify(boundsgath)
    if(associated(genpargath)) nullify(genpargath)
  end subroutine compute_fct_gdt_nullify

  function compute_fct_gdt_reg(grad,f,resigath) result(stat)
    real, dimension(:) ::      grad
    double precision   ::           f
    type(TraceSpace), dimension(:) :: resigath
    integer            ::                              stat
    integer            :: n1,ntotaltraces

    double precision                :: ftmp,ftmp1,scaling
    real, dimension(:), allocatable :: gtmp

    real, dimension(:), pointer :: gtmp1
    real, dimension(:), pointer :: gtmp2

    allocate(gtmp(mod%nx*mod%ny*mod%nz))

    gtmp=0.
    ftmp=0.
    ftmp1=0.
    f=0.

!    gtmp1=>gtmp(1:mod%nx*mod%ny*mod%nz)
!    if (invparam%nparam.eq.2) gtmp2=>gtmp(1+mod%nx*mod%ny*mod%nz:)

    stat=compute_fct_gdt(grad,f,resigath)

    scaling=mod%nx*mod%ny*mod%nz
    if (sparseparam%compute_eps) then
       if (invparam%eval.eq.0) then
          sparseparam%eps_log=0.
          sparseparam%eps=0.
       end if
    end if

    if (invparam%wantreg) then

       call SparseRegularization_apply(sparseparam,mod%vel,gtmp,ftmp)
       
       ftmp=ftmp/scaling
       gtmp=gtmp/sngl(scaling)
       
       if (sparseparam%compute_eps) then
          if (invparam%eval.eq.0) then
             if (sparseparam%ratio.eq.0.) then
                sparseparam%eps=0.
             else
                sparseparam%eps=abs(f/(sparseparam%ratio*ftmp))
             end if
             write(0,*) 'INFO:'
             write(0,*) 'INFO: ---- Regularization: setting epsilon ---'
             write(0,*) 'INFO: For ratio=',sparseparam%ratio,' eps=',sparseparam%eps
             write(0,*) 'INFO:'
          end if
       end if
       
       grad(1:mod%nx*mod%ny*mod%nz)=grad(1:mod%nx*mod%ny*mod%nz)-gtmp*sparseparam%eps

!       call srite('grad_reg',gtmp,4*size(gtmp))
!       call to_history('n1',mod%nz,'grad_reg')
!       call to_history('n2',mod%nx,'grad_reg')
!       call to_history('n3',mod%ny,'grad_reg')
    else
       sparseparam%eps=0.
    end if

    if (invparam%wantlog) then
       gtmp=0.

       call LogisticRegularization_apply(sparseparam,mod%vel,gtmp,ftmp1)
       
       ftmp1=ftmp1/scaling
       gtmp=gtmp/sngl(scaling)
       
       if (sparseparam%compute_eps) then
          if (invparam%eval.eq.0) then
             if (sparseparam%ratio_log.eq.0.) then
                sparseparam%eps_log=0.
             else
                sparseparam%eps_log=abs(f/(sparseparam%ratio_log*ftmp1))
             end if
             write(0,*) 'INFO:'
             write(0,*) 'INFO: ---- Logistic Reg.: setting epsilon ---'
             write(0,*) 'INFO: For ratio=',sparseparam%ratio_log,' eps=',sparseparam%eps_log
             write(0,*) 'INFO:'
          end if
       end if

       grad(1:mod%nx*mod%ny*mod%nz)=grad(1:mod%nx*mod%ny*mod%nz)-gtmp*sparseparam%eps_log
!       call srite('grad_reg_log',gtmp,4*size(gtmp))
!       call to_history('n1',mod%nz,'grad_reg_log')
!       call to_history('n2',mod%nx,'grad_reg_log')
!       call to_history('n3',mod%ny,'grad_reg_log')

    else
       sparseparam%eps_log=0.
    end if

    call mult_grad_mask(grad,invparam%modmask)
    f=f+ftmp*sparseparam%eps+ftmp1*sparseparam%eps_log
    
    deallocate(gtmp)

  end function compute_fct_gdt_reg

  function compute_fct_gdt(grad,f,resigath) result(stat)
    
    real, dimension(:) ::  grad
    double precision   ::       f
    type(TraceSpace), dimension(:) :: resigath
    integer            ::                          stat
    
    type(DataSpace)    :: dat

    type(TraceSpace),  dimension(:), allocatable :: dmodgath
    type(ModelSpace_elevation), dimension(:), allocatable :: elevgath
    type(WaveSpace), target                      :: wfld_fwd

    real, dimension(:,:,:,:,:), allocatable :: gradthread
    real, dimension(:,:,:,:),   allocatable :: illuthread
    double precision, dimension(:),       allocatable :: fthread

    real    :: d1
    integer :: i,j,k,l,m,n1,begi,endi
    integer :: ntsnap,ntotaltraces,index,indexp
    double precision :: memory_needed,gist,scaling

    real, dimension(:), allocatable          :: illu

    ntotaltraces=invparam%ntotaltraces
    allocate(elevgath(size(shotgath)))  ! Each shot has an elevation file

    n1=invparam%n1
    d1=genpar%dt
    genpar%ntsnap=int(genpar%nt/genpar%snapi)
    genpar%dt2=(genpar%dt*genpar%snapi)**2

    endi=0
    begi=0

    genpar%verbose=.false.
    genpar%optim=.false.
    memory_needed=0.

    genpar%nthreads=min(genpar%nthreads,size(shotgath))
    call omp_set_num_threads(genpar%nthreads)

    ! image and illumination for each thread
    allocate(gradthread(mod%nz,mod%nx,mod%ny,invparam%nparam,genpar%nthreads)); gradthread=0
    allocate(illuthread(mod%nz,mod%nx,mod%ny,genpar%nthreads));                 illuthread=0
    allocate(fthread(genpar%nthreads));fthread=0.
    
    f=0.
    grad=0.

    ! Convert impedances back to density for propagator
    if ((invparam%nparam.eq.2).and.(invparam%vprho_param.eq.1)) mod%rho=mod%imp/mod%vel

    !$OMP PARALLEL DO PRIVATE(i,k,j,dmodgath,begi,endi,wfld_fwd)
    do i=1,size(shotgath)
!
       call copy_window_vel_gath(mod,modgath(i),genpar,genpargath(i),boundsgath(i),i)
       genpargath(i)%ntsnap=int(genpargath(i)%nt/genpargath(i)%snapi)
!
       allocate(elevgath(i)%elev(boundsgath(i)%nmin2:boundsgath(i)%nmax2, boundsgath(i)%nmin3:boundsgath(i)%nmax3))
       elevgath(i)%elev=0.
!
       allocate(dmodgath(size(shotgath(i)%gathtrace)))
       do j=1,size(shotgath(i)%gathtrace)
          allocate(dmodgath(j)%trace(n1,1))
       end do
       dmodgath=shotgath(i)%gathtrace
       do j=1,size(shotgath(i)%gathtrace)
          dmodgath(j)%trace=0.
       end do

       ! Forward: modeling
        call AMOD_to_memory(modgath(i),genpargath(i),dat,boundsgath(i),elevgath(i),dmodgath,sourcegath,wfld_fwd,i) 

       begi=shotgath(i)%begi

       call Compute_OF_RES_ADJ(begi,invparam,shotgath(i)%gathtrace,dmodgath,mutepar%maskgath(i)%gathtrace,resigath,fthread(omp_get_thread_num ()+1))

       ! Backward: imaging
       call AGRAD_to_memory(modgath(i),genpargath(i),dat,boundsgath(i),elevgath(i),dmodgath,wfld_fwd,invparam%nparam)
!       !write(0,*) minval(modgath(i)%imagesmall_nparam),maxval(modgath(i)%imagesmall_nparam)
       ! Copy to final image space and convert gradient to match parameterization
       call mod_copy_image_nparam(modgath(i),gradthread(:,:,:,:,omp_get_thread_num()+1),illuthread(:,:,:,omp_get_thread_num()+1),invparam%vprho_param)

       ! Deallocate arrays
       deallocate(modgath(i)%imagesmall_nparam)
       deallocate(modgath(i)%illumsmall)
       call deallocateModelSpace_elev(elevgath(i))       
       do k=1,shotgath(i)%ntraces
          call deallocateTraceSpace(dmodgath(k))
       end do
       deallocate(modgath(i)%vel)
       deallocate(dmodgath)
       if (allocated(modgath(i)%rho)) deallocate(modgath(i)%rho)
       if (allocated(modgath(i)%imp)) deallocate(modgath(i)%imp)
       if (allocated(modgath(i)%rho2)) deallocate(modgath(i)%rho2)

    end do
    !$OMP END PARALLEL DO
    
    allocate(illu(mod%nz*mod%nx*mod%ny)); illu=0.
    ! Add all images for each thread to final image
    do i=1,genpar%nthreads

       f   =f   +fthread(i)
       !$OMP PARALLEL DO PRIVATE(m,k,j,l,index)
       do m=1,invparam%nparam
          do k=1,mod%ny
             do j=1,mod%nx
                do l=1,mod%nz
                   index=l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx+(m-1)*mod%nz*mod%nx*mod%ny
                   grad(index)=grad(index)+gradthread(l,j,k,m,i)            
                end do
             end do
          end do
       end do
       !$OMP END PARALLEL DO 

       !$OMP PARALLEL DO PRIVATE(k,j,l)
          do k=1,mod%ny
             do j=1,mod%nx
                do l=1,mod%nz
                   illu(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)=illu(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)+illuthread(l,j,k,i)
                end do
             end do
          end do
       !$OMP END PARALLEL DO 

    end do

    if ((invparam%nparam.eq.2).and.(invparam%vprho_param.eq.1)) then
       !$OMP PARALLEL DO PRIVATE(k,j,l,index,indexp)
       do k=1,mod%ny
          do j=1,mod%nx
             do l=1,mod%nz
                index=l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx
                indexp=index+mod%nz*mod%nx*mod%ny
                grad(index)=grad(index)-(mod%rho(l,j,k)/mod%vel(l,j,k))*grad(indexp) ! grad(v)=grad(v)-Ip/vp^2*grad(rho)
                grad(indexp)=(1/mod%vel(l,j,k))*grad(indexp)                         ! grad(Ip)=1/vp*grad(rho)
             end do
          end do
       end do
       !$OMP END PARALLEL DO 
    end if

    illu=illu**invparam%illupow

    ! Scaling
    scaling=dble(2*n1*ntotaltraces)
    f=f/scaling
    grad=2*grad/sngl(scaling)
    illu=(illu+maxval(illu)/10000)/sqrt(sum(dprod(illu,illu))/size(illu))

    grad(1:mod%nz*mod%nx*mod%ny)=invparam%sigma(1)*grad(1:mod%nz*mod%nx*mod%ny)/illu
    if (invparam%nparam.eq.2) then
       grad(1+mod%nz*mod%nx*mod%ny:)=invparam%sigma(2)*grad(1+mod%nz*mod%nx*mod%ny:)/illu
    end if

    call mult_grad_mask(grad,invparam%modmask)
    call triangle2(smoothpar,grad(1:mod%nz*mod%nx*mod%ny))
    if (invparam%nparam.eq.2) then
       call triangle2(smoothpar,grad(1+mod%nz*mod%nx*mod%ny:))
    end if
    call mult_grad_mask(grad,invparam%modmask)
    deallocate(gradthread,illu,illuthread,fthread)
    deallocate(elevgath)

    stat=0

  end function compute_fct_gdt

  subroutine mult_grad_mask(grad,mask)
    real, dimension(:) ::   grad
    real, dimension(:,:) ::      mask

    integer::l,k,j,i

    !$OMP PARALLEL DO PRIVATE(i,j,k,l)
    do l=1,invparam%nparam
       do k=1,mod%ny
          do j=1,mod%nx
             do i=1,mod%nz
                grad(i+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx+(l-1)*mod%nz*mod%nx*mod%ny)=grad(i+(j-1)*mod%nz+(k-1)* &
               &             mod%nz*mod%nx+(l-1)*mod%nz*mod%nx*mod%ny)*mask(i+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx,l)
             end do
          end do
       end do
    end do
  end subroutine mult_grad_mask

  function compute_mig() result(stat)
    
    integer            :: stat
    type(DataSpace)    :: dat

    type(TraceSpace),  dimension(:), allocatable :: dmodgath
    type(ModelSpace_elevation), dimension(:), allocatable :: elevgath
    type(WaveSpace), target                      :: wfld_fwd

    real, dimension(:,:,:,:,:), allocatable :: imagthread
    real, dimension(:,:,:,:),   allocatable :: illuthread

    real    :: d1
    integer :: i,j,k,l,n1
    integer :: ntsnap,ntotaltraces,nprocessed
    double precision :: memory_needed,gist,scaling

    real, dimension(:), allocatable          :: illu,imag

    allocate(elevgath(size(shotgath)))  ! Each shot has an elevation file

    d1=genpar%dt
    genpar%ntsnap=int(genpar%nt/genpar%snapi)
    genpar%dt2=(genpar%dt*genpar%snapi)**2


    genpar%verbose=.false.
    genpar%optim=.false.
    memory_needed=0.

    genpar%nthreads=min(genpar%nthreads,size(shotgath))
    call omp_set_num_threads(genpar%nthreads)

    ! image and illumination for each thread
    allocate(imagthread(mod%nz,mod%nx,mod%ny,1,genpar%nthreads)); imagthread=0
    allocate(illuthread(mod%nz,mod%nx,mod%ny,genpar%nthreads)); illuthread=0
    
    imag=0.
    nprocessed=0

    !$OMP PARALLEL DO PRIVATE(i,k,j,dmodgath,wfld_fwd)
    do i=1,size(shotgath)

       call copy_window_vel_gath(mod,modgath(i),genpar,genpargath(i),boundsgath(i),i)
       genpargath(i)%ntsnap=int(genpargath(i)%nt/genpargath(i)%snapi)

       allocate(elevgath(i)%elev(boundsgath(i)%nmin2:boundsgath(i)%nmax2, boundsgath(i)%nmin3:boundsgath(i)%nmax3))
       elevgath(i)%elev=0.

       allocate(dmodgath(size(shotgath(i)%gathtrace)))
       do j=1,size(shotgath(i)%gathtrace)
          allocate(dmodgath(j)%trace(n1,1))
       end do
       dmodgath=shotgath(i)%gathtrace
       do j=1,size(shotgath(i)%gathtrace)
          dmodgath(j)%trace=0.
       end do

       ! Forward: modeling
       call AMOD_to_memory(modgath(i),genpargath(i),dat,boundsgath(i),elevgath(i),dmodgath,sourcegath,wfld_fwd,i) 
       ! Backward: imaging
       call AGRAD_to_memory(modgath(i),genpargath(i),dat,boundsgath(i),elevgath(i),shotgath(i)%gathtrace,wfld_fwd,1)
      ! Copy to final image space
       call mod_copy_image_nparam(modgath(i),imagthread(:,:,:,:,omp_get_thread_num()+1),illuthread(:,:,:,omp_get_thread_num()+1),0)

       ! Deallocate arrays
       deallocate(modgath(i)%imagesmall_nparam)
       deallocate(modgath(i)%illumsmall)
       call deallocateModelSpace_elev(elevgath(i))       
       do k=1,shotgath(i)%ntraces
          call deallocateTraceSpace(dmodgath(k))
       end do
       deallocate(modgath(i)%vel)
       deallocate(dmodgath)
       
       !$OMP CRITICAL
       nprocessed=nprocessed+1
       if (modulo(nprocessed,10).eq.0) write(0,*) 'INFO: Migrated',nprocessed*100/size(shotgath),'% of shots'
       !$OMP END CRITICAL

    end do
    !$OMP END PARALLEL DO

    allocate(illu(mod%nz*mod%nx*mod%ny)); illu=0.
    allocate(imag(mod%nz*mod%nx*mod%ny)); imag=0.
    ! Add all images for each thread to final image
    do i=1,genpar%nthreads

       !$OMP PARALLEL DO PRIVATE(k,j,l)
       do k=1,mod%ny
          do j=1,mod%nx
             do l=1,mod%nz
                imag(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)=imag(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)+imagthread(l,j,k,1,i)               
                illu(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)=illu(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)+illuthread(l,j,k,i)
             end do
          end do
       end do
       !$OMP END PARALLEL DO 
       
    end do

    illu=(illu+maxval(illu)/10000)/sqrt(sum(dprod(illu,illu))/size(illu))
    imag=imag/(illu)

    write(0,*) 'INFO:'
    write(0,*) 'INFO: RTM'
    write(0,*) 'INFO: ---'
    write(0,*) 'INFO: copying image to disk'
    write(0,*) 'INFO:'
    call srite('image',imag,4*size(imag))
    call to_history('n1',mod%nz,'image')
    call to_history('n2',mod%nx,'image')
    call to_history('n3',mod%ny,'image')

    deallocate(imagthread,illu,imag,illuthread)
    deallocate(elevgath)

    stat=0

  end function compute_mig

  function compute_mod() result(stat)
    
    integer            ::                 stat
    
    type(DataSpace)    :: dat

    type(TraceSpace),  dimension(:), allocatable :: resigath
    type(TraceSpace),  dimension(:), allocatable :: dmodgath
    type(ModelSpace_elevation), dimension(:), allocatable :: elevgath
    type(WaveSpace), target                      :: wfld_fwd

    real    :: d1
    integer :: i,j,k,l,n1,begi,endi
    integer :: ntsnap,ntotaltraces,nprocessed
    double precision :: memory_needed,gist

    real, dimension(:), allocatable          :: illu

    call from_aux('coordfile','n2',ntotaltraces)

    allocate(elevgath(size(shotgath)))  ! Each shot has an elevation file
    allocate(resigath(ntotaltraces))

    n1=genpar%nt
    d1=genpar%dt
    genpar%ntsnap=int(genpar%nt/genpar%snapi)
    genpar%dt2=(genpar%dt*genpar%snapi)**2

    do i=1,ntotaltraces
       allocate(resigath(i)%trace(n1,1))
       resigath(i)%trace=0.
    end do
    endi=0
    begi=0

    genpar%verbose=.false.
    genpar%optim=.false.
    memory_needed=0.

    genpar%nthreads=min(genpar%nthreads,size(shotgath))
    call omp_set_num_threads(genpar%nthreads)

    nprocessed=0

    !$OMP PARALLEL DO PRIVATE(i,k,j,dmodgath,begi,endi,wfld_fwd)
    do i=1,size(shotgath)

       call copy_window_vel_gath(mod,modgath(i),genpar,genpargath(i),boundsgath(i),i)
       genpargath(i)%ntsnap=int(genpargath(i)%nt/genpargath(i)%snapi)

       allocate(elevgath(i)%elev(boundsgath(i)%nmin2:boundsgath(i)%nmax2, boundsgath(i)%nmin3:boundsgath(i)%nmax3))
       elevgath(i)%elev=0.

       allocate(dmodgath(shotgath(i)%ntraces))
       do j=1,shotgath(i)%ntraces
          allocate(dmodgath(j)%trace(n1,1))
          dmodgath(j)%coord=shotgath(i)%gathtrace(j)%coord
          dmodgath(j)%dcoord=shotgath(i)%gathtrace(j)%dcoord
          dmodgath(j)%icoord=shotgath(i)%gathtrace(j)%icoord
       end do
       do j=1,shotgath(i)%ntraces
          dmodgath(j)%trace=0.
       end do

       ! Forward: modeling
       call AMOD_to_memory(modgath(i),genpargath(i),dat,boundsgath(i),elevgath(i),dmodgath,sourcegath,wfld_fwd,i) 

       call deallocateWaveSpace(wfld_fwd)
       call deallocateModelSpace_elev(elevgath(i))
       begi=shotgath(i)%begi
       endi=begi+shotgath(i)%ntraces-1
       k=0
       do j=begi,endi
          k=k+1
          resigath(j)%trace=dmodgath(k)%trace
          call deallocateTraceSpace(dmodgath(k))
       end do
       deallocate(modgath(i)%vel)
       deallocate(dmodgath)
       !$OMP CRITICAL
       nprocessed=nprocessed+1
       if (modulo(nprocessed,10).eq.0) write(0,*) 'INFO: Modeled',nprocessed*100/size(shotgath),'% of shots'
       !$OMP END CRITICAL

    end do
    !$OMP END PARALLEL DO

    write(0,*) 'INFO: Done with Processing'
    write(0,*) 'INFO:'

    do i=1,ntotaltraces
       call srite('data',resigath(i)%trace(:,1),4*n1)
    end do

    call to_history('n1',n1,'data')
    call to_history('n2',ntotaltraces,'data')
    call to_history('d1',d1,'data')
    call to_history('d2',1.,'data')
    call to_history('o1',0.,'data')
    call to_history('o2',0.,'data')

    do i=1,ntotaltraces
       call deallocateTraceSpace(resigath(i))
    end do
    deallocate(resigath)
    deallocate(elevgath)

    stat=0

  end function compute_mod

end module compute_function_gradient
