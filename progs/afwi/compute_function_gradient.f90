module compute_function_gradient

  use sep
  use omp_lib
  use ReadParam_mod
  use DataSpace_types
  use Readsouvelrho_mod
  use ExtractPadModel_mod

  use Mute_gather
  use ModelSpace_types
  use GeneralParam_types

  use to_disk_to_memory_AMOD_mod
  use to_disk_to_memory_AGRAD_mod

  implicit none

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

  subroutine compute_fct_gft_dmute_init(mutepar_in)
    type(MuteParam), target  ::         mutepar_in

    mutepar=>mutepar_in
  end subroutine compute_fct_gft_dmute_init

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
    if(associated(mutepar))    nullify(mutepar)
    if(associated(genpar))     nullify(genpar)
    if(associated(modgath))    nullify(modgath)
    if(associated(mod)) nullify(mod)
    if(associated(bounds))     nullify(bounds)
    if(associated(shotgath))   nullify(shotgath)
    if(associated(sourcegath)) nullify(sourcegath)
    if(associated(boundsgath)) nullify(boundsgath)
    if(associated(genpargath)) nullify(genpargath)
  end subroutine compute_fct_gdt_nullify

  function compute_fct_gdt(grad,f) result(stat)
    
    real, dimension(:) ::  grad
    double precision   ::       f
    integer            ::                 stat
    
    type(DataSpace)    :: dat

    type(TraceSpace),  dimension(:), allocatable :: resigath
    type(TraceSpace),  dimension(:), allocatable :: dmodgath
    type(ModelSpace_elevation), dimension(:), allocatable :: elevgath
    type(WaveSpace), target                      :: wfld_fwd

    real, dimension(:,:,:,:), allocatable :: illuthread,gradthread
    double precision, dimension(:),       allocatable :: fthread

    real    :: d1
    integer :: i,j,k,l,n1,begi,endi
    integer :: ntsnap,ntotaltraces
    double precision :: memory_needed,gist

    real, dimension(:), allocatable          :: illu

    call from_aux('traces','n2',ntotaltraces)

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

    ! image and illumination for each thread
    allocate(gradthread(mod%nz,mod%nx,mod%ny,genpar%nthreads)); gradthread=0
    allocate(illuthread(mod%nz,mod%nx,mod%ny,genpar%nthreads)); illuthread=0
    allocate(fthread(genpar%nthreads));fthread=0.
    
    f=0.
    grad=0.

    !$OMP PARALLEL DO PRIVATE(i,k,j,dmodgath,begi,endi,wfld_fwd)
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

       ! Compute residual and of
       do j=1,size(shotgath(i)%gathtrace)
          ! rd=W(L(m)-d)
          dmodgath(j)%trace=mutepar%maskgath(i)%gathtrace(j)%trace*(shotgath(i)%gathtrace(j)%trace-dmodgath(j)%trace)
!          dmodgath(j)%trace=(shotgath(i)%gathtrace(j)%trace-dmodgath(j)%trace)
          ! f=rd'rd
          fthread(omp_get_thread_num ()+1)=fthread(omp_get_thread_num ()+1)+sum(dprod(dmodgath(j)%trace,dmodgath(j)%trace))
          ! rd=W'W(L(m)-d)
          dmodgath(j)%trace=mutepar%maskgath(i)%gathtrace(j)%trace*dmodgath(j)%trace
       end do

       ! Backward: imaging
       call AGRAD_to_memory(modgath(i),genpargath(i),dat,boundsgath(i),elevgath(i),dmodgath,wfld_fwd)

       ! Copy to final image space 
       call mod_copy_image(modgath(i),gradthread(:,:,:,omp_get_thread_num()+1),illuthread(:,:,:,omp_get_thread_num()+1))     
       deallocate(modgath(i)%imagesmall)
       deallocate(modgath(i)%illumsmall)
       !
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

    end do
    !$OMP END PARALLEL DO

!    write(0,*) 'INFO: Done with Processing'
!    write(0,*) 'INFO:'

    do i=1,ntotaltraces
       call srite('data',resigath(i)%trace(:,1),4*n1)
    end do

    call to_history('n1',n1,'data')
    call to_history('n2',ntotaltraces,'data')
    call to_history('d1',d1,'data')
    call to_history('d2',1.,'data')
    call to_history('o1',0.,'data')
    call to_history('o2',0.,'data')

    allocate(illu(mod%nz*mod%nx*mod%ny)); illu=0.
    ! Add all images for each thread to final image
    do i=1,genpar%nthreads
       f   =f   +fthread(i)

       !$OMP PARALLEL DO PRIVATE(k,j,l)
       do k=1,mod%ny
          do j=1,mod%nx
             do l=1,mod%nz
                grad(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)=grad(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)+gradthread(l,j,k,i)               
                illu(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)=illu(l+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx)+illuthread(l,j,k,i)
             end do
          end do
       end do
       !$OMP END PARALLEL DO 

    end do

    ! Scaling
    f=f/(2*n1*ntotaltraces)
    grad=grad/(n1*ntotaltraces)
    illu=(illu+maxval(illu)/10000)/sqrt(sum(dprod(illu,illu))/size(illu))

    grad=grad/illu

    call srite('grad',grad,4*mod%nx*mod%ny*mod%nz)
    deallocate(gradthread,illu,illuthread,fthread)

    do i=1,ntotaltraces
       call deallocateTraceSpace(resigath(i))
    end do
    deallocate(resigath)
    deallocate(elevgath)

    stat=0

  end function compute_fct_gdt

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
