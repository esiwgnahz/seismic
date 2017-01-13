program TWODAFWI

  use sep
  use ReadParam_mod
  use DataSpace_types
  use Readsouvelrho_mod
  use ExtractPadModel_mod

  use ModelSpace_types
  use GeneralParam_types

  use to_disk_to_memory_AMOD_mod
  use to_disk_to_memory_AGRAD_mod

  implicit none

  type(GeneralParam) :: genpar
  type(ModelSpace)   :: mod
  type(DataSpace)    :: dat
  type(FDbounds)     :: bounds
  type(ModelSpace_elevation) :: elev
  type(ModelSpace_elevation), dimension(:), allocatable :: elevgath

  type(ModelSpace),  dimension(:), allocatable :: modgath
  type(GatherSpace), dimension(:), allocatable :: shotgath
  type(TraceSpace),  dimension(:), allocatable :: resigath
  type(TraceSpace),  dimension(:), allocatable :: dmodgath
  type(TraceSpace),  dimension(:), allocatable :: source
  type(TraceSpace),  dimension(:), allocatable :: sourcegath
  type(FDbounds),    dimension(:), allocatable :: boundsgath
  type(GeneralParam),dimension(:), allocatable :: genpargath
  type(WaveSpace), target                      :: wfld_fwd

  real, dimension(:,:,:), allocatable :: illu,grad

  integer :: i,j,k,n1,begi,endi
  integer :: ntsnap,ntotaltraces
  double precision :: memory_needed

  call sep_init()

  mod%veltag='vel'
  mod%waFtag='wave_fwd'

  call from_aux('traces','n2',ntotaltraces)

  allocate(resigath(ntotaltraces));

  if (.not.exist_file(mod%veltag)) call erexit('ERROR: need velocity file')
!  if (.not.exist_file('data')) call erexit('ERROR: need dat file for output')
  call read_3D_params(genpar)
  call readsou(source,genpar)
  call readgathercoords(shotgath,sourcegath,genpar)
  call readgathertraces(shotgath,source)
  call copysou2sougath(source,sourcegath)
  allocate(modgath(size(shotgath)))   ! Each shot has one model
  allocate(boundsgath(size(shotgath)))! Each shot has different bounds
  allocate(genpargath(size(shotgath)))! Each shot has different parameters
  allocate(elevgath(size(shotgath)))  ! Each shot has an elevation file
  do i=1,size(shotgath)
     call extract_coord_source_receiver_patch(shotgath(i)%gathtrace,sourcegath(i),modgath(i),genpar)
  end do
  call read_vel(mod,genpar)

  n1=genpar%nt
  genpar%ntsnap=int(genpar%nt/genpar%snapi)
  do i=1,ntotaltraces
     allocate(resigath(i)%trace(n1,1))
     resigath(i)%trace=0.
  end do
  endi=0
  begi=0

  allocate(grad(mod%nz,mod%nx,mod%ny)); grad=0.
  allocate(illu(mod%nz,mod%nx,mod%ny)); illu=0.

  genpar%verbose=.false.
  genpar%optim=.false.

  do i=1,size(shotgath)
     memory_needed=memory_needed+dble(modgath(i)%nxw)*dble(modgath(i)%nz)*dble(modgath(i)%nyw)*1e-9*4*dble(genpargath(i)%ntsnap)
     write(0,*) 'INFO:'
     write(0,*) 'INFO: Total Memory needed to keep wavefields in memory for all shots=',memory_needed,'Gb'
     write(0,*) 'INFO:'
  end do

  !$OMP PARALLEL DO PRIVATE(i,k,j,dmodgath,begi,endi,wfld_fwd)
  do i=1,size(shotgath)

     !write(0,*) 'INFO: Processing shot',i,'/',size(shotgath)

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

     call AMOD_to_memory(modgath(i),genpargath(i),dat,boundsgath(i),elevgath(i), &
     &                   dmodgath,sourcegath(i),wfld_fwd)   
     call AGRAD_to_memory(modgath(i),genpargath(i),dat,boundsgath(i),elevgath(i),&
     &                    dmodgath,sourcegath(i),wfld_fwd)

     !$OMP CRITICAL
     call mod_copy_image(modgath(i),grad,illu)
     !$OMP END CRITICAL
     
     deallocate(modgath(i)%imagesmall)
     deallocate(modgath(i)%illumsmall)

     call deallocateModelSpace_elev(elevgath(i))
     begi=endi+1
     endi=begi+shotgath(i)%ntraces-1
     k=0
     do j=begi,endi
        k=k+1
        resigath(j)%trace=shotgath(i)%gathtrace(k)%trace-dmodgath(k)%trace
        call deallocateTraceSpace(dmodgath(k))
     end do
     deallocate(dmodgath)

  end do   
  !$OMP END PARALLEL DO

  write(0,*) 'INFO: Done with Processing'
  write(0,*) 'INFO:'

  do i=1,ntotaltraces
     call srite('data',resigath(i)%trace(:,1),4*n1)
  end do

  call to_history('n1',n1,'data')
  call to_history('n2',ntotaltraces,'data')
  call to_history('d1',source(1)%dimt%dt,'data')
  call to_history('d2',1.,'data')
  call to_history('o1',0.,'data')
  call to_history('o2',0.,'data')

  call srite('grad',grad,4*mod%nx*mod%ny*mod%nz)
  call srite('illu',illu,4*mod%nx*mod%ny*mod%nz)

  call to_history('n1',mod%nz,'grad')
  call to_history('n2',mod%nx,'grad')
  call to_history('n3',mod%ny,'grad')
  call to_history('n1',mod%nz,'illu')
  call to_history('n2',mod%nx,'illu')
  call to_history('n3',mod%ny,'illu')


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

!  call deallocateModelSpace(mod)
!
!  write(0,*) 'INFO:'
!  write(0,*) 'INFO: -- Born Modeling End -- '
!  write(0,*) 'INFO:'
!
end program TWODAFWI
