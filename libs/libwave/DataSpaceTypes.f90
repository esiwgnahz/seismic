! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module DataSpace_types

  use FD_types
  use GeneralParam_types

  implicit none

  type DataSpace
     integer :: it
     integer :: nt
     real    :: dt
     real    :: ot
  end type DataSpace

  ! Wavespace type for x,y,z,t,nc (# components)
  type WaveSpace
     real, allocatable :: wave(:,:,:,:,:)

     type(DataSpace) :: dimt

     integer :: counter
     integer :: nx
     integer :: ny
     real    :: ox
     real    :: oy

  end type WaveSpace
  
  ! Wavefields needed for FD schemes. Time index given by number
  type USpace
     real, pointer :: utmp(:,:,:)
     real, pointer :: u_1(:,:,:)
     real, pointer :: u0(:,:,:)
     real, pointer :: u1(:,:,:)
     real, pointer :: u2(:,:,:)
     real, pointer :: u3(:,:,:)
  end type USpace

  ! A collection of traces with their values and coordinates
  type TraceSpace
     type(DataSpace)      :: dimt
     real, allocatable    :: trace(:,:) ! value per component
     real,    dimension(3):: coord      ! physical coordinate 
     real,    dimension(3):: dcoord     ! delta for extraction 
     integer, dimension(3):: icoord     ! index
  end type TraceSpace

  ! A collection of traces belonging to the same gather
  type GatherSpace
     type(TraceSpace), allocatable :: gathtrace(:)
     integer          :: ntraces
     integer          :: begi ! Index of first trace in trace file
  end type GatherSpace

contains

  subroutine deallocateWaveSpace(wav)
    type(WaveSpace)  :: wav
    if (allocated(wav%wave)) deallocate(wav%wave)
  end subroutine deallocateWaveSpace

  subroutine deallocateTraceSpace(dat)
    type(TraceSpace) :: dat
    if (allocated(dat%trace)) deallocate(dat%trace)
  end subroutine deallocateTraceSpace

  subroutine deallocateGatherSpace(gath)
    type(GatherSpace)  :: gath
    integer :: nshots,ntraces,i,j

    ntraces=gath%ntraces
    do j=1,ntraces
       call deallocateTraceSpace(gath%gathtrace(j))
    end do

  end subroutine deallocateGatherSpace

  subroutine deallocateUSpace(grid)
    type(USpace) :: grid
    if (allocated(grid%u_1))  deallocate(grid%u_1)
    if (allocated(grid%u0))   deallocate(grid%u0)
    if (allocated(grid%u1))   deallocate(grid%u1)
    if (allocated(grid%u2))   deallocate(grid%u2)
    if (allocated(grid%u3))   deallocate(grid%u3)
    if (associated(grid%utmp)) nullify(grid%utmp)
  end subroutine deallocateUSpace

  subroutine allocateUSpace(grid,genpar,bounds)
    type(USpace)      ::    grid
    type(GeneralParam)::         genpar
    type(FDbounds)    ::                bounds

    allocate(grid%u_1(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))   
    allocate( grid%u0(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))    
    allocate( grid%u1(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))    
    allocate( grid%u2(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))    
    allocate( grid%u3(bounds%nmin1-4:bounds%nmax1+4, bounds%nmin2-4:bounds%nmax2+4,bounds%nmin3-genpar%nbound:bounds%nmax3+genpar%nbound))   

  end subroutine allocateUSpace

  subroutine zeroUSpace(grid)
    type(USpace)  ::    grid
    grid%u_1=0.
    grid%u0=0.
    grid%u1=0.
    grid%u2=0.
    grid%u3=0.
  end subroutine zeroUSpace

end module DataSpace_types
