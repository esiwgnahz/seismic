module DataSpace_types
  implicit none

  type DataSpace
     integer :: it
     integer :: nt
     real    :: dt
     real    :: ot
  end type DataSpace

  type WaveSpace
     real, allocatable :: wave(:,:,:,:,:)

     type(DataSpace) :: dimt

     integer :: counter
     integer :: nx
     integer :: ny
     real    :: ox
     real    :: oy

  end type WaveSpace
  
  type TraceSpace
     type(DataSpace)      :: dimt
     real, allocatable    :: trace(:,:) ! value per component
     real,    dimension(3):: coord      ! physical coordinate 
     real,    dimension(3):: dcoord     ! delta for extraction 
     integer, dimension(3):: icoord     ! index
  end type TraceSpace

contains

  subroutine deallocateWaveSpace(wav)
    type(WaveSpace)  :: wav
    if (allocated(wav%wave)) deallocate(wav%wave)
  end subroutine deallocateWaveSpace

  subroutine deallocateTraceSpace(dat)
    type(TraceSpace) :: dat
    if (allocated(dat%trace)) deallocate(dat%trace)
  end subroutine deallocateTraceSpace

end module DataSpace_types
