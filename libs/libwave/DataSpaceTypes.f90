module DataSpace_types
  implicit none

  type DataSpace
     real, allocatable :: source(:)
     real, allocatable :: data(:,:,:)
     real, allocatable :: wave(:,:,:,:)

     integer :: it
     integer :: nt
     real    :: dt

     integer :: nx
     integer :: ny

     real    :: ox
     real    :: oy
  end type DataSpace

contains

  subroutine deallocateDataSpace(dat)
    type(DataSpace) :: dat
    if (allocated(dat%source)) deallocate(dat%source)
    if (allocated(dat%data)) deallocate(dat%data)
    if (allocated(dat%wave)) deallocate(dat%wave)
  end subroutine deallocateDataSpace

end module DataSpace_types
