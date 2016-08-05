module DataSpace_types
  implicit none

  type DataSpace
     real, allocatable :: source(:)
     real, allocatable :: data(:,:,:)

     integer :: it
     integer :: nt
     integer :: dt
  end type DataSpace

end module DataSpace_types
