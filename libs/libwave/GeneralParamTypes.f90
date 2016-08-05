module GeneralParam_types
  implicit none

  type GeneralParam
     integer :: lsinc
     integer :: nbound
     real    :: dt
     integer :: ntaper

  end type GeneralParam

  type HigdonParam
     real, allocatable :: gx(:,:,:)
     real, allocatable :: gy(:,:,:)
     real, allocatable :: gz(:,:,:)

  end type HigdonParam
end module GeneralParam_types
