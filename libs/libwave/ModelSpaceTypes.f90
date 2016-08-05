module ModelSpace_types
  implicit none

  type ModelSpace
     real, allocatable :: vel(:,:,:)
     real, allocatable :: rho(:,:,:)
     real, allocatable :: rho2(:,:,:) 
     real, allocatable :: vel2(:,:,:)      
     real, allocatable :: elev(:,:)     ! Free surface elevation
     real, allocatable :: elev_sou(:,:) ! Source elevation
     real, allocatable :: elev_rec(:,:) ! Receiver elevation
     real, allocatable :: image(:,:,:)  ! Resulting image 
     real, allocatable :: illum(:,:,:)  ! illumination

     integer :: nx
     integer :: ny
     integer :: nz

     real    :: dx
     real    :: dy
     real    :: dz

     real    :: shot_z
     real    :: dshot_z
     real    :: o1model
     integer :: ishot_x
     integer :: ishot_y
     integer :: ishot_z
     
     integer :: surf_type
     integer :: shot_type
     integer :: rcvr_type

  end type ModelSpace

contains
  
  subroutine ModelSpace_compute_shotz_positions(mod)
    type(ModelSpace)  ::                        mod

    !
    ! Shot elevation
    !
    mod%ishot_z = nint( (mod%elev_sou(mod%ishot_x, mod%ishot_y) + mod%shot_z - &
    &  mod%o1model) / mod%dz ) + 1
    mod%dshot_z = mod%elev_sou(mod%ishot_x, mod%ishot_y) + mod%shot_z - &
    &  (mod%o1model + float(mod%ishot_z-1)*mod%dz)

  end subroutine ModelSpace_compute_shotz_positions
  
end module ModelSpace_types
