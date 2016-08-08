module ModelSpace_types

  use FD_types

  implicit none

  type ModelSpace
     real, allocatable :: vel(:,:,:)
     real, allocatable :: rho(:,:,:)
     real, allocatable :: rho2(:,:,:) 
     real, allocatable :: vel2(:,:,:)   
     real, allocatable :: image(:,:,:)  ! Resulting image 
     real, allocatable :: illum(:,:,:)  ! illumination

     integer :: nx
     integer :: ny
     integer :: nz

     real    :: dx
     real    :: dy
     real    :: dz
  end type ModelSpace

  type ModelSpace_elevation   
     real,    allocatable :: elev(:,:)     ! Free surface elevation
     real,    allocatable :: elev_sou(:,:) ! Source elevation
     real,    allocatable :: elev_rec(:,:) ! Receiver elevation
     integer, allocatable :: irec_z(:,:)
     integer, allocatable :: ielev_z(:,:)
     real,    allocatable :: drec_z(:,:)
     real   , allocatable :: delev_z(:,:)
     real    :: rec_z
     real    :: shot_z
     real    :: dshot_z
     real    :: o1model
     integer :: ishot_x
     integer :: ishot_y
     integer :: ishot_z
     real    :: dz
  end type ModelSpace_elevation

contains
  
  subroutine deallocateModelSpace(mod)
    type(ModelSpace) :: mod
    if (allocated(mod%vel))   deallocate(mod%vel)
    if (allocated(mod%vel2))  deallocate(mod%vel2)
    if (allocated(mod%rho))   deallocate(mod%rho)
    if (allocated(mod%rho2))  deallocate(mod%rho2)
    if (allocated(mod%image)) deallocate(mod%image)
    if (allocated(mod%illum)) deallocate(mod%illum)
  end subroutine deallocateModelSpace

  subroutine deallocateModelSpace_elev(elev)
    type(ModelSpace_elevation) :: elev
    if (allocated(elev%elev))     deallocate(elev%elev)
    if (allocated(elev%elev_sou)) deallocate(elev%elev_sou)
    if (allocated(elev%elev_rec)) deallocate(elev%elev_rec)
    if (allocated(elev%ielev_z))  deallocate(elev%ielev_z)
    if (allocated(elev%irec_Z))   deallocate(elev%irec_z)
    if (allocated(elev%drec_z))   deallocate(elev%drec_z)
    if (allocated(elev%delev_z))  deallocate(elev%delev_z)
  end subroutine deallocateModelSpace_elev
  
  subroutine ModelSpace_compute_shotz_positions(elev)
    type(ModelSpace_elevation)  ::              elev

    !
    ! Shot elevation
    !
    elev%ishot_z = nint( (elev%elev_sou(elev%ishot_x, elev%ishot_y) + elev%shot_z - &
    &  elev%o1model) / elev%dz ) + 1
    elev%dshot_z = elev%elev_sou(elev%ishot_x, elev%ishot_y) + elev%shot_z - &
    &  (elev%o1model + float(elev%ishot_z-1)*elev%dz)

  end subroutine ModelSpace_compute_shotz_positions

  subroutine ModelSpace_elevation_parameters(elev,bounds,genpar)
    type(ModelSpace_elevation)  ::           elev
    type(FDbounds)              ::                bounds
    type(GeneralParam)          ::                       genpar
    
    integer :: k,j

    do k=bounds%nmin3,bounds%nmax3
       do j=bounds%nmin2,bounds%nmax2
          elev%irec_z(j,k) = nint( (elev%elev_rec(j,k) + elev%rec_z - &
          &  elev%o1model) / elev%dz ) + 1
          elev%drec_z(j,k) = elev%elev_rec(j,k) + elev%rec_z - &
          &  (elev%o1model + float(elev%irec_z(j,k)-1)*elev%dz)
       end do
    end do
    if (genpar%surf_type.ne.0) then
       do k=bounds%nmin3,bounds%nmax3
          do j=bounds%nmin2,bounds%nmax2
             elev%ielev_z(j,k) = int( (elev%elev(j,k) - &
             &  elev%o1model) / elev%dz ) + 1
             elev%delev_z(j,k) = elev%elev(j,k) - &
             &  (elev%o1model + float(elev%ielev_z(j,k)-1)*elev%dz)
          end do
       end do
    endif
    
  end subroutine ModelSpace_elevation_parameters
  
end module ModelSpace_types
