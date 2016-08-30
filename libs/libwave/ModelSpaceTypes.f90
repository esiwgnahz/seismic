module ModelSpace_types

  use DataSpace_types
  use FD_types

  implicit none

  type ModelSpace

     character(len=3)  :: veltag
     character(len=3)  :: rhotag
     character(len=3)  :: reftag
     character(len=8)  :: waFtag
     character(len=8)  :: waBtag

     real, allocatable :: vel(:,:,:)
     real, allocatable :: rho(:,:,:)
     real, allocatable :: rho2(:,:,:) 
     real, allocatable :: vel2(:,:,:) 
  
     real, allocatable :: image(:,:,:)  ! Resulting image 

     real, allocatable :: illum(:,:,:)  ! illumination

     real, allocatable :: imagesmall(:,:,:)  ! Resulting image 

     real, allocatable :: illumsmall(:,:,:)  ! illumination

     type(WaveSpace), pointer :: wvfld

     integer :: nx  ! dimensions original model space
     integer :: ny
     integer :: nz

     integer :: nxw ! dimensions windowed model space according to
     integer :: nyw ! source/receier/aperture parameters

     real    :: ox  ! origins original model space
     real    :: oy
     real    :: oz

     real    :: oxw ! origins windowed model space according to
     real    :: oyw ! source/receier/aperture parameters

     real    :: dx
     real    :: dy
     real    :: dz

     real    :: endx
     real    :: endy

     integer :: counter

  end type ModelSpace

  type ModelSpace_elevation   
     real,    allocatable :: elev(:,:)     ! Free surface elevation
     integer, allocatable :: ielev_z(:,:)
     real   , allocatable :: delev_z(:,:)
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
    if (allocated(mod%imagesmall)) deallocate(mod%imagesmall)
    if (allocated(mod%illumsmall)) deallocate(mod%illumsmall)
    if (associated(mod%wvfld)) call deallocateWaveSpace(mod%wvfld)
  end subroutine deallocateModelSpace

  subroutine deallocateModelSpace_elev(elev)
    type(ModelSpace_elevation) :: elev
    if (allocated(elev%elev))          deallocate(elev%elev)
    if (allocated(elev%ielev_z))       deallocate(elev%ielev_z)
    if (allocated(elev%delev_z))       deallocate(elev%delev_z)
  end subroutine deallocateModelSpace_elev
  
  subroutine ModelSpace_compute_array_xyz_positions(genpar,vec)
    type(TraceSpace), dimension(:)   ::                   vec
    type(GeneralParam)               ::            genpar
    integer :: i

    do i=1,size(vec)
       call ModelSpace_compute_xyz_positions(genpar,vec(i))
    end do
    
  end subroutine ModelSpace_compute_array_xyz_positions

  subroutine ModelSpace_compute_xyz_positions(genpar,sou)
    type(GeneralParam)  ::                    genpar
    type(TraceSpace)    ::                           sou
    integer :: i
    
    !
    ! Shot elevation
    !   
    do i=1,3
       call find_i_d(sou%icoord(i),sou%dcoord(i),sou%coord(i),genpar%omodel(i),genpar%delta(i))
    end do

  end subroutine ModelSpace_compute_xyz_positions

  subroutine ModelSpace_elevation_parameters(elev,bounds,genpar)
    type(ModelSpace_elevation)  ::           elev
    type(FDbounds)              ::                bounds
    type(GeneralParam)          ::                       genpar
    
    integer :: k,j

    ! Surface elevation parameters
    if (genpar%surf_type.ne.0) then
       do k=bounds%nmin3,bounds%nmax3
          do j=bounds%nmin2,bounds%nmax2
             call find_i_d(elev%ielev_z(j,k),elev%delev_z(j,k),elev%elev(j,k),genpar%omodel(1),genpar%delta(1))
          end do
       end do
    endif
    
  end subroutine ModelSpace_elevation_parameters
  
  subroutine find_i_d(ielev,delev,elev,orig,delta)
    integer ::        ielev
    real    ::              delev,elev,orig,delta

    ielev= nint((elev- orig)/delta)+1           ! Index
    delev=       elev-(orig+float(ielev-1)*delta) ! Differential

  end subroutine find_i_d

end module ModelSpace_types
