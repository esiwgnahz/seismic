! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module DataSpace_types_flt

  implicit none

  type cube
     real,    allocatable, dimension(:) :: d
     real,    allocatable, dimension(:) :: o
     integer, allocatable, dimension(:) :: n

     real,    allocatable, dimension(:) :: dat
  end type cube

contains

  subroutine cube_deallocate(c)
    type(cube) :: c
    if (allocated(c%d)) deallocate(c%d)
    if (allocated(c%n)) deallocate(c%n)
    if (allocated(c%o)) deallocate(c%o)
    if (allocated(c%dat)) deallocate(c%dat)
  end subroutine cube_deallocate
end module DataSpace_types_flt
