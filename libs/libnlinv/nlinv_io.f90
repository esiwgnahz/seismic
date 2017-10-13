module nlinv_io_mod

  use sep
  implicit none

  type sepfile_type
     real,    dimension(:), allocatable :: d
     real,    dimension(:), allocatable :: o
     integer, dimension(:), allocatable :: n

     integer :: ndim

     character(len=1024) :: tag
     
     real, dimension(:), allocatable :: array
  end type sepfile_type

contains

  subroutine deallocate_sepfile(file)
    type(sepfile_type) ::       file

    if(allocated(file%array)) deallocate(file%array)
    if(allocated(file%n))     deallocate(file%n)
    if(allocated(file%d))     deallocate(file%d)
    if(allocated(file%o))     deallocate(file%o)
  end subroutine deallocate_sepfile

  subroutine read_sepfile(file)
    type(sepfile_type) :: file
    integer            :: i
    character(len=1024):: label

    label=" "

    write(1010,*) '  '
    write(1010,"(A13,A10)") 'Reading file ',file%tag

    file%ndim=sep_dimension(file%tag)
    allocate(file%n(file%ndim))
    allocate(file%o(file%ndim))
    allocate(file%d(file%ndim))
    
    write(1010,*) '-- Dimensions=',file%ndim

    do i=1,file%ndim
       call sep_get_data_axis_par(file%tag,i,file%n(i),file%o(i),file%d(i),label)
       write(1010,*) '-- axis(',i,'),n,o,d=',file%n(i),file%o(i),file%d(i)
    end do
    
    if (.not.allocated(file%array)) allocate(file%array(product(file%n)))
    call sreed(file%tag,file%array,4*product(file%n))
    call auxclose(file%tag)

  end subroutine read_sepfile

  subroutine write_sepfile(file)
    type(sepfile_type) ::  file
    integer            ::  sep_dim,i
    character(len=1024)::  label
    
    label=" "

    do i=1,file%ndim
       call sep_put_data_axis_par(file%tag,i,file%n(i),file%o(i),file%d(i),label)
    end do
    
    call srite(file%tag,file%array,4*product(file%n))
    call auxclose(file%tag)

  end subroutine write_sepfile

end module nlinv_io_mod
