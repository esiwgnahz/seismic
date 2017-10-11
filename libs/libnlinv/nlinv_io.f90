module nlinv_io_mod

  implicit none

  use sep

contains

  type sepfile_type
     real, dimension(:), allocatable :: d
     real, dimension(:), allocatable :: o
     real, dimension(:), allocatable :: n

     character(len=*) :: tag
     
     real, dimension(:), allocatable :: array
  end type sepfile_type

  subroutine deallocate_sepfile(file)
    type(sepfile_type) ::       file

    if(allocated(file%array)) deallocate(file%array)
    if(allocated(file%n))     deallocate(file%n)
    if(allocated(file%d))     deallocate(file%d)
    if(allocated(file%o))     deallocate(file%o)
  end subroutine deallocate_sepfile

  ! Read 4*n bytes float array single precision from file tag
  subroutine io_read_sngl(tag,array,n)
    character(len=*)   :: tag
    integer            ::           n
    real, dimension(n) ::     array

    ! SEPLib call
    call sreed(tag,array,4*n)

  end subroutine io_read_sngl

  ! Read 8*n bytes float array double precision from file tag
  subroutine io_read_dble(tag,array,n)
    character(len=*) ::   tag
    integer          ::             n
    double precision, dimension(n) ::     array

    ! SEPLib call
    call sreed(tag,array,8*n)

  end subroutine io_read_dble

  ! Write 4*n bytes float array single precision from file tag
  subroutine io_write_sngl(tag,array,n)
    character(len=*)   ::  tag
    integer            ::            n
    real, dimension(n) ::      array

    ! SEPLib call
    call srite(tag,array,4*n)

  end subroutine io_write_sngl

  ! Write 8*n bytes float array double precision from file tag
  subroutine io_write_dble(tag,array,n)
    character(len=*) ::    tag
    integer         ::               n
    double precision, dimension(n) ::     array

    ! SEPLib call
    call srite(tag,array,8*n)

  end subroutine io_write_dble

  subroutine read_sepfile(file)
    type(sepfile_type) :: file
    integer            :: sep_dim,i

    sep_dim=sep_dimension(file%tag)
    allocate(file%n(sep_dim))
    allocate(file%o(sep_dim))
    allocate(file%d(sep_dim))
    
    do i=1,sep_dim
       call sep_get_data_axis_par(file%tag,i,file%n(i),file%o(i),file%d(i),label)
    end do
    
    allocate(file%array(product(file%n)))
    call sreed(file%tag,file%array,4*product(file%n))
    call auxclose(file%tag)

  end subroutine read_sepfile

  subroutine write_sepfile(file)
    type(sepfile_type) ::  file
    integer            ::  sep_dim,i
    
    do i=1,sep_dim
       call sep_put_data_axis_par(file%tag,i,file%n(i),file%o(i),file%d(i),label)
    end do
    
    call srite(file%tag,file%array,4*product(file%n))
    call auxclose(file%tag)

  end subroutine read_sepfile

end module nlinv_io_mod
