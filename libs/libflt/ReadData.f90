module ReadData_mod

  use sep
  use DataSpace_types_flt

  implicit none

contains

  subroutine InfoData_dim(tag,data,ndim)
    character(len=*)  ::  tag
    type(cube)        ::      data
    integer           ::           ndim
    integer           :: i

    write(0,*) 'INFO:'
    write(0,*) 'INFO: ------------------------------------'
    write(0,*) 'INFO:       Dimensions ',tag,' Data        '
    write(0,*) 'INFO: ------------------------------------'
    write(0,*) 'INFO:'
    do i=1,ndim
       write(0,*) 'INFO: n(',i,')=',data%n(i)
       write(0,*) 'INFO: o(',i,')=',data%o(i)
       write(0,*) 'INFO: d(',i,')=',data%d(i)
    end do
    write(0,*) 'INFO:'
    
  end subroutine InfoData_dim

  subroutine ReadData_dim(tag,data,ndim)
    character(len=*)  ::  tag
    type(cube)        ::      data
    integer           ::           ndim
    integer           :: i
    character(len=1024) :: label

    label=" "
    allocate(data%n(3),data%o(3),data%d(3))

    data%n=1; data%o=0. ; data%d=0.
    do i=1,ndim
       call sep_get_data_axis_par(tag,i,data%n(i),data%o(i),data%d(i),label)
    end do
    
  end subroutine ReadData_dim

  subroutine WriteData_dim(tag,data,ndim)
    character(len=*)  ::   tag
    type(cube)        ::       data
    integer           ::            ndim
    integer           :: i
    character(len=1024) :: label
    label=" "
    do i=1,ndim
       call sep_put_data_axis_par(tag,i,data%n(i),data%o(i),data%d(i),label)
    end do
    
  end subroutine WriteData_dim

  subroutine ReadData_cube(tag,data)
    character(len=*)  ::   tag
    type(cube)        ::       data
    integer           :: i

    allocate(data%dat(data%n(1)*data%n(2)*data%n(3)))

    do i=1,data%n(3)
       call sreed(tag,data%dat(1+(i-1)*data%n(1)*data%n(2):i*data%n(1)*data%n(2)),4*data%n(1)*data%n(2))
    end do
    
  end subroutine ReadData_cube

  subroutine WriteData_cube(tag,data)
    character(len=*)  ::   tag
    type(cube)        ::       data
    integer           :: i

    do i=1,data%n(3)
       call srite(tag,data%dat(1+(i-1)*data%n(1)*data%n(2):i*data%n(1)*data%n(2)),4*data%n(1)*data%n(2))
    end do
    
  end subroutine WriteData_cube

  logical function hdrs_are_consistent(data1,data2)
    type(cube)        ::               data1,data2
    integer           :: i

    hdrs_are_consistent=.true.

    do i=1,3
       if (data1%o(i).ne.data2%o(i)) hdrs_are_consistent=.false.
       if (data1%d(i).ne.data2%d(i)) hdrs_are_consistent=.false.
       if (data1%n(i).ne.data2%n(i)) hdrs_are_consistent=.false.
    end do
   
  end function hdrs_are_consistent

end module ReadData_mod
