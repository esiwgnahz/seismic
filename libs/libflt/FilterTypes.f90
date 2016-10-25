module Filter_types
  
  use sep
  use helix
  use nhelix
  use ReadData_mod
  use DataSpace_types_flt

  implicit none

  type NSfilter
     type (nfilter)                       ::nmatch ! filter
     integer,       pointer, dimension(:) ::nfilt  ! no. of coefficients in each filter
     integer,       pointer, dimension(:) ::npatch ! no. of patches, dimension(ndim)
     integer,       pointer, dimension(:) ::idx    ! index counter for filter
     integer,       pointer, dimension(:) ::pch    ! ns filter indexing
     integer,       pointer, dimension(:) ::psize  ! non-overlapping patch-size
     integer                              ::ncoef
  end type NSfilter

contains
 
  subroutine NSfilter_deallocate(NSflt)
    type(NSfilter) :: NSflt

    call ndeallocate(NSflt%nmatch)
    
  end subroutine NSfilter_deallocate

  subroutine NSfilter_write_to_file(tag,tagpch,flt,data)
    type(NSfilter)             ::              flt
    character(len=*)           ::   tag       
    character(len=*)           ::       tagpch
    type(cube)                 ::                  data
    integer :: i
    character(len=1024) :: label

    label=''

    call to_history('nfilt',flt%nfilt,tag)
    call to_history('npatch',flt%npatch,tag)
    call to_history('ncoef',flt%ncoef,tag)
    call to_history('psize',flt%psize,tag)
    call to_history('ncoef',flt%ncoef,tag)

    call to_history('lag',flt%nmatch%hlx(1)%lag,tag)

    do i=1,product(flt%npatch)
       call srite(tag,flt%nmatch%hlx(i)%flt,4*flt%ncoef)
    end do

    call to_history('n1',flt%ncoef,tag)
    call to_history('n2',product(flt%npatch),tag)
    call to_history('label1','Filter coefs',tag)
    call to_history('label2','Patches',tag)

    do i=1,size(data%n)
       call sep_put_data_axis_par(tagpch,i,data%n(i),data%o(i),data%d(i),label)
    end do
    
    do i=1,data%n(3)
       call srite(tagpch,real(flt%pch(1+(i-1)*data%n(1)*data%n(2):i*data%n(1)*data%n(2))),4*data%n(1)*data%n(2))
    end do

  end subroutine NSfilter_write_to_file

  subroutine NSfilter_read_param_from_file(tag,flt,ndim)
    type(NSfilter)             ::              flt
    character(len=*)           ::          tag       
    integer :: i,ndim
    character(len=1024) :: label

    label=''

    allocate(flt%npatch(ndim))
    allocate(flt%nfilt(ndim))

    call from_aux(tag,'nfilt',flt%nfilt)
    call from_aux(tag,'npatch',flt%npatch)
    call from_aux(tag,'ncoef',flt%ncoef)

  end subroutine NSfilter_read_param_from_file

  subroutine NSfilter_read_from_file(tag,tagpch,flt,data,ndim)
    type(NSfilter)             ::               flt
    character(len=*)           ::    tag    
    character(len=*)           ::        tagpch  
    type(cube)                 ::                   data,fltpch
    integer :: i,ndim

    call from_aux(tag,'lag',flt%nmatch%hlx(1)%lag)
    call from_aux(tag,'psize',flt%psize)

    do i=1,product(flt%npatch)
       call sreed(tag,flt%nmatch%hlx(i)%flt,4*flt%ncoef)
    end do
    
    call ReadData_dim(tagpch,fltpch,ndim)
    if (.not.hdrs_are_consistent(fltpch,data)) call erexit('ERROR: patches and input file have different sizes, exit now')

    do i=1,data%n(3)
       call sreed(tagpch,flt%pch(1+(i-1)*data%n(1)*data%n(2):i*data%n(1)*data%n(2)),4*data%n(1)*data%n(2))
    end do

    call cube_deallocate(fltpch)

  end subroutine NSfilter_read_from_file

end module Filter_types
