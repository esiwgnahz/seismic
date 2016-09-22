module Filter_types
  
  use helix

  implicit none

  type NSfilter
     type (filter), pointer, dimension(:) ::nmatch ! filter
     integer,       pointer, dimension(:) ::nfilt  ! no. of coefficients in each filter
     integer,       pointer, dimension(:) ::npatch ! no. of patches, dimension(ndim)
     integer,       pointer, dimension(:) ::idx    ! index counter for filter
     integer,       pointer, dimension(:) ::pch    ! ns filter indexing
     integer,       pointer, dimension(:) ::psize  ! non-overlapping patch-size
  end type NSfilter

contains

  subroutine NSfilter_deallocate(NSflt)
    type(NSfilter) :: NSflt
    integer        :: i

    if (allocated(NSflt%nfilt)) deallocate(NSflt%nfilt)
    if (allocated(NSflt%npatch)) deallocate(NSflt%npatch)
    if (allocated(NSflt%idx)) deallocate(NSflt%idx)
    if (allocated(NSflt%pch)) deallocate(NSflt%pch)
    if (allocated(NSflt%psize)) deallocate(NSflt%psize)
    do i=1,size(NSflt%nmatch)
       call deallocatehelix(NSflt%nmatch(i))
    end do
    
  end subroutine NSfilter_deallocate
end module Filter_types
