module Taper_mod

  use ModelSpace_types

  implicit none

contains

  subroutine compute_taper(mod)
    type(ModelSpace) ::    mod
    
    real    ::pi=3.14159265359

    allocate(mod%taperx(mod%nxw))
    allocate(mod%tapery(mod%nyw))

    mod%taperx=1.
    mod%tapery=1.
    
    call cos_sin_taper(mod%taperx,min(mod%nxw/2,10),mod%nxw)
    if (mod%nyw.gt.1) call cos_sin_taper(mod%tapery,min(mod%nyw/2,10),mod%nyw)

  end subroutine compute_taper

  subroutine cos_sin_taper(array,ntaper,n)
    real, dimension(:) ::  array
    integer ::                   ntaper
    integer ::                          n
    integer :: i
    real    :: pi=3.14159265359

    do i=1,min(ntaper,n)
       array(i)=sin(pi*(i-1)/(2*(min(ntaper,n)-1)))**2
    end do
    do i=max(1,n-ntaper),n
       array(i)=cos(pi*(i-max(1,n-ntaper))/(2*(ntaper)))**2
    end do
    
  end subroutine cos_sin_taper

  subroutine compute_taper_close(mod)    
    type(ModelSpace) ::    mod
    
    if(allocated(mod%taperx)) deallocate(mod%taperx)   
    if(allocated(mod%tapery)) deallocate(mod%tapery)
  end subroutine compute_taper_close

end module Taper_mod
