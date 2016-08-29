module io_mod

  use sep

  implicit none

contains

  pure integer function npos(i1, i2, i3, n1, n2)
    integer, intent(in) :: i2, i3, i1, n1, n2
    npos = (i3-1)*n2*n1 + (i2-1)*n1 + (i1-1)
  end function npos

  subroutine iread(tag,array,i1,i2,i3,n1,n2,blocksize)
    integer :: blocksize
    real    :: array(blocksize)
    integer :: i2,i3,i1,n1,n2
    character(len=*) :: tag

    call sseek_block(tag,npos(i1,i2,i3,n1,n2),4*blocksize,0)
    call sreed(tag,array,4*blocksize)
  end subroutine iread

  subroutine iwrite(tag, array, i1, i2, i3, n1, n2, blocksize)
    integer :: blocksize
    real    :: array(blocksize)
    integer :: i2, i3, i1, n1, n2
    character(len=*) :: tag

    call sseek_block(tag, npos(i1,i2,i3,n1,n2), 4*blocksize, 0)
    call srite(tag, array, 4*blocksize)
  end subroutine iwrite

  subroutine iread_smp(tag,array,blocksize)
    integer :: blocksize
    real    :: array(blocksize)   
    character(len=*) :: tag
    call sreed(tag, array, 4*blocksize)
  end subroutine iread_smp

  subroutine iwrite_smp(tag,array,blocksize)
    integer :: blocksize
    real    :: array(blocksize)   
    character(len=*) :: tag
    call srite(tag, array, 4*blocksize)
  end subroutine iwrite_smp

end module io_mod
