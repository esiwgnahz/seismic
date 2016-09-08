module Laplace_mod

  implicit none

contains
  
  !----------------------------------------------------------------
  subroutine Laplace(nx, ny, nz, n4, n5, dx, dy, dz)
    !
    integer :: nx, ny, nz, n4, n5
    integer :: pad, pady
    integer :: i, j, k, kk, ii, jj, k4, k5
    real :: c(-1:1), dx, dy, dz
    real :: dzz, dxx, dyy, fdum, f3, f2
    real, allocatable :: value(:,:,:), slice(:,:), fscale(:,:,:)
    !
    pad = 1
    if (ny .gt. 1) then
       pady = pad
    else
       pady = 0
    end if
    !
    allocate(value(1-pad:nz+pad, 1-pad:nx+pad, 1-pady:1+pady))
    allocate(fscale(-pad:pad, -pad:pad, -pady:pady))
    allocate(slice(nz, nx))
    !
    ! Initialize
    !
    if (ny .le. 1 ) then
       kk=0
       do jj=-pad,pad
          do ii=-pad,pad
             fscale(ii,jj,kk) = -1.
          end do
       end do
       fscale(0,0,kk) = 8.
    else
       do ii=-1,1,2
          fscale(ii,-1,-1) = -1.
          fscale(ii,-1,1) = -1.
          fscale(ii,1,-1) = -1.
          fscale(ii,1,1) = -1.
          fscale(ii,0,0) = -1.
          fscale(ii,0,-1) = 0.
          fscale(ii,0,1) = 0.
          fscale(ii,1,0) = 0.
          fscale(ii,-1,0) = 0.
       end do
       ii=0
       fscale(ii,-1,-1) = 0.
       fscale(ii,-1,1) = 0.
       fscale(ii,1,-1) = 0.
       fscale(ii,1,1) = 0.
       fscale(ii,0,0) = 14.
       fscale(ii,0,-1) = -1.
       fscale(ii,0,1) = -1.
       fscale(ii,1,0) = -1.
       fscale(ii,-1,0) = -1.
    end if
    !
    !
    !
    do k5=1,n5
       do k4=1,n4
          value = 0.
          !if (n4.gt.1) write(0,*) k4*k5,' out of ',n4*n5
          if (n4.gt.1 .and. mod(k4-1,max(1,n4/5)).eq.0) &
          write(0,*) '  Percentage completed ', &
          &  nint(10.0*float(k4)/float(n4))*10
          !
          ! Read first slab
          !
          do k=1,1+pady
             call sreed('in',slice,4*nz*nx)
             do j=1,nx
                do i=1,nz
                   value(i,j,k) = slice(i,j)
                end do
             end do
             call Padding(value, nx, ny, nz, pad, pad, pady, k)
          end do
          if (pady .gt. 0) then
             do k=1-pady,0
                do j=1-pad,nx+pad
                   do i=1-pad,nz+pad
                      value(i,j,k) = value(i,j,2-k)
                   end do
                end do
             end do
          end if
          !
          if (ny .gt. 1) then
             do k=1,ny
                do j=1,nx
                   do i=1,nz
                      fdum = 0.
                      do kk=-pady,pady
                         do jj=-pad,pad
                            do ii=-pad,pad
                               fdum = fdum + value(i+ii,j+jj,1+kk) * &
                               & fscale(ii,jj,kk)
                            end do
                         end do
                      end do
                      slice(i,j) = fdum
                   end do
                end do
                !if (n4.le.1) write(0,*) k,' out of ',ny
                if (mod(k-1,max(1,ny/5)).eq.0) &
                write(0,*) '  Percentage completed ', &
                &  nint(10.0*float(k)/float(ny))*10

                call srite('out',slice,4*nz*nx)
                !
                ! Read next slice
                !
                do kk=1-pady,pady
                   do j=1-pad,nx+pad
                      do i=1-pad,nz+pad
                         value(i,j,kk) = value(i,j,kk+1)
                      end do
                   end do
                end do
                if (k+pady .lt. ny) then
                   call sreed('in',slice,4*nz*nx)
                   kk = 1+pady
                   do j=1,nx
                      do i=1,nz
                         value(i,j,kk) = slice(i,j)
                      end do
                   end do
                   call Padding(value, nx, ny, nz, pad, pad, pady, kk)
                else if (k .lt. ny) then
                   kk = 1+pady
                   do j=1-pad,nx+pad
                      do i=1-pad,nz+pad
                         value(i,j,kk) = value(i,j,ny-k)
                      end do
                   end do
                end if
             end do
             !
          else
             kk = 1
             do j=1,nx
                do i=1,nz
                   fdum = 0.
                   do jj=-pad,pad
                      do ii=-pad,pad
                         fdum = fdum + value(i+ii,j+jj,kk) * &
                         & fscale(ii,jj,kk-1)
                      end do
                   end do
                   slice(i,j) = fdum
                end do
             end do
             call srite('out',slice,4*nz*nx)
          end if
       end do
    end do
    !
    deallocate(value, slice, fscale)

  end subroutine Laplace

  subroutine Padding(value, nx, ny, nz, rect1, rect2, rect3, k)
    !
    implicit none
    integer :: rect1, rect2, rect3
    integer :: nx, ny, nz, i, j, k
    real :: value(1-rect1:nz+rect1, 1-rect2:nx+rect2, 1-rect3:1+rect3)
    !
    do j=1,nx
       do i=1-rect1,0
          value(i,j,k) = value(2-i,j,k)
       end do
       do i=nz+1,nz+rect1
          value(i,j,k) = value(2*nz-i,j,k)
       end do
    end do
    do j=1-rect2,0
       do i=1-rect1,nz+rect1
          value(i,j,k) = value(i,2-j,k)
       end do
    end do
    do j=nx+1,nx+rect2
       do i=1-rect1,nz+rect1
          value(i,j,k) = value(i,2*nx-j,k)
       end do
    end do
    !
  end subroutine Padding

end module Laplace_mod
