module Imaging_mod
  
  use sep

  use GeneralParam_types
  use ModelSpace_types
  use DataSpace_types
  use Interpolate_mod

  use ExtractPadModel_mod

  implicit none

contains

  subroutine Imaging_condition_from_disk(model,genpar)
    type(ModelSpace) :: model
    type(GeneralParam) :: genpar
    integer :: i,j,k,l,counter
    real, dimension (:,:,:), allocatable :: fwd,bwd,tmpim,tmpil
 
    integer :: ierr,blocksize,index

    allocate(tmpim(model%nz,model%nxw,model%nyw))
    allocate(tmpil(model%nz,model%nxw,model%nyw))

    tmpim=0.
    tmpil=0.

    allocate(fwd(model%nz,model%nxw,model%nyw))
    allocate(bwd(model%nz,model%nxw,model%nyw))

    blocksize=model%nz*model%nxw*model%nyw
    call sseek(model%waFtag,0,0)
    call sseek(model%waBtag,0,0)
    do l=1,genpar%ntsnap

       fwd=0.
       bwd=0.

       index=genpar%ntsnap-l
       call sseek_block(model%waFtag,index,4*blocksize,0)
       call sreed(model%waBtag,bwd,4*blocksize)
       call sreed(model%waFtag,fwd,4*blocksize)
       
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=1,model%nyw
          do j=1,model%nxw
             do i=1,model%nz
                tmpim(i,j,k)=tmpim(i,j,k)+fwd(i,j,k)*bwd(i,j,k)
                tmpil(i,j,k)=tmpil(i,j,k)+fwd(i,j,k)*fwd(i,j,k)
             end do
          end do
       end do
       !$OMP END PARALLEL DO

    end do
    deallocate(fwd,bwd)
    
    allocate(model%image(model%nz,model%nx,model%ny))
    allocate(model%illum(model%nz,model%nx,model%ny))
    
    model%image=0.
    model%illum=0.

    call mod_window_pad(.false.,model%image,tmpim,model)
    call mod_window_pad(.false.,model%illum,tmpil,model)

    deallocate(tmpim,tmpil)

  end subroutine Imaging_condition_from_disk

end module Imaging_mod
