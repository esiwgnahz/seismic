module BuildAdaptiveFilter_mod
  
  use ncnbound  
  use cartesian
  use lapfac3_mod
  use lapfac_mod
  use Filter_types
  use GenParam_types_flt
  use DataSpace_types_flt

  implicit none

contains

  !! Calculate non-overlapping patch-size given a certain number of patches
  subroutine psize_init(dat,ndim,nsmatch)
    type(cube)      ::  dat
    integer         ::      ndim
    type(NSfilter)  ::            nsmatch
    integer                            :: i
    allocate(nsmatch%psize(ndim))
    do i=1,ndim
       nsmatch%psize(i)=1+(dat%n(i)-1)/nsmatch%npatch(i)
    end do
  end subroutine psize_init

  !! Divides data-space between rectangular patches
  subroutine pch_init(dat,nsmatch)
    type(cube)      ::dat
    type(NSfilter)  ::    nsmatch

    integer, dimension(:), pointer     :: ii
    integer :: i,j

    allocate(nsmatch%pch(product(dat%n)),ii(size(dat%n)))

    do i=1,product(dat%n)
       nsmatch%pch(i)=0
       call line2cart(dat%n,i,ii)
       do j=1,size(dat%n)
          ii(j)=(1+(ii(j)-1)/nsmatch%psize(j))
       end do
       call cart2line(nsmatch%npatch,ii,nsmatch%pch(i))
    end do
    deallocate(ii)

    write(0,*) 'INFO:'
    write(0,*) 'INFO: ------------------------------------'
    write(0,*) 'INFO:        NS filters parameters        '
    write(0,*) 'INFO: ------------------------------------'
    write(0,*) 'INFO: Number of patches in dim 1 =',nsmatch%npatch(1)
    write(0,*) 'INFO: Number of patches in dim 2 =',nsmatch%npatch(2)
    write(0,*) 'INFO: Number of patches in dim 3 =',nsmatch%npatch(3)
    write(0,*) 'INFO: - - - - - - - -'
    write(0,*) 'INFO: Size patches in dim 1 =',nsmatch%psize(1)
    write(0,*) 'INFO: Size patches in dim 2 =',nsmatch%psize(2)
    write(0,*) 'INFO: Size patches in dim 3 =',nsmatch%psize(3)
    write(0,*) 'INFO: - - - - - - - -'
    write(0,*) 'INFO: Size NS filters in dim 1  =',nsmatch%nfilt(1)
    write(0,*) 'INFO: Size NS filters in dim 2  =',nsmatch%nfilt(2)
    write(0,*) 'INFO: Size NS filters in dim 3  =',nsmatch%nfilt(3)
    write(0,*) 'INFO: ------------------------------------'
    write(0,*) 'INFO:'

  end subroutine pch_init

  subroutine create_match_filter(nfilt,filt,ndim,n)
    type(filter)         ::            filt
    integer              ::                 ndim
    integer, dimension(:)::      nfilt,          n

    integer, dimension(:), allocatable :: idx
    integer                        :: ilag,idim

    allocate(idx(ndim))
    call allocatehelix(filt,product(nfilt))
    do ilag=1,product(nfilt)
       call line2cart(nfilt,ilag,idx)
       filt%lag(ilag)=idx(1)-nfilt(1)/2-1
       do idim=2,ndim
          filt%lag(ilag)=filt%lag(ilag)+(idx(idim)-nfilt(idim)/2-1)*product(n(1:idim-1))
       end do
    end do
    deallocate(idx)
    
  end subroutine create_match_filter

  subroutine create_nsmatch_filter(dat,ndim,nsmatch)
    type(cube)      ::             dat
    integer         ::                 ndim
    type(NSfilter)  ::                       nsmatch

    type(filter)    ::tmp
    integer         ::ipatch

    integer, dimension(:), allocatable :: nh

    allocate(nh(product(nsmatch%npatch)))           !! nh is the length of each filter in table
    nh=product(nsmatch%nfilt)                       !! Set them all constant
    
    call nallocate(nsmatch%nmatch,nh,nsmatch%pch)   
    call create_match_filter(nsmatch%nfilt,tmp,ndim,dat%n)
    
    do ipatch=1,product(nsmatch%npatch)
       nsmatch%nmatch%hlx(ipatch)%lag=tmp%lag
    end do

    nsmatch%ncoef=size(tmp%lag)

    deallocate(nh)
    call deallocatehelix(tmp)
    call ncnboundn(1,dat%n,nsmatch%nfilt,nsmatch)

  end subroutine create_nsmatch_filter

  subroutine create_lap_3d(rough,nlaplac)
    
    type(NSfilter)       :: rough
    integer              :: npts(3),i,j,nlaplac
    integer, allocatable :: nh(:),map(:)
    
    npts=(/rough%npatch(1),rough%npatch(2),rough%npatch(3)/)
    allocate(nh(product(npts)))
    allocate(map(rough%ncoef*product(npts)))
    
    nh=rough%ncoef
    map=1 ! This means that we have one patch only for all coefficients
    
    nh=min(nlaplac,max(1,npts(1)/2))
    call nallocate(rough%nmatch,nh,map) 

    if (npts(1).ge.2) then
       call nallocate(rough%nmatch,nh,map) 
       ! we only compute the filter for one patch=the whole thing
       rough%nmatch%hlx(1)=lapfac3d_mod(.001,npts(1),npts(2),nh(1))
       ! bad 2D smoothing
!       call allocatehelix(rough%nmatch%hlx(1),3)
!       rough%nmatch%hlx(1)%lag(1)=1
!       rough%nmatch%hlx(1)%lag(2)=npts(2)
!       rough%nmatch%hlx(1)%lag(3)=npts(1)*npts(2)
!       rough%nmatch%hlx(1)%flt(1)=-0.9
!       rough%nmatch%hlx(1)%flt(2)=-0.8
!       rough%nmatch%hlx(1)%flt(3)=-0.7
    else
       ! bad 2D smoothing
       call allocatehelix(rough%nmatch%hlx(1),2)
       rough%nmatch%hlx(1)%lag(1)=1
       rough%nmatch%hlx(1)%lag(2)=npts(2)
       rough%nmatch%hlx(1)%flt(1)=-0.9
       rough%nmatch%hlx(1)%flt(2)=-0.8
    end if
    rough%nmatch%hlx(1)%lag=rough%nmatch%hlx(1)%lag*rough%ncoef
       
    write(0,*) 'INFO: Regularization filter parameters:'
    do i=1,size(rough%nmatch%hlx(1)%flt)
       write(0,*) 'INFO: Filter ',i,'=',rough%nmatch%hlx(1)%flt(i),' at lag ',rough%nmatch%hlx(1)%lag(i)
    end do
    write(0,*) 'INFO:'

    deallocate(nh,map)
    
  end subroutine create_lap_3d
 
end module BuildAdaptiveFilter_mod
