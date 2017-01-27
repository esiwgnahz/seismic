module Bandpass_mod

  use sep
  use DataSpace_types

  implicit none

  type BandPassParam
     real :: fhi
     real :: flo
     integer :: phase
     integer :: nphi
     integer :: nplo
     
     real    :: dt
     real    :: nt
  end type BandPassParam

contains

  subroutine Init_BandPassParam(bpparam)
    type(BandPassParam) ::      bpparam

    call from_param('fhi',bpparam%fhi,20.)
    call from_param('flo',bpparam%flo,0.)
    call from_param('nphi',bpparam%nphi,14)
    call from_param('nplo',bpparam%nplo,14)
    call from_param('phase',bpparam%phase,0)

    write(0,*) 'INFO:--------------------------'
    write(0,*) 'INFO: Data bandpass parameters '
    write(0,*) 'INFO:--------------------------'
    write(0,*) 'INFO:'
    write(0,*) 'INFO:  fhi    = ',bpparam%fhi
    write(0,*) 'INFO:  flo    = ',bpparam%flo
    write(0,*) 'INFO:  nphi   = ',bpparam%nphi
    write(0,*) 'INFO:  nplo   = ',bpparam%nplo
    write(0,*) 'INFO:  phase  = ',bpparam%phase
    write(0,*) 'INFO:'
    write(0,*) 'INFO:-------------------------'

    
  end subroutine Init_BandPassParam

  subroutine BandpassSouTraces(bpparam,shotgath,sourcegath)
    type(BandPassParam)             :: bpparam
    type(GatherSpace), dimension(:) :: shotgath
    type(TraceSpace),  dimension(:) :: sourcegath

    real, dimension(:), allocatable :: srctmp,tracetmp
    integer :: nshots,i,j,nt,ntraces

    nshots=size(shotgath)
    allocate(srctmp(sourcegath(1)%dimt%nt),tracetmp(sourcegath(1)%dimt%nt))
    bpparam%dt=sourcegath(1)%dimt%dt
    bpparam%nt=sourcegath(1)%dimt%nt

    write(0,*) 'INFO:----------------------------------------------'
    write(0,*) 'INFO: Starting bandassing'
    do i=1,nshots

       srctmp=0.
       call bandpass(1,sourcegath(i)%trace(:,1),srctmp,bpparam)
       sourcegath(i)%trace(:,1)=srctmp

       ntraces=shotgath(i)%ntraces
       if (modulo(i,10).eq.0) write(0,*) 'INFO: Starting bandassing',ntraces,'traces for shot',i

       do j=1,ntraces
          tracetmp=0.
          call bandpass(1,shotgath(i)%gathtrace(j)%trace(:,1),tracetmp,bpparam)
          shotgath(i)%gathtrace(j)%trace(:,1)=tracetmp
       end do

    end do
    write(0,*) 'INFO: Done bandassing'
    write(0,*) 'INFO:----------------------------------------------'

    deallocate(srctmp,tracetmp)
  end subroutine BandpassSouTraces

  subroutine bandpass(ntr,    in,ou,bpparam)                
    integer ::        ntr 
    type(BandPassParam) ::      bpparam

    real, dimension(bpparam%nt,ntr):: in,ou
    real ::                         fhi_tmp,flo_tmp
    integer :: nplo_i, nphi_i
    integer :: phase_i
    integer :: n1pad, npad
    integer :: i,j
    real    :: pi=3.1415926535

    real, dimension(:), allocatable :: data,newdata,tempdata

    flo_tmp=bpparam%flo*bpparam%dt
    fhi_tmp=bpparam%fhi*bpparam%dt
    npad=10

    n1pad=bpparam%nt+2*npad
    allocate(data(n1pad),newdata(n1pad),tempdata(n1pad))

    if (flo_tmp.lt.0.0001.and.fhi_tmp.gt.0.4999) call erexit("must secify flo or fhi")

    do i=1,ntr

       data=0.
       newdata=0.
       tempdata=0.

       ! nphi and nplo are changed in the C subroutines and need
       ! to be reset for each trace
       ! --------------------------
       nphi_i=bpparam%nphi
       nplo_i=bpparam%nplo
       phase_i=bpparam%phase

!       data(1:npad)=in(1,i)
       data(npad+1:bpparam%nt+npad)=in(:,i)
!       data(n1+npad+1:)=in(n1,i)

       do j=1,npad
          data(j)=in(1,i)*sin((j-1)*pi/(2*(npad-1)))**2
          data(n1pad-npad+j)=in(bpparam%nt,i)*sin((npad-j)*pi/(2*(npad-1)))**2
       end do
 
       if (flo_tmp.gt.0.000001) then
          call lowcut(flo_tmp,nplo_i,phase_i,n1pad,data,newdata,tempdata)         
       end if
      
       if (fhi_tmp.lt.4.999999)  then
          call highcut(fhi_tmp,nphi_i,phase_i,n1pad,data,newdata,tempdata)         
       end if

       ou(:,i)=data(npad+1:bpparam%nt)

       call free_array()
    end do
    
    deallocate(data,newdata,tempdata)

  end subroutine bandpass

end module Bandpass_mod
