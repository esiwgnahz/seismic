! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module Inversion_mod

  use sep
  use Inversion_types
  use DataSpace_types
  use ModelSpace_types
  
  implicit none
  
contains
  
  !
  ! BFGS routine based on Nocedal's BFGS solver with Limited memory
  ! LINK: http://www.ece.northwestern.edu/~nocedal/lbfgs.html
  !
  ! Authors
  ! Jorge Nocedal
  !
  ! References
  ! * J. Nocedal. Updating Quasi-Newton Matrices with Limited Storage (1980), 
  !               Mathematics of Computation 35, pp. 773-782.
  ! * D.C. Liu and J. Nocedal. On the Limited Memory Method for Large Scale 
  !                            Optimization (1989), Mathematical Programming B, 45, 3, pp. 503-528.
  !
  subroutine l_bfgs(invparam,fctgdt,modin)
    
    interface
       function fctgdt(g,f,res) result (stat)
         use DataSpace_types 
         integer                ::  stat
         real, dimension(:)     ::  g
         double precision       ::  f
         type(TraceSpace), dimension(:) :: res
       end function fctgdt
    end interface

    type(InversionParam)                                         :: invparam
    type(ModelSpace)                                             :: modin 
    integer                                                      :: count

    real, dimension(:),     allocatable :: g


    double precision, dimension(:), allocatable :: gd
    double precision, dimension(:), allocatable :: xd
    double precision, dimension(:), allocatable :: xsave
    double precision, dimension(:), allocatable :: xdsave
    double precision, dimension(:), allocatable :: wd
    double precision, dimension(:), allocatable :: diagd
 
    type(TraceSpace), dimension(:), allocatable :: resigath
   
    real             :: f,g0
    double precision :: fd,f0

    integer,dimension(2) :: iprint
    double precision     :: EPS,XTOL

    logical              :: cont,found
    integer              :: i,iter,info
    integer              :: stat2,iflag,myinfo

    DOUBLE PRECISION GTOL,STPMIN,STPMAX
    DOUBLE PRECISION, dimension(invparam%nparam) :: XMIN,XMAX
    INTEGER MP,LP,NDIM,NWORK,MSAVE

    external :: LB2
    COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX

    NDIM=size(modin%vel)*invparam%nparam

    MSAVE=5
    NWORK=NDIM*(2*MSAVE+1)+2*MSAVE

    iprint(1)=1
    iprint(2)=0
    XTOL=epsilon(XTOL)
    EPS=1.0D-20

    allocate(wd(NWORK))

    allocate(diagd(NDIM))
    allocate(gd(NDIM))
    allocate( g(NDIM))
    allocate(xd(NDIM))
    allocate(xsave(NDIM))
    allocate(xdsave(NDIM))
    XMIN=dble(invparam%parmin)
    XMAX=dble(invparam%parmax)
!
    allocate(resigath(invparam%ntotaltraces))
    do i=1,invparam%ntotaltraces
       allocate(resigath(i)%trace(invparam%n1,1))
       resigath(i)%trace=0.
    end do
!
    g=0.
    gd=0.
    xd=0.
    xsave=0.
    iter=0
    iflag=0
    myinfo=0
    info=0
    count=0
    cont=.True.
    found=.False.
    ! We keep the starting guess in memory
    call mod_to_x(.true.,modin,xsave,invparam)    
    xdsave=xsave

    do while (cont.and.(.not.found))

      stat2=fctgdt(g=g,f=fd,res=resigath)

       if (invparam%eval.eq.0) then
          g0=maxval(xsave)/(100*maxval(abs(g)))
          f0=dble(maxval(xsave)**2)/abs(fd)
          if(exist_file('function')) call srite('function',sngl(f0*fd),4)
       end if

       g=g0*g 
       fd=f0*fd

       call mod_to_x(.true.,modin,xd,invparam)

       gd=dble(g)
       xdsave=xd

       call LBFGS(NDIM,MSAVE,xd,fd,-gd,&
       &          .False.,diagd,iprint,EPS,&
       &          XTOL,wd,iflag,myinfo,XMIN,XMAX,invparam%nparam,invparam%modmask,invparam%const_type,invparam%freeze_soft,invparam%modinit)

       call mod_to_x(.false.,modin,xd,invparam)
       
       info=myinfo

       if (myinfo.eq.1) then
          if (invparam%const_type.eq.2) call FREEZE_X(XMIN,XMAX,invparam%nparam,invparam%modmask,invparam%freeze_soft,XDSAVE,NDIM,invparam%modinit)
          xsave=xdsave
          myinfo=0
          invparam%iter=invparam%iter+1
          count=count+1
          if(count.ge.invparam%niter) then
             found=.true.
             write(0,*) 'INFO: iter reached max_iter - Iterations will stop'
          end if
          if(exist_file('model'))    call srite('model',sngl(xdsave),4*NDIM)
          if(exist_file('function')) call srite('function',sngl(fd),4)
          if(exist_file('gradient')) call srite('gradient',-g,4*NDIM)
          if (exist_file('residual')) then
             do i=1,invparam%ntotaltraces
                call srite('residual',resigath(i)%trace(:,1),4*invparam%n1)
             end do           
          end if
       end if

       if(iflag.le.0) then
          cont=.false.
          write(0,*) 'INFO: LINE SEARCH STOPPED IFLAG = ',iflag
          if (iflag.eq.-1) then
             if (info.eq.0) write(0,*) 'INFO:       IMPROPER INPUT PARAMETERS'
             if (info.eq.2) write(0,*) 'INFO:       INTERVAL OF UNCERTAINTY ISSUE'
             if (info.eq.3) write(0,*) 'INFO:       MORE THAN 20 FCT EVALS'
             if (info.eq.4) write(0,*) 'INFO:       STEP TOO SMALL'
             if (info.eq.4) write(0,*) 'INFO:       STEP TOO LARGE'
             if (info.eq.6) write(0,*) 'INFO:       ROUNDING ERRORS'
          end if
       end if
       invparam%eval=invparam%eval+1
       if(invparam%eval.ge.invparam%neval) then
          write(0,*) 'INFO: iter reached neval - Iterations will stop'
          cont=.false.
       end if
    end do
!    ! We didn't really find a new model vector, the line-search
    ! was still going. We take the previous best solution.
    if (.not.cont) then
       call mod_to_x(.false.,modin,xsave,invparam)
    end if

    if (allocated(xd))     deallocate(xd)
    if (allocated(xsave))  deallocate(xsave)
    if (allocated(xdsave)) deallocate(xdsave)
    if (allocated(g))      deallocate(g)
    if (allocated(wd))     deallocate(wd)
    if (allocated(gd))     deallocate(gd)
    if (allocated(diagd))  deallocate(diagd)

    do i=1,invparam%ntotaltraces
       call deallocateTraceSpace(resigath(i))
    end do
    deallocate(resigath)

  end subroutine l_bfgs
  
  subroutine copy_array(adj,mod,array1,array3)
    type(ModelSpace) ::     mod
    logical          :: adj
    double precision, dimension(:)::        array1
    real, dimension(:,:,:) ::          array3
    integer          :: i,j,k,index

    if (adj) then 
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=1,mod%ny
          do j=1,mod%nx
             do i=1,mod%nz
                index=i+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx
                array1(index)=dble(array3(i,j,k))
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    else 
       !$OMP PARALLEL DO PRIVATE(k,j,i)
       do k=1,mod%ny
          do j=1,mod%nx
             do i=1,mod%nz
                index=i+(j-1)*mod%nz+(k-1)*mod%nz*mod%nx
                array3(i,j,k)=sngl(array1(index))
             end do
          end do
       end do
       !$OMP END PARALLEL DO
    end if

  end subroutine copy_array

  subroutine mod_to_x(adj,mod,x,invparam)
    type(InversionParam)     :: invparam
    type(ModelSpace)         :: mod
    double precision, dimension(:):: x
    logical          :: adj

    call copy_array(adj,mod,x(1:size(mod%vel)),mod%vel)
       if (invparam%nparam.eq.2) then
          if (invparam%vprho_param.eq.0) then
             call copy_array(adj,mod,x(1+size(mod%vel):),mod%rho)
          else
             call copy_array(adj,mod,x(1+size(mod%vel):),mod%imp)
          end if
       end if

     end subroutine mod_to_x
    
end module Inversion_mod
