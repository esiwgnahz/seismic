! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program QcompZO_Iterative_inversion

  use sep
  use Futterman
  use futterman_fwd_adj
  use solver_smp_mod2
  use cgstep_mod2

  use identity_mod
  use hycdstep_mod
  use hycdsolver_reg_mod
  use omp_lib
  use dottest

  implicit none

  real, dimension(:,:), pointer :: futs,futsadj,futunit,futsinv
  real, dimension(:), allocatable :: model,data

  real, dimension(:,:),allocatable :: mmovie,rmovie

  integer, dimension(:), allocatable :: n
  real,    dimension(:), allocatable :: d,o

  integer :: ndim,i,j,k
  character(len=3) :: task
  character(len=1024) :: label
  real :: Qinv, dt_shift
  integer :: t2,t1,clock_rate,clock_max
  integer :: ntmax,stat,num_threads,n2dim,niter
  
  real :: threshm,threshd,eps

  real, dimension (2) :: dot1, dot2 

  call sep_init(SOURCE)
  ndim=sep_dimension()

  call from_param('num_threads',num_threads,omp_get_num_procs())
  call omp_set_num_threads(num_threads)
  write(0,*) num_threads

  allocate(n(ndim),o(ndim),d(ndim))

  do i=1,ndim
     call sep_get_data_axis_par('in',i,n(i),o(i),d(i),label)
  end do

  if (ndim.eq.1) then
     n2dim=1
  else
     n2dim=product(n(2:))
  end if

  call from_param('dt_shift',dt_shift,0.15)
  call futterman_init(n(1),d(1),dt_shift,0.375)
  call from_param('Qinv',Qinv,40.)

  allocate(futs(n(1),n(1)))
  allocate(futsadj(n(1),n(1)))
  allocate(futunit(n(1),n(1)))
  allocate(futsinv(n(1),n(1)))
  write(0,*) 'INFO: Computing Futterman wavelets with Qinv=',Qinv,'for',n(1),'time samples'

  write(0,*) 'INFO: Build futterman matrix'
  call makefuts(Qinv,futs,1.,1.,0.)
  call makefuts(Qinv,futsadj,1.,-1.,0.)
  call makefuts(Qinv,futunit,0.,-1.,0.)
  call makefuts(Qinv,futsinv,-1.,-1.,40.)

  call srite('futs',futs,4*n(1)*n(1))
  call sep_put_data_axis_par('futs',1,n(1),0.,d(1),label)
  call sep_put_data_axis_par('futs',2,n(1),0.,d(1),label)

  write(0,*) 'INFO: Applying Q compensation with Qinv=',Qinv,'using Futterman wavelets'
  allocate(data(n(1)*n2dim))
  allocate(model(n(1)*n2dim))

  call from_param('task',task,'inv')

  if (task.eq.'fwd') then
     call sreed('in',model,4*product(n))
  else
     call sreed('in',data,4*product(n))
  end if
  
  if (task.eq.'phs') then    
     call futterman_fwd_adj_init(futs,futunit,n2dim,n(1))
  else if (task.eq.'psd') then  
     call futterman_fwd_adj_init(futs,futsinv,n2dim,n(1))
  else
     futsadj=transpose(futs)
     call futterman_fwd_adj_init(futs,futsadj,n2dim,n(1))
  end if

  !call dot_test(futterman_fwd_adj_lop,n(1)*n2dim,n(1)*n2dim,dot1,dot2)

  select case(task)
  case('inv')

     call from_param('niter',niter,100)
     call from_param('threshd',threshd,maxval(data))
     call from_param('threshm',threshm,maxval(data)/100)
     call from_param('eps',eps,1.)
     call identity_init(threshd,threshm)
     allocate(mmovie(n(1)*n2dim,niter))
     allocate(rmovie(n(1)*n2dim,niter))

     write(0,*) 'INFO:'
     write(0,*) 'INFO: --------------------'
     write(0,*) 'INFO: Inversion parameters'
     write(0,*) 'INFO: --------------------'
     write(0,*) 'INFO:'
     write(0,*) 'INFO:  niter   =', niter
     write(0,*) 'INFO:  threshd =', threshd
     write(0,*) 'INFO:  threshm =', threshm
     write(0,*) 'INFO:  eps     =', eps
     write(0,*) 'INFO: --------------------'
     write(0,*) 'INFO:'

     call hycdsolver_reg(model,data,futterman_fwd_adj_lop,identity_lop,identityd_lop,&
     &  identitym_lop,hycdstep2,product(n),niter,eps,verb=.true.,mmov=mmovie,rmov=rmovie)

     do i=1,niter
        call srite('rmovie',rmovie(:,i),4*n(1)*n2dim) 
        call srite('mmovie',mmovie(:,i),4*n(1)*n2dim)  
     end do
     
  case('fwd')

     write(0,*) 'INFO:'
     write(0,*) 'INFO: --------------------------------------'
     write(0,*) 'INFO: Forward Futterman operator: Applying Q'
     write(0,*) 'INFO: --------------------------------------'
     write(0,*) 'INFO:'
     
     stat=futterman_fwd_adj_lop(.false.,.false.,model,data)

  case('adj')
     
     write(0,*) 'INFO:'
     write(0,*) 'INFO: ---------------------------'
     write(0,*) 'INFO: Adjoint Futterman operator '
     write(0,*) 'INFO: ---------------------------'
     write(0,*) 'INFO:'
     stat=futterman_fwd_adj_lop(.true.,.false.,model,data)

  case('phs')

     write(0,*) 'INFO:'
     write(0,*) 'INFO: --------------------------------------'
     write(0,*) 'INFO: Phase only inverse Futterman operator '
     write(0,*) 'INFO: --------------------------------------'
     write(0,*) 'INFO:'    
     stat=futterman_fwd_adj_lop(.true.,.false.,model,data)

  case('psd')

     write(0,*) 'INFO:'
     write(0,*) 'INFO: --------------------------------------'
     write(0,*) 'INFO: Pseudo inverse Futterman operator '
     write(0,*) 'INFO: --------------------------------------'
     write(0,*) 'INFO:'    
     stat=futterman_fwd_adj_lop(.true.,.false.,model,data)

  end select

  if (task.eq.'fwd') then
     call srite('out',data,4*product(n))
  else
     call srite('out',model,4*product(n))
  end if
 
  write(0,*) 'INFO: Done with Q compensation with iterative inversion'
  deallocate(futs,data,model,futsadj,futunit,futsinv)

  if (task.eq.'inv') then
     do i=1,ndim
        call sep_put_data_axis_par('rmovie',i,n(i),o(i),d(i),label)
        call sep_put_data_axis_par('mmovie',i,n(i),o(i),d(i),label)
     end do
     call sep_put_data_axis_par('rmovie',ndim+1,niter,0,1,'Iteration')
     call sep_put_data_axis_par('mmovie',ndim+1,niter,0,1,'Iteration')
  end if

  do i=1,ndim
     call sep_put_data_axis_par('out',i,n(i),o(i),d(i),label)
  end do
  
end program QcompZO_Iterative_inversion
