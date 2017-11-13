program blocky_model

  use sep
  use helix
  use print
  use Norm_mod
  use helicon_mod
  use createhelixmod

  real, dimension(:), allocatable   :: m, g, tmp
  real                              :: f
  integer                           :: ndim,npanels,nsamples,nthreads
  integer, dimension(:), allocatable:: n0,x,y,z,center,gap,npef
  integer, dimension(:), allocatable:: n
  real,    dimension(:), allocatable:: o
  real,    dimension(:), allocatable:: d
  character(len=1024)               :: label
  type(filter)                      :: xx,zz
  integer                           :: stat,norm
  real                              :: thresh

  call sep_init(SOURCE)
  ndim=sep_dimension()
  
  if (ndim.eq.2) ndim=3

  allocate(n(ndim),o(ndim),d(ndim))
  n=1
  
  do i=1,ndim
     call sep_get_data_axis_par("in",i,n(i),o(i),d(i),label)
  end do

  npanels=product(n(3:))
  nsamples=n(2)*n(1)
  allocate(m(nsamples))
  allocate(g(nsamples))
  allocate(tmp(nsamples))

  allocate(x(2),y(2),z(2),center(2),gap(2),n0(2),npef(2))
  center=1; gap=0; npef=(/n1,n2/); n0=npef
  z=(/2,1/)
  x=(/1,2/) 
  
  zz=createhelix(npef,center,gap,z)
  xx=createhelix(npef,center,gap,x)            
  zz%flt=-1
  xx%flt=-1

  write(0,*) 'INFO: Filter z:'
  call printn(n0,center,z,zz)
  write(0,*) 'INFO: Filter x:' 
  call printn(n0,center,x,xx)

  call from_param('num_threads',nthreads,4)
  write(0,*) 'INFO: nthreads=',nthreads
  call omp_set_num_threads(nthreads)
  call from_param('norm',norm,2)
  if (norm.eq.2) then
     write(0,*) 'INFO: using L2 norm'
  else if (norm.eq.1) then
     write(0,*) 'INFO: using L1 norm'
  else if (norm.eq.12) then
     write(0,*) 'INFO: using Huber norm'
     call from_param('thresh',thresh)
     write(0,*) 'INFO: with threshold',thresh
  else if (norm.eq.3) then
     write(0,*) 'INFO: using Cauchy norm'
     call from_param('thresh',thresh)
     write(0,*) 'INFO: with hyperparameter',thresh
  else if (norm.eq.4) then
     write(0,*) 'INFO: using Hyperbolic functional'
     call from_param('thresh',thresh)
     write(0,*) 'INFO: with hyperparameter',thresh
  end if
  call helicon_mod_init(zz)

  f=0.
  do i=1,npanels
     if (mod(i,5)) write(0,*) 'INFO: Processing panel',i,'/',npanels
     m=0.;g=0.;tmp=0.
     if (mod(i,5)) write(0,*) 'INFO:    read panel ',i
     call sreed('in',m,4*nsamples)
     if (mod(i,5)) write(0,*) 'INFO:    convolve panel ',i,' with Dz'
     stat= helicon_mod_lop(.false.,.false.,m,tmp)
     if (mod(i,5)) write(0,*) 'INFO:    compute fct for panel ',i
     f=f+fct_compute(norm,tmp,nsamples,thresh)
     if (mod(i,5)) write(0,*) 'INFO:    compute gdt for panel ',i
     stat=gdt_compute(norm,tmp,nsamples,thresh)
     if (mod(i,5)) write(0,*) 'INFO:    compute adjt for panel ',i
     stat= helicon_mod_lop(.true.,.false.,g,tmp)
     if (mod(i,5)) write(0,*) 'INFO:    write panel ',i
     call srite('out',g,4*nsamples)
  end do
  write(0,*) 'INFO: Now write objective function'
  call srite('fct',f,4)
  call to_history('n1',1,'fct')
  call to_history('n2',1,'fct')
  call to_history('n3',1,'fct')

  deallocate(m,g,tmp)
  write(0,*) 'INFO: Done with fct/gdt computation'
end program blocky_model
