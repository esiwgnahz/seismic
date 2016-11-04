program Ddot_gg_rr_ss_OC

  use sep
  use occgmod_mod
  use ReadData_mod
  use DataSpace_types_flt

  implicit none

  type(cube)       :: gg,rr,ss
  double precision :: gdg,sds,gdr,gds,sdr,rdr
  logical          :: add
  integer :: ndim,iter
  call sep_init(SOURCE)

  ndim=sep_dimension('gg')

  call ReadData_dim('gg',gg,ndim)
  call ReadData_dim('rr',rr,ndim)
  call ReadData_dim('ss',ss,ndim)
  if (.not.ns_are_consistent(ndim,gg,rr)) call erexit('ERROR: dimensions gg and rr are not the same, exist now')
  if (.not.ns_are_consistent(ndim,gg,ss)) call erexit('ERROR: dimensions gg and ss are not the same, exist now')

  call ReadData_cube('gg',gg)
  call auxclose('gg')
  call ReadData_cube('rr',rr)
  call auxclose('rr')
  call ReadData_cube('ss',ss)
  call auxclose('ss')

  call from_param('add',add,.false.)
  call from_param('iter',iter)

  write(0,*) 'INFO: add=',add
  write(0,*) 'INFO: iter=',iter

  if(.not.exist_file('dot_product')) call erexit('ERROR: need dot_product file to save value')

  if (add) then
     call sreed('dot_product',gdg,8)
     call sreed('dot_product',sds,8)
     call sreed('dot_product',gdr,8)
     call sreed('dot_product',gds,8)
     call sreed('dot_product',sdr,8)
     call sreed('dot_product',rdr,8)
     call auxclose('dot_product')
  else
     gdg=0.d0
     sds=0.d0
     gdr=0.d0
     gds=0.d0
     sdr=0.d0
     rdr=0.d0
     call to_history('n1',6,'dot_product')
     call to_history('esize',8,'dot_product')
  end if

  gdg=gdg+dot_product(gg%dat,gg%dat)
  sds=sds+dot_product(ss%dat,ss%dat)
  gdr=gdr-dot_product(gg%dat,rr%dat)
  gds=gds+dot_product(gg%dat,ss%dat)
  sdr=sdr-dot_product(ss%dat,rr%dat)
  rdr=rdr+dot_product(rr%dat,rr%dat)

  call srite('dot_product',gdg,8)
  call srite('dot_product',sds,8)
  call srite('dot_product',gdr,8)
  call srite('dot_product',gds,8)
  call srite('dot_product',sdr,8)
  call srite('dot_product',rdr,8)
  
  call to_history('n1',6,'dot_product')
  call to_history('esize',8,'dot_product')

end program Ddot_gg_rr_ss_OC
