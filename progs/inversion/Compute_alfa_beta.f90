! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
program Compute_alfa_beta

  use sep
  use occgmod_mod

  double precision :: gdg,sds,gdr,gds,sdr
  double precision :: alfa,beta
  integer          :: iter

  call sep_init(SOURCE)
  call from_param('iter',iter)

  write(0,*) 'INFO: iter=',iter
  if(.not.exist_file('dot_product')) call erexit('ERROR: need dot_product file to save value')

  call sreed('dot_product',gdg,8)
  call sreed('dot_product',sds,8)
  call sreed('dot_product',gdr,8)
  call sreed('dot_product',gds,8)
  call sreed('dot_product',sdr,8)

  call compute_steps(iter,gdg,sds,gds,gdr,sdr,alfa,beta)
  
  call srite('alfabeta',alfa,8)
  call srite('alfabeta',beta,8)
  call to_history('n1',2,'alfabeta')
  call to_history('esize',8,'alfabeta')

end program Compute_alfa_beta
