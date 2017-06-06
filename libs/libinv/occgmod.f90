! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module occgmod_mod
  
  implicit none

contains

  double precision function dot_product(dat1,dat2)
    real, dimension(:) :: dat1,dat2

    dot_product=sum(dprod(dat1,dat2))

  end function dot_product

  subroutine compute_steps(iter,gdg,sds,gds,gdr,sdr,alfa,beta)
    integer          ::    iter
    double precision ::         gdg,sds,gds,gdr,sdr,alfa,beta
    double precision :: determ

    if (iter.eq.1) then
       beta = 0.d0 
       alfa = gdr/gdg
    else
       determ = gdg * sds * max (1.d0 - (gds/gdg)*(gds/sds), 1.d-12)
       alfa = ( sds * gdr - gds * sdr ) / determ
       beta = (-gds * gdr + gdg * sdr ) / determ
    end if
    
  end subroutine compute_steps

  subroutine update_array_step(alfa,beta,dat1,dat2)
    double precision  ::       alfa,beta
    real, dimension(:)::                 dat1,dat2

    dat1=alfa*dat2+beta*dat1

  end subroutine update_array_step

end module occgmod_mod
