! 
! -----------------------------------------------
! Copyright (c) 2016-2017 Bellevue Geophysics LLC
! -----------------------------------------------
! 
module solver_smp_mod2
  ! 0 = W (F J m - d)
  use omp_lib
  use chain0_mod 
  use  solver_report_mod
  implicit none
  logical, parameter, private  :: T = .true., F = .false.
  contains
  subroutine solver_smp2( m,d, Fop, stepper, niter   ,             Wop&
    &,Jop,m0,err,resd,mmov,rmov,verb)
    optional :: Wop,Jop,m0,err,resd,mmov,rmov,verb
    interface  
!-------------------------- begin definitions -----------
      integer function Fop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function 
      integer function Wop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function 
      integer function Jop(adj,add,m,d)
        real::m(:),d(:)
        logical,intent(in)::adj,add
      end function 
      integer function stepper(forget,m,g,rr,gg)
        real, dimension(:)  ::        m,g,rr,gg    
        logical             :: forget
      end function
    end interface 
    real, dimension(:),    intent(in)          :: d, m0
    integer,               intent(in)          :: niter
    logical,               intent(in)          :: verb
    real, dimension(:),    intent(out)         :: m,err, resd
    real, dimension(:,:),  intent(out)         ::        rmov, mmov   
    real, dimension(size(m))                   :: g
    real, dimension(size(d)),         target   :: rr, gg
    real, dimension(size(d)+size(m)), target   :: tt
    real, dimension(:), pointer                :: rd, gd, td
    real, dimension(:), pointer                :: rm, gm, tm
    integer                                    :: iter, stat
    logical                                    :: forget

    write(0,*) 'starting solver'
    rd => rr(1:size(d))
    gd => gg(1:size(d))
    td => tt(1:size(d))
    tm => tt(1+size(d):)
    if (present( Wop)) then
      stat=Wop(F,F,-d,rd) ! begin initialization --------
    else
      rd = -d                                           !Rd = -W  d
    end if 
    if (present( m0)) then
      m=m0                                 !m  =     m0
      if (present( Wop)) then
        call chain0(Wop,Fop,F,T,m,rd,td)
      else
        stat =           Fop(F,T,m,rd   ) !Rd+=  WF m0
      end if
    else
      m=0
    end if 
    forget = T
!-------------------------- begin iterations ------------
    do iter = 1,niter  
       write(0,*) 'here solver'
      if (present(Wop)) then
        call chain0(Wop,Fop,T,F,g,rd,td)
      else
        stat =          Fop(T,F,g,rd   )   !g  = (WF)'Rd
      end if 
       write(0,*) 'here solver 1'
      if (present(Jop)) then
        tm=g
        stat  = Jop(F,F,tm, g  )
      end if
!g  =   J  g
       write(0,*) 'here solver 2'
      if (present(Wop)) then
        call chain0(Wop,Fop,F,F,g,gd,td)
      else
        stat =          Fop(F,F,g,gd   )   !Gd = (WF) g
      end if 
       write(0,*) 'here solver3'
      stat = stepper(forget, m,g, rr,gg)   
       write(0,*) 'here solver4'               
      !m+=dm; R+=dR
      if (stat .eq.1) then
        exit ! got stuck descending
      end if
      if (present( mmov)) then
        mmov(:,iter) =  m(:size(mmov,1)) ! report -----
      end if
      if (present( rmov)) then
        rmov(:,iter) = rd(:size(rmov,1))
      end if
      if (present( err )) then
        err(  iter) = dot_product(rd,rd)
      end if
      if (present( verb)) then
        if (verb) then
          call solver_report(iter,m,g,rd)
        end if
      end if
      forget=F
    end do 
    if (present( resd)) then
      resd = rd
    end if
  end subroutine 
end module 
