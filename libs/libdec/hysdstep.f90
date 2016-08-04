module hycdstep_mod 

  implicit none
  
  integer, private :: n,nreg
  
contains

  subroutine hycdstep_init(n_init)
    integer :: n_init
    n=n_init
  end subroutine hycdstep_init

  subroutine hycdstep_nreg_init(nreg_init)
    integer :: nreg_init
    nreg=nreg_init
  end subroutine hycdstep_nreg_init

  integer function hycdstep(alfa_keep, rr, gg, newton_steps) result(stat)
    real, dimension (n) ::             rr, gg
    double precision    ::  alfa_keep
    integer, intent(out) :: newton_steps
    
    real, dimension (:), allocatable :: hp,hs
    double precision     :: alfa, h, hm, hsave
    integer              :: nsteps,k,l
    real                 :: cst
    logical              :: cont, icont, giveup
    
    newton_steps=0
    nsteps=10
    allocate (hp (size (rr)))
    allocate (hs (size (rr)))
    
    k=0
    alfa_keep = 0.d0    ! gradient direction scaling
    alfa = 0.d0
    
    cst  = 0.001
    hm  = sum(sqrt(1.d0+rr**2)-1.d0)
    h   = 0.
    cont=.true.
    giveup=.false.

    do while(cont.and.(.not.giveup))
       hp = rr/sqrt(1+rr**2)    
       hs = 1./(1.+rr**2)**1.5
       if( dot_product(hs*gg, gg) == 0 ) call erexit('hycdstep: grad vanishes identically')
       alfa = - sum( dprod( gg, hp)) / sum( dprod( hs*gg, gg))
       
       h  = sum(sqrt(1.d0+(rr+alfa*gg)**2)-1.d0)
       if (abs(alfa).lt.epsilon(abs(alfa))) exit         
       
       hsave = h
       icont=.false.
       if (h.gt.hm) icont=.true.
       do while (icont.and.(.not.giveup)) 
          alfa=alfa/2
          h  = sum(sqrt(1.d0+(rr+alfa*gg)**2)-1.d0)
          if (hsave.eq.h) giveup=.true.
          if (h.le.hm) icont=.false.
          hsave = h
       end do
       
       newton_steps=newton_steps+1
       rr = rr + alfa*gg
       alfa_keep=alfa_keep+alfa 
       
       k=k+1
       if (h.gt.hm) cont=.false.
       if (k.ge.nsteps) giveup=.true.
       
       hm = h
       
    end do
    
    if (giveup) then
       stat=2
    else
       stat=1
    end if

    deallocate(hs,hp)

  end function hycdstep
  
  integer function hycdstep_reg(alfa_keep, rr, gg, regt, dut, eps, newton_steps) result(stat)
    real, dimension (n) ::                 rr, gg
    real, dimension (nreg) ::                      regt, dut
    real                ::                                    eps
    double precision    ::      alfa_keep
    integer, intent(out)::                                         newton_steps
    
    real, dimension (:), allocatable :: hp,hs,tmp
    double precision     :: alfa, h, hm, hreg, hdat, hsave
    integer              :: nsteps,k,l
    real                 :: cst
    logical              :: cont, icont, giveup
    
    newton_steps=0
    nsteps=10
    allocate (hp (size (rr)))
    allocate (hs (size (rr)))
    allocate (tmp(size (regt)))
    
    k=0
    alfa_keep = 0.d0    ! gradient direction scaling
    alfa = 0.d0
    
    tmp = regt
    cst = 0.00
    hreg= eps*sum(regt**2)/2.d0
    hdat= sum(sqrt(1.d0+rr**2)-1.d0)
    hm  = hdat+hreg
    h   = 0.
    cont=.true.
    giveup=.false.

    do while(cont.and.(.not.giveup))
       hp = rr/sqrt(1+rr**2)    
       hs = 1./(1.+rr**2)**1.5
       if( dot_product(hs*gg, gg) == 0 ) call erexit('hycdstep: grad vanishes identically')
       alfa = - (sum( dprod( hp,    gg))+eps*sum(dprod(tmp,dut))) / &
       &        (sum( dprod( hs*gg, gg))+eps*sum(dprod(dut,dut)))
       
       hreg= eps*sum((tmp+alfa*dut)**2)/2.d0
       hdat= sum(sqrt(1.d0+(rr+alfa*gg)**2)-1.d0)
       h   = hdat+hreg
       if (abs(alfa).lt.epsilon(abs(alfa))) exit         
       
       hsave = h
       icont=.false.
       if (h.gt.hm) icont=.true.
       do while (icont.and.(.not.giveup)) 
          alfa=alfa/2
          hreg= eps*sum((tmp+alfa*dut)**2)/2.d0
          hdat= sum(sqrt(1.d0+(rr+alfa*gg)**2)-1.d0)
          h   = hdat+hreg
          if (hsave.eq.h) giveup=.true.
          if (h.le.hm) icont=.false.
          hsave = h
       end do
       
       newton_steps=newton_steps+1
       rr = rr + alfa*gg
       tmp=tmp + alfa*dut
       alfa_keep=alfa_keep+alfa 
       
       k=k+1
       if (h.gt.hm) cont=.false.
       if (k.ge.nsteps) giveup=.true.
       
       hm = h
       
    end do
    
    if (giveup) then
       stat=2
    else
       stat=1
    end if

    deallocate(hs,hp,tmp)

  end function hycdstep_reg
  
end module hycdstep_mod
