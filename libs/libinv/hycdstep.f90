module hycdstep_mod
  use hpf_mod
  implicit none
  contains
  integer function hycdstep2(sd,alfa_keep,beta_keep,ss,rd,gg,rm,eps)
    real, dimension (:):: ss,rd,gg
    real, dimension (:), optional :: rm
    double precision :: alfa_keep, beta_keep
    real, optional :: eps
    logical :: sd ! steepest descent or not
    real, dimension (:), allocatable :: hp,hs,tmpd,tmpm
    double precision :: alfa, beta, h, hm, hsave
    double precision :: sds, gdg, gds, determ, gdr, sdr
    integer :: nsteps,k,ntries,l,nd,nm
    real :: cst
    logical :: cont, icont, reg

    reg=present(rm)
    nsteps=5
    ntries=5
    nm=0
    nd=size(rd)
    if (reg) then
      nm=size(rm)
    end if
    allocate (hp(nd+nm))
    allocate (hs(nd+nm))
    allocate (tmpd(nd))
    if (reg) then
      allocate (tmpm(nm))
    end if
    k=0
    alfa_keep = 0.d0 ! gradient direction scaling
    beta_keep = 0.d0 ! previous direction scaling
    alfa = 0.d0
    beta = 0.d0
    cst = 0.001
    hm = hpf(rd)
    if (reg) then
      hm= hm+eps*hpf(rm)
    end if
    h = 0.
    cont=.true.

    if (sd) then
! Line search
      do while (cont)
        hp(1:nd) = hpfp(rd)
        if (reg) then
          hp(nd+1:)= eps*hpfp(rm)
        end if
        hs(1:nd) = hpfs(rd)
        if (reg) then
          hs(nd+1:)= eps*hpfs(rm)
        end if
        if ( dot_product(hs*gg, gg) .eq. 0 ) then
          call erexit('hycdstep: grad vanishes identically')
        end if
        alfa = - sum( dprod( gg, hp)) / sum( dprod( hs*gg, gg))
        tmpd=rd+alfa*gg(1:nd)
        if (reg) then
          tmpm=rm+alfa*gg(nd+1:)
        end if
        h = hpf(tmpd)
        if (reg) then
          h= h+eps*hpf(tmpm)
        end if
        if (abs(alfa).lt.epsilon(abs(alfa))) then
          exit
        end if
        hsave = h
        icont=.false.
        if (h.gt.hm) then
          icont=.true.
        end if
        do while (icont)
          alfa=alfa/2
          tmpd=rd+alfa*gg(1:nd)
          if (reg) then
            tmpm=rm+alfa*gg(nd+1:)
          end if
          h = hpf(tmpd)
          if (reg) then
            h= h+eps*hpf(tmpm)
          end if
          if (hsave.eq.h) then
            icont=.false.
          end if
          if (h.le.hm) then
            icont=.false.
          end if
          hsave = h
        end do
!write(0,*) 'h',h,'hm',hm
        rd = rd + alfa*gg(1:nd)
        if (reg) then
          rm = rm + alfa*gg(nd+1:)
        end if
        alfa_keep=alfa_keep+alfa
        k=k+1
        if (h.gt.hm) then
          cont=.false.
        end if
        if (k.ge.nsteps) then
          cont=.false.
        end if
        hm = h
      end do
    else
      do while (cont)
        hp(1:nd) = hpfp(rd)
        if (reg) then
          hp(nd+1:)= eps*hpfp(rm)
        end if
        hs(1:nd) = hpfs(rd)
        if (reg) then
          hs(nd+1:)= eps*hpfs(rm)
        end if
        gdg = sum( dprod( hs*gg, gg))
        ! --- ---
        sds = sum( dprod( hs*ss, ss))
        ! \ H_i''|(gi) (gi si)|| alfa |= -\ H_i'|gi|
        gds = sum( dprod( hs*gg, ss))
        ! /__ |(si) || beta | /__ |si|
        gdr = - sum( dprod( gg, hp))
        ! i i
        sdr = - sum( dprod( ss, hp))
        if ( gdg.eq.0. .or. sds.eq.0.) then
          hycdstep2 = 1
          return
        end if
        determ = gdg * sds * max (1.d0 - (gds/gdg)*(gds/sds), 1.d-12)
        alfa = ( sds * gdr - gds * sdr ) / determ
        beta = (-gds * gdr + gdg * sdr ) / determ
        tmpd=rd+alfa*gg(1:nd)+beta*ss(1:nd)
        tmpm=rm+alfa*gg(nd+1:)+beta*ss(nd+1:)
        h = hpf(tmpd)
        if (reg) then
          h= h+eps*hpf(tmpm)
        end if
        hsave = h
        icont=.false.
        if (h.gt.hm) then
          icont=.true.
        end if
        do while (icont)
          alfa=alfa/2
          beta=beta/2
          tmpd=rd+alfa*gg(1:nd)+beta*ss(1:nd)
          if (reg) then
            tmpm=rm+alfa*gg(nd+1:)+beta*ss(nd+1:)
          end if
          h = hpf(tmpd)
          if (reg) then
            h= h+eps*hpf(tmpm)
          end if
          if (hsave.eq.h) then
            icont=.false.
          end if
          if (h.le.hm) then
            icont=.false.
          end if
          hsave = h
        end do
        rd = rd+alfa*gg(1:nd)+beta*ss(1:nd)
        if (reg) then
          rm = rm+alfa*gg(nd+1:)+beta*ss(nd+1:)
        end if
        alfa_keep=alfa_keep+alfa
        beta_keep=beta_keep+beta
        k=k+1
        if (h.gt.hm) then
          cont=.false.
        end if
        if (k.gt.nsteps) then
          cont=.false.
        end if
        hm = h
      end do
    end if
    hycdstep2=0
    deallocate(hp,hs,tmpd)
    if (reg) then
      deallocate(tmpm)
    end if
  end function
end module
