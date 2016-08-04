module chain_mod

  implicit none
  logical, parameter, private :: T = .true., F = .false.
!!$
!!$=head1 NAME
!!$
!!$chain  - create a chain of 2,3, or 4 operators
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call chain(oper1,oper2,adj,add,model,data,temp1)>
!!$
!!$C<call chain(oper1,oper2,oper3,adj,add,model,data,temp1,temp2)>
!!$
!!$C<call chain(oper1,oper2,oper3,oper4,adj,add,model,data,temp1,temp2,temp3)>
!!$
!!$=head1 PARAMETERS
!!$
!!$=over 4
!!$
!!$=item oper1,oper2,oper3,oper4 - func  
!!$
!!$      standard operators
!!$
!!$=item adj,add,model,data      -
!!$
!!$      standard operator interface
!!$
!!$=item temp1,temp2,temp3       - C<real(:)>  
!!$
!!$      temporary arrays
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Create a chain of multiple operators
!!$
!!$head1 COMMENTS
!!$
!!$temp1 should be the output of oper1(model)
!!$temp2 should be the output of oper2(oper1(model))
!!$temp3 should be the output of oper3(oper2(oper1(model)))
!!$
!!$=head1 SEE ALSO
!!$
!!$L<array>,L<solver>
!!$
!!$=head1 LIBRARY
!!$
!!$B<geef90>
!!$
!!$=cut
!!$

  interface chain
     module procedure chain2
     module procedure chain3
     module procedure chain4
  end interface

  interface array_chain
     module procedure array_chain2
     module procedure array_chain3
  end interface

contains

  subroutine doub (oper, mod, dat)
    interface
       function oper (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper
    end interface
    real, dimension (:)          :: mod, dat

    real, dimension (size (mod)) :: temp, tp1, tp2
    integer                      :: stat1, stat2

    stat1 = oper (T, F, temp, mod)
    stat2 = oper (F, F, temp, tp1)

    stat2 = oper (F, F, mod, temp)
    stat1 = oper (T, F, tp2, temp)

    dat = tp1 + tp2
  end subroutine doub

  subroutine asym (oper, adj, mod, dat)
    interface
       function oper (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper
    end interface
    logical, intent (in)          :: adj
    real, dimension (:)           :: mod, dat

    real, dimension (size (mod))  :: temp1, temp2
    integer                       :: stat1, stat2

    if (adj) then
       stat1 = oper (T, F, temp1, dat)
       stat2 = oper (F, F, dat, temp2)
       mod = mod + temp1 - temp2
    else
       stat2 = oper (F, F, mod, temp2)
       stat1 = oper (T, F, temp1, mod)
       dat = dat + temp2 - temp1
    end if
  end subroutine asym

  subroutine normal (oper, add, mod, dat, tmp)
    interface
       function oper (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper
    end interface
    logical, intent (in)  :: add
    real, dimension (:)   :: mod, dat, tmp

    integer               :: stat1, stat2

    stat2 = oper (T, F,   tmp, mod)
    stat1 = oper (F, add, tmp, dat)
  end subroutine normal

  subroutine anormal (oper, add, mod, dat, tmp)
    interface
       function oper (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper
    end interface
    logical, intent (in)  :: add
    real, dimension (:)   :: mod, dat, tmp

    integer               :: stat1, stat2

    stat2 = oper (F, F,   mod, tmp)
    stat1 = oper (T, add, dat, tmp)
  end subroutine anormal

  subroutine chain2 (oper1, oper2, adj, add, mod, dat, tmp1)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper2
    end interface
    logical, intent (in)  :: adj, add
    real, dimension (:)   :: mod, dat, tmp1

    integer               :: stat1, stat2

    if (adj) then
       stat1 = oper1 (T,   F, tmp1, dat)
       stat2 = oper2 (T, add, mod, tmp1)
    else
       stat2 = oper2 (F,   F, mod, tmp1)
       stat1 = oper1 (F, add, tmp1, dat)
    end if
  end subroutine chain2

  subroutine chain3 (oper1, oper2, oper3, adj, add, mod, dat, tmp1, tmp2)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper2
       function oper3 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper3
    end interface
    logical, intent (in)  :: adj, add
    real, dimension (:)   :: mod, dat, tmp1, tmp2

    integer               :: stat1, stat2, stat3

    if (adj) then
       stat1 = oper1 (T,   F, tmp2, dat)
       stat2 = oper2 (T,   F, tmp1, tmp2)
       stat3 = oper3 (T, add, mod,  tmp1)
    else
       stat3 = oper3 (F,   F, mod,  tmp1)
       stat2 = oper2 (F,   F, tmp1, tmp2)
       stat1 = oper1 (F, add, tmp2, dat)
    end if
  end subroutine chain3

  subroutine chain4 &
  (oper1, oper2, oper3, oper4, adj, add, mod, dat, tmp1, tmp2, tmp3)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper2
       function oper3 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper3
       function oper4 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper4
    end interface
    logical, intent (in)  :: adj, add
    real, dimension (:)   :: mod, dat, tmp1, tmp2, tmp3

    if (adj) then
       call chain2 (oper1, oper2, T,   F, tmp2, dat, tmp3)
       call chain2 (oper3, oper4, T, add, mod, tmp2, tmp1)
    else
       call chain2 (oper3, oper4, F,   F, mod, tmp2, tmp1)
       call chain2 (oper1, oper2, F, add, tmp2, dat, tmp3)
    end if
  end subroutine chain4

  subroutine chain_square (oper1, oper2, adj, add, mod, dat, tmp1, tmp2)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer               :: stat
         logical, intent (in)  :: adj, add
         real, dimension (:)   :: mod, dat
       end function oper2
    end interface
    logical, intent (in)  :: adj, add
    real, dimension (:)   :: mod, dat, tmp1, tmp2

    integer               :: stat1, stat2, stat3

    if (adj) then
       stat1 = oper1 (T,   F, tmp2, dat)
       stat2 = oper2 (T,   F, tmp1, tmp2)
       stat3 = oper2 (T, add, mod,  tmp1)
    else
       stat3 = oper2 (F,   F, mod,   tmp1)
       stat2 = oper2 (F,   F, tmp1, tmp2)
       stat1 = oper1 (F, add, tmp2, dat)
    end if
  end subroutine chain_square
!!$=head1 NAME
!!$
!!$array  - create an array operator
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call array(oper1,oper2,adj,add,mod,dat1,dat2)>
!!$
!!$=head1 PARAMETERS
!!$
!!$=over 4
!!$
!!$=item oper1,oper2  -  func
!!$
!!$      standard operator
!!$
!!$=item adj,add,mod -
!!$
!!$      standard operator interface
!!$
!!$=item  dat1,dat2  -  real(:)
!!$
!!$       portions of data
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Create an opertor which actually a compound operator , eg:
!!$
!!$  | dat1  |    | oper1|
!!$
!!$  | ----  |  = | ---- |   model
!!$
!!$  | dat2  |    |oper2 |
!!$
!!$
!!$=head1 SEE ALSO
!!$
!!$L<chain>,L<solver>
!!$
!!$=head1 LIBRARY
!!$
!!$B<geef90>
!!$
!!$=cut


  subroutine array (oper1, oper2, adj, add, mod, dat1, dat2)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper2
    end interface
    logical, intent (in)  :: adj, add
    real, dimension (:)   :: mod, dat1, dat2

    integer               :: stat1, stat2

    if (adj) then
       stat1 = oper1 (T, add, mod, dat1)
       stat2 = oper2 (T,   T, mod, dat2)
    else
       stat1 = oper1 (F, add, mod, dat1)
       stat2 = oper2 (F, add, mod, dat2)
    end if
  end subroutine array

  subroutine chain_array (oper1, oper2, oper3, adj, add, mod, dat1, dat2, tmp1)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper2
       function oper3 (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper3
    end interface
    logical, intent (in)  :: adj, add
    real, dimension (:)   :: mod, dat1, dat2, tmp1

    integer               :: stat

    if (adj) then
       call chain (oper1, oper2, T, add, mod, dat1, tmp1)
       stat =             oper3 (T,   T, mod, dat2)
    else
       call chain (oper1, oper2, F, add, mod, dat1, tmp1)
       stat =             oper3 (F, add, mod, dat2)
    end if
  end subroutine chain_array

  subroutine row (oper1, oper2, adj, add, mod1, mod2, dat)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer    :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper2
    end interface
    logical, intent (in)  :: adj, add
    real, dimension (:)   :: mod1, mod2, dat

    integer               :: stat1, stat2

    if (adj) then
       stat1 = oper1 (T, add, mod1, dat)
       stat2 = oper2 (T, add, mod2, dat)
    else
       stat1 = oper1 (F, add, mod1, dat)
       stat2 = oper2 (F,   T, mod2, dat)
    end if
  end subroutine row

  subroutine array_chain2 (oper1, oper2, oper3, adj, add, mod, dat1, dat2, tmp1)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer    :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper2
       function oper3 (adj, add, mod, dat) result (stat)
         integer    :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper3
    end interface
    logical, intent (in)   :: adj, add
    real, dimension (:)    :: mod, dat1, dat2, tmp1

    if (adj) then
       call chain (oper1, oper3, T, add, mod, dat1, tmp1)
       call chain (oper2, oper3, T,   T, mod, dat2, tmp1)
    else
       call chain (oper1, oper3, F, add, mod, dat1, tmp1)
       call chain (oper2, oper3, F, add, mod, dat2, tmp1)
    end if
  end subroutine array_chain2

 subroutine array_chain3 &
 (oper1, oper2, oper3, oper4, adj, add, mod, dat1, dat2, tmp1, tmp2)
    interface
       function oper1 (adj, add, mod, dat) result (stat)
         integer              :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper1
       function oper2 (adj, add, mod, dat) result (stat)
         integer    :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper2
       function oper3 (adj, add, mod, dat) result (stat)
         integer    :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper3
       function oper4 (adj, add, mod, dat) result (stat)
         integer    :: stat
         logical, intent (in) :: adj, add
         real, dimension (:)  :: mod, dat
       end function oper4
    end interface
    logical, intent (in)   :: adj, add
    real, dimension (:)    :: mod, dat1, dat2, tmp1, tmp2

    if (adj) then
       call chain (oper1, oper2, oper4, T, add, mod, dat1, tmp1, tmp2)
       call chain (oper3, oper4,        T,   T, mod, dat2, tmp1)
    else
       call chain (oper1, oper2, oper4, F, add, mod, dat1, tmp1, tmp2)
       call chain (oper3, oper4,        F, add, mod, dat2, tmp1)
    end if
  end subroutine array_chain3

end module chain_mod




