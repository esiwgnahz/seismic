!!$=head1 NAME
!!$
!!$adj_mod - Simple adjnull function
!!$
!!$=head1 SYNOPSIS
!!$
!!$C<call adjnull(adj,add,x,y)>
!!$
!!$=head1 PARAMETERS
!!$
!!$=over 4
!!$
!!$=item adj,add,model,data -
!!$
!!$      Standard operator parameters
!!$
!!$=back
!!$
!!$=head1 DESCRIPTION
!!$
!!$Erases output if add=.false.
!!$
!!$=head1 LIBRARY
!!$
!!$B<geef90>
!!$
!!$=cut
!!$
module adj_mod
contains
  subroutine adjnull(       adj, add, x, y)
    logical, intent (in) :: adj, add
    real, dimension (:)  :: x, y
    if( .not. add) then 
       if( adj) then
          x = 0.
       else 		    
          y = 0.
       end if
    end if
  end subroutine adjnull
end module adj_mod
