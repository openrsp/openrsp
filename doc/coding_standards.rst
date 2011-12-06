

Coding standards
================


Indentation
-----------

2 blanks after module, subroutine, function, program (and contains)

The corresponding 'contains' at same indentation.

3 blanks after if, do, type

No tabs, really. Tab enthusiasts will tell you that they are great
but they are a big pain if you copy paste code or want to comment code out.

Don't use "<<<<" or ">>>>" anywhere in the code
-----------------------------------------------

These are conflict resolution markers.


Private vs. public
------------------

Everything should be private except names explicitly declared as public

This speeds up compilation and avoids bugs and name clashes. It gives more modular
and thus better code.


Write end module, end type, end function, end subroutine, end if, end do, without name
--------------------------------------------------------------------------------------

subroutine foo
end subroutine
instead of
subroutine foo
end subroutine foo

Reason: It is unnecessary and annoying when you want to rename things, you have to rename
them in two places.

Also no single-word variants endmodule, endtype, endfunction, endsubroutine, endif or enddo.


Never use F77 style common blocks
---------------------------------

They are like nuclear waste.


Don't put several commands on one line
--------------------------------------

Hard to read and inconvenient for debugging (when you need to uncomment one of the commands):
Ef=-Ef; Eff=-Eff; Effff=-Effff; Egff=-Egff


It is also inconvenient from the version control point of view: more commands
on one line increase risk of conflicts during merges.


Use module-wide implicit none
-----------------------------

Example::

  module birefring
  
    use this
    use that
  
    implicit none
  
    ...

Never use implicit.h. After you spend a day hunting a bug created
by using implicit you know why.


Do not leave commented code behind
----------------------------------

Also if you replace a routine by a better routine, remove the less better routine.


No program-specific CPP filters
-------------------------------

They mean that the library is not general. Program specificity has to go to the
interfaces.


Use self-explaining variable names
----------------------------------

Also never ever reuse variable names just to save a declaration.


Module naming
-------------

Modules and files containing modules should have the same name.

Also we should not let Dalton (or DIRAC) restrict our naming. We had examples
where we chose a less ideal name because the better name could upset some
developers of Dalton.  It is our library and we can dictate naming. Name
conflicts can be resolved in interfaces.
