

# Style guide


## Indentation

Use spaces, not tabs.


## Do not use "<<<<" or ">>>>" anywhere in the code

These are conflict resolution markers.


## Private vs. public

Everything should be private except names explicitly declared as public

This speeds up compilation and avoids bugs and name clashes. It gives more modular
and thus better code.


## Never use F77 style common blocks

They are like nuclear waste.


## Never put several commands on one line

Hard to read and inconvenient for debugging (when you need to uncomment one of the commands):

```
Ef=-Ef; Eff=-Eff; Effff=-Effff; Egff=-Egff
```

It is also inconvenient from the version control point of view: more commands
on one line increase risk of conflicts during merges.


## Use module-wide implicit none

Example:

```
module birefring

  use this
  use that

  implicit none

  ...
```

Never use `implicit.h`. After you spend a day hunting a bug created
by using implicit you know why.


## Do not leave commented code behind

Also if you replace a routine by a better routine, remove the less better routine.


## No program-specific CPP filters

They mean that the library is not general. Program specificity has to go to the
interfaces.


## Use self-explaining variable names

Also never ever reuse variable names just to save a declaration.


## Module naming

Fortran modules and files containing modules should have the same name.

Also we should not let Dalton (or DIRAC) restrict our naming. We had examples
where we chose a less ideal name because the better name could upset some
developers of Dalton.  It is our library and we can dictate naming. Name
conflicts can be resolved in interfaces.
