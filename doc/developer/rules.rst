
Tentative Rules for Developers
==============================

Short-version
-------------

#. First analyze the problem, then design the code (data and algorithm
   structures), prepare test suite. Last, write the code.

#. Make everything as simple as possible.

#. Do your best to prepare a readable document.

Long-version
------------

#. First of all, please write explicitly what you would like to implement in
   ``doc``! Describe your idea using formulas and/or words. Then translate them
   into algorithms and data structure. Please do write what
   objects/types/variables you will define and their corresponding public and
   private functions (including detailed descriptions of the input and output
   arguments). It would be better if you could write down the framework of your
   implementation using figures. Please also write down the limitations or
   risks of your code, for instance, does it stable or have some numerical
   error? If yes, how to prevent or how to know if the results are reasonable?

   In this stage, you may refer to some rules in object-oriented programming
   (OOP). For instance, when you define a module/class etc.:

   #. it should be open for extension but closed for modification (Open Closed
      Principle, OCP),
   #. subclasses should be substitutable for their base classes (Liskov
      Substitution Principle, LSP),
   #. depend upon abstractions, do not depend upon concretions (Dependency
      Inversion Principle, DIP),
   #. many client specific interfaces are better than one general purpose
      interface (Interface Segregation Principle, ISP),
   #. In other words: low coupling, high cohesion, open for extension, and
      closed for changes (from "Developing Chemical Information Systems: An
      Object-Oriented Approach Using Enterprise Java", Fan Li).

#. Write the codes. During this stage, we would be happy if you could:

   #. write comments (in english, one line for each 10-20 line of codes at
      least),
   #. try to use descriptive names for your classes and methods,
   #. do your best to avoid global variables,
   #. try to re-use code and try to use libraries,
   #. use ``fixme`` or ``FIXME`` to identify the codes needed to be modified or
      fixed later.

#. This is very important, and should be considered and implemented during the
   aforementioned two steps:

   Always provide a test suite for each function/subroutine/module
   etc., unless you are 100% sure what you did is right. Integration testing
   will also be required in some cases.
