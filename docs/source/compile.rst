Running and Compiling Duo
=========================


Downalod Duo from github_

.. _github: https://github.com/Trovemaster/Duo


Duo is provided both as source code and as compiled executables for Linux and Microsoft Windows.
Using the executables is the easiest way to run Duo. Duo works from the command line (also known as the _terminal_
or _command prompt_).
A Duo input is a plain text input file which is fed into the program with a command of the kind
::

     ./duo.exe < ./inputs/my_input.inp > my_output.txt

If you accidentally start Duo without specifying any input nothing will happen and you will be temporarely stuck;
to terminate Duo press C while holding down the Ctrl key.

Please note that Duo is still in active development and new versions with bug fixes and new functionalities are expected to
appear regularly. If you have found a bug or you would like to make a comment please do not hesitate to contact the authors
(contact details are reported in the first page of this manual).





To compile, edit the makefile and compile by running 
::

     make

Prepare an input file `filename.inp` (see example on github) and run using 
::

     exocross.exe <filename.inp >filename.out





