Running and Compiling Duo
=========================


Downalod Duo from github_

.. _github: https://github.com/Trovemaster/Duo


Duo is provided both as source code and as compiled executables for Linux and Microsoft Windows.
Using the executables is the easiest way to run Duo. Duo works from the command line (also known as the ``terminal``
or ``command prompt``).
A Duo input is a plain text input file which is fed into the program with a command of the kind
::

     ./duo.exe < ./inputs/my_input.inp > my_output.txt

If you accidentally start Duo without specifying any input nothing will happen and you will be temporarely stuck;
to terminate Duo press C while holding down the Ctrl key.

Please note that Duo is still in active development and new versions with bug fixes and new functionalities are expected to
appear regularly. If you have found a bug or you would like to make a comment please do not hesitate to contact the authors
(contact details are reported in the first page of this manual).

The Duo examples can be obtained from github here 
git clone https://github.com/Trovemaster/Duo/tree/MOLPRO/examples

On Windows 
^^^^^^^^^^

On Windows, download the executable file duo_win.exe_ and the batch file run_duo.bat_
from the githubexe_ repository. To run from a window, modify run_duo.bat (change the names of the input and output files) 
and double click on it. Form CMD, navigate to the folder with the input file and execute
::
     
     duo_win.exe < input_file  > output_file

     
A Windows (DOS) batch file has the following format:
::
    
    REM example of a Windows batch file
    duo_win.exe < BeH_Koput_01.inp > BeH_Koput_01.inp

.. _duo_win.exe: https://github.com/Trovemaster/Duo/blob/MOLPRO/executables/duo_win.exe


.. _run_duo.bat: https://github.com/Trovemaster/Duo/blob/MOLPRO/executables/run_duo.bat


Running on Linux (or Mac)
^^^^^^^^^^^^^^^^^^^^^^^^^

Copy the executable file into into a working folder. A precompiled file can be obtained from githubexe_.

.. _githubexe: https://github.com/Trovemaster/Duo/blob/MOLPRO/executables/

Navigate to the folder with the input_file and execute from the command line: 
::
     
     ./duo_linux_64.exe < input_file  > output_file



Compilation
^^^^^^^^^^^

You may need to re-compile Duo if the provided executables do not work on your system or, for example, if you want to make modifications to the program.
Duo makes use of some Fortran 2003 features and therefore requires a compiler with (at least partial) support for Fortran 2003.
At the time of writing (July 2015) there are two freely-available Fortran compilers (for Windows, Linux and OS X) which can be used for compiling Duo, 
namely gfortran and g95.
Lists of Fortran compilers can be found on the Internet_.
The Fortran 95 Windows compiler Silverfrost_ FTN95 v.7.20 is also available for free for personal use, 
but its support for Fortran 2003 is very incomplete and Duo will not work with this compiler at this time.

.. _Internet:  http://fortranwiki.org/fortran/show/compilers


.. _Silverfrost: http://www.silverfrost.com/32/ftn95/ftn95\_personal\_edition.aspx 


Duo has been tested with the Intel Fortran Compiler v. 12.1 under Linux and Windows 8, with the Portland Group Fortran compiler v. 13.1, 
with the NAG Fortran compiler v 5.2 under Linux and Windows 8, with g95 v. 0.94 under Windows 8 and
with gfortran v.4.9.2 under Windows 8, Linux and OS X.





To compile, edit the makefile and compile by running 
::

     make



