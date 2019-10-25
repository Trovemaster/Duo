Computing energy levels
=======================

In the following we present all keywords and options relevant to the calculations of energy levels.

Calculation setup
^^^^^^^^^^^^^^^^^

* Atoms:  defines the chemical symbols of the two atoms.
Example:
::

    atoms Na-23 H-2

specifies the 23NaD diatomic. Duo includes an extensive database of atomic properties (atomic masses,
nuclear spins, isotopic abundances and other quantities) and will use the appropriate values
when required. The database should cover all naturally-occurring nuclei as well as all radioactive ones
with a half-life greater than one day and is based on the AME2012 and
NUBASE2012 databases. Each atom should be specified by its chemical symbol,
a hyphen (minus sign) and the atomic mass number, like in the example above.
Atomic masses will be used, which is generally the most appropriate choice unless one is explicitely including
non-adiabatic corrections. The hydrogen isotopes deuterium and tritium can also
be optionally specified by the symbols `D` and `T`.


The atomic mass number can be omitted, like in the following example:
::

     atoms Li F

In this case Duo will use the most-abundant isotopes (7Li and 19F in the example above) or, for
radioactive nuclei not naturally found, the longest lived one. For example 
::

     atoms Tc H

selects for technetium the isotope 97Tc, which is the longest-lived one. A few nuclides in the database are nuclear metastable isomers,
i.e. long-lived excited states of nuclei; these can be specified with a notation of the kind
::

     atoms Sb-120m H


In the example above the radioactive isotope of antimony 120mSb is specified (and hydrogen).
Another example
::

     atoms Sc-44m3 H

specifies the scandium radioactive isotope 44m3Sc (and hydrogen).






