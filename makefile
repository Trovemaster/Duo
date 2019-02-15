goal:   duo

tarball:
	tar cfz duo$(PLAT).tar.gz makefile *.f90
        
checkin:
	ci -l Makefile *.f90


#FPATH = 

EXE = j-duo-v218.v1.x

FOR  = ifort   # Fortran compiler 
##FOR = gfortran  

# Lorenzo Lodi  ---- meaning of some flags used by the Intel fortran compiler
#  see (e.g.) file:///opt/intel/composer_xe_2011_sp1.8.273/Documentation/en_US/compiler_f/main_for/index.htm
# -O0 removes compiler optimization of the code. Code is slower in execution but compilation is much faster.
# -fltconsistency disables some optimisations which may lead to very small numerical roundoff errors
# -stand f95 gives warnings if non-Fortran95-compliant features are used
# -stand f03 gives warnings if non-Fortran03-compliant features are used
# -check all enables the following run-time checks:
#                    arg_temp_created  => Determines whether checking occurs for actual arguments before routine calls.
#                    bounds            => checks if accessing an out-of-bounds array element (e.g., x(0) when x is defined as x(1:NMAX)
#                    format            => mismatch in the output format, e.g. write(*,'(F10.5)') i , and i is an integer
#                    output_conversion => output error, e.g. the number won't fit in the given spaces
#                    pointer           => checks for disassociated/uninitialized pointers
#                    uninit            => checks for uninitialized variables
# -warn all enables the following compile-time diagnostic checks:
#                    alignments        => data that is not naturally aligned
#                    declarations      => undeclared names
#                    errors            => warnings are changed to errors
#                    general           => various things
#                    ignore_loc        => %LOC is stripped from an actual argument. (?)
#                    interfaces        => compiler checks the interfaces
#                    stderrors         => Fortran standard violations are changed to errors.
#                    truncated_source  => source exceeds the maximum column width
#                    uncalled          => statement function is never called
#                    unused            => declared variables that are never used.
#                    usage             => questionable programming practices.
# -traceback ***GOLDEN OPTION***  when a run-time error occurs it'll tell you the line in the source code responsible
# -fp-stack-check  stops the program immediately when a NaN (Not-a-numeber) error occurs
# Other options:
# -O3 most agressive general optimization setting
# -ip enables additional interprocedural optimizations
# -g generate full debugging information in the object file
#  -C compile only (=make .o files), but do not link .o files to produce executable
#  -prof-value-profiling=all   All value profile types are enabled and value profiling is performed.

#NOTE: -fpe0 will stop on floating-point exceptions. Do not use this flag because LAPACK makes use of divide by zero etc.
#
##FFLAGS = -O0 -fpe0  -fltconsistency -stand f03 -check all -warn all -traceback -fp-stack-check  # debugging options

FFLAGS = -O3 -ip -openmp -mkl=parallel # -xHost -fast

##FFLAGS = -C -check bounds -g  -gen-interfaces -warn interfaces  -check arg_temp_created -prof-value-profiling=all -warn all
##FFLAGS = -O3 -ip -openmp # no optimization -- fast compilation
##FFLAGS = -W -Wall -fbounds-check -pedantic-errors -std=f2003 -Wunderflow -O0 -fbacktrace -g -Wextra


#ARPACK =  ~/libraries/ARPACK/libarpack_omp_64.a

LAPACK = -mkl=parallel -static
#LAPACK = -mkl=parallel -static

LIB     =   $(LAPACK)

###############################################################################

OBJ = grids.o accuracy.o lapack.o timer.o input.o diatom.o refinement.o functions.o  symmetry.o dipole.o header_info.o atomic_and_nuclear_data.o  Lobatto.o

diatom.o: symmetry.o functions.o input.o lapack.o Lobatto.o timer.o atomic_and_nuclear_data.o accuracy.o
dipole.o: timer.o accuracy.o diatom.o symmetry.o
duo.o: header_info.o diatom.o accuracy.o refinement.o timer.o dipole.o
functions.o: accuracy.o timer.o
grids.o: accuracy.o Lobatto.o
header_info.o: accuracy.o
lapack.o: accuracy.o timer.o
Lobatto.o: accuracy.o timer.o
refinement.o: timer.o accuracy.o diatom.o
symmetry.o: accuracy.o
timer.o: accuracy.o

# clear internal suffix rules
.SUFFIXES:
# specify our own suffix rules
.SUFFIXES: .f90 .o
.f90.o:
	$(FOR) -c -o $@ $< $(FFLAGS)

duo:	$(OBJ) duo.o
	$(FOR) -o $(EXE) $(OBJ) $(FFLAGS) duo.o $(LIB)

clean:
	rm -f $(OBJ) *.mod *__genmod.f90 duo.o duo_test_0* eigen_vectors.chk eigen_vib.chk Bob-Rot_centrifugal_functions.dat  _Lp__functions.dat       Spin-Orbit.dat               Spin-spin_functions.dat Dipole_moment_functions.dat        Potential_functions.dat  Spin-rotation_functions.dat Spin-spin-o__non-diagonal__functions.dat

test:  test_duo.f90
	$(FOR) -O0 test_duo.f90 -o duo_test.exe
