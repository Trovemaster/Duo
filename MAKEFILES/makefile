goal:   duo.x

tarball:
	tar cfz duo$(PLAT).tar.gz makefile *.f90
        
checkin:
	ci -l Makefile *.f90


#FPATH = 

PLAT = 2210_gf

##FOR  = ifort   # Fortran compiler 
FOR = gfortran

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
####FFLAGS = -O0 -fpe0  -fltconsistency -stand f03 -check all -warn all -traceback -fp-stack-check  # debugging options

''FFLAGS = -O3 -ip -qopenmp -mkl=parallel ##  -xHost -fast

##FFLAGS = -C -check bounds -g  -gen-interfaces -warn interfaces  -check arg_temp_created -prof-value-profiling=all -warn all
##FFLAGS = -O3 -ip -openmp # no optimization -- fast compilation

FFLAGS = -W  -fbounds-check -pedantic-errors -std=f2003 -Wunderflow -O0 -fbacktrace -g -Wextra

#-Wall

#ARPACK =  ~/libraries/ARPACK/libarpack_omp_64.a

LAPACK = 

#LAPACK = -mkl=parallel 
#LAPACK = -mkl=parallel -static

LIB     =   $(LAPACK)

###############################################################################

OBJ = atomic_and_nuclear_data.o grids.o accuracy.o lapack.o timer.o input.o diatom.o refinement.o functions.o  symmetry.o dipole.o header_info.o Lobatto.o
#compilation_details.o 

duo.x:	$(OBJ) duo.o
	$(FOR) -o j-duo$(PLAT).x $(OBJ) $(FFLAGS) duo.o $(LIB)

duo.o:	duo.f90 $(OBJ) 
	$(FOR) -c duo.f90 $(FFLAGS)

grids.o:	grids.f90 accuracy.o input.o Lobatto.o
	$(FOR) -c grids.f90 $(FFLAGS)

diatom.o:	diatom.f90 accuracy.o input.o lapack.o functions.o symmetry.o atomic_and_nuclear_data.o Lobatto.o
	$(FOR) -c diatom.f90 $(FFLAGS)

refinement.o:	refinement.f90 accuracy.o input.o lapack.o diatom.o
	$(FOR) -c refinement.f90 $(FFLAGS)

functions.o:	functions.f90 accuracy.o input.o lapack.o
	$(FOR) -c functions.f90 $(FFLAGS)

dipole.o:	dipole.f90 accuracy.o input.o lapack.o diatom.o 
	$(FOR) -c dipole.f90 $(FFLAGS)

accuracy.o:  accuracy.f90
	$(FOR) -c accuracy.f90 $(FFLAGS)

symmetry.o:  symmetry.f90
	$(FOR) -c symmetry.f90 $(FFLAGS)

lapack.o:  lapack.f90 accuracy.o timer.o 
	$(FOR) -c lapack.f90 $(FFLAGS)

timer.o:  timer.f90
	$(FOR) -c timer.f90 $(FFLAGS)

input.o:  input.f90
	$(FOR) -c input.f90 $(FFLAGS)

#compilation_details.o: compilation_details.f90
#	$(FOR) -c compilation_details.f90 $(FFLAGS)

header_info.o:  accuracy.o
	$(FOR) -c header_info.f90 $(FFLAGS)

atomic_and_nuclear_data.o:  atomic_and_nuclear_data.f90
	$(FOR) -c atomic_and_nuclear_data.f90 $(FFLAGS)

Lobatto.o: Lobatto.f90 timer.o
	$(FOR) -c Lobatto.f90 $(FFLAGS)

clean:
	rm -f $(OBJ) *.mod *__genmod.f90 duo.o duo_test_0* eigen_vectors.chk eigen_vib.chk Bob-Rot_centrifugal_functions.dat  _Lp__functions.dat       Spin-Orbit.dat               Spin-spin_functions.dat Dipole_moment_functions.dat Potential_functions.dat  Spin-rotation_functions.dat Spin-spin-o__non-diagonal__functions.dat

