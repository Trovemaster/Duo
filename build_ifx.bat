@echo off
setlocal EnableExtensions EnableDelayedExpansion

rem Build Duo with Intel ifx in an isolated directory.  Keeping all .mod files
rem here prevents an older module in the source directory from being selected.
set "SOURCE_DIR=%~dp0"
set "BUILD_DIR=%SOURCE_DIR%build\ifx"

where ifx >nul 2>nul
if errorlevel 1 (
    set "ONEAPI_SETVARS=%ProgramFiles(x86)%\Intel\oneAPI\setvars.bat"
    if exist "!ONEAPI_SETVARS!" call "!ONEAPI_SETVARS!" intel64 >nul
)

where ifx >nul 2>nul
if errorlevel 1 (
    echo ERROR: ifx was not found.
    echo Run this file from an Intel oneAPI command prompt, or install the
    echo Intel oneAPI HPC Toolkit and ensure that setvars.bat is available.
    exit /b 1
)

if not exist "%BUILD_DIR%" mkdir "%BUILD_DIR%"
if errorlevel 1 exit /b 1

rem Remove only products made by this script.  Source-tree .mod files are
rem deliberately ignored because compilation takes place inside BUILD_DIR.
del /q "%BUILD_DIR%\*.obj" "%BUILD_DIR%\*.mod" "%BUILD_DIR%\duo.exe" 2>nul

rem Override these flags before invoking the script if required, for example:
rem   set DUO_IFX_FLAGS=/Od /debug:full /traceback /check:all /Qopenmp
if not defined DUO_IFX_FLAGS set "DUO_IFX_FLAGS=/O3 /Qopenmp /traceback"
rem Explicit libraries are used instead of /Qmkl because some Windows ifx
rem installations accept that shortcut but do not pass the MKL libraries to
rem link.exe.  Sequential MKL is deterministic and Duo still uses OpenMP.
if not defined DUO_MKL_LIBS set "DUO_MKL_LIBS=mkl_intel_lp64.lib mkl_sequential.lib mkl_core.lib"

rem This is a topological order of the Fortran modules.  In particular, the
rem merged intensity pipeline requires:
rem intensity_filters feeds dipole and quadrupole independently.  The rank-1
rem implementation still requires magnetic_dipole before dipole.
set "SOURCES=accuracy timer lapack input symmetry Lobatto atomic_and_nuclear_data me_numer functions grids diatom refinement F1_hyperfine F1_intensity intensity_filters magnetic_dipole dipole quadrupole header_info RWF intensity duo"

pushd "%BUILD_DIR%"
if errorlevel 1 exit /b 1

set "OBJECTS="
for %%S in (%SOURCES%) do (
    echo [ifx] Compiling %%S.f90
    ifx /nologo /c %DUO_IFX_FLAGS% /module:"%BUILD_DIR%" /I"%BUILD_DIR%" /object:"%BUILD_DIR%\%%S.obj" "%SOURCE_DIR%%%S.f90"
    if errorlevel 1 (
        echo ERROR: compilation stopped at %%S.f90
        popd
        exit /b 1
    )
    set "OBJECTS=!OBJECTS! "%BUILD_DIR%\%%S.obj""
)

echo [ifx] Linking duo.exe
ifx /nologo %DUO_IFX_FLAGS% /exe:"%BUILD_DIR%\duo.exe" !OBJECTS! %DUO_MKL_LIBS%
if errorlevel 1 (
    echo ERROR: link failed.
    popd
    exit /b 1
)

popd
echo Build completed: "%BUILD_DIR%\duo.exe"
exit /b 0
