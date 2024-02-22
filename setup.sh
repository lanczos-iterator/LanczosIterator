#!/bin/sh

#for unix-like systems  (Linux, macOS, etc.)
#maybe can be adapted for unix-like environment in Windows (MSYS2?)

#choices for these for some common systems are given below,
# otherwise specify these here
FC=               # fortran compiler
CC=               # C compiler
OPT=              # fortran compiler options
NO_FORTRAN_MAIN=  # needed by some Fortran compilers if the main program is C 
DEBUG=            # for compiling a version useful for debugging.

#uncomment for debug versions 
#DEBUG="yes"

if [ -z "$DEBUG" ] ; then
    LIBNAME=lanczos
else
    LIBNAME=lanczos_debug
fi

#preset options for some systems
GNU=
MACOS=
INTEL=
INTEL_CLASSIC=
NVIDIA=

#uncomment one of these that matches your compiler and Lapack library
#GNU="yes"             #gfortran, gcc, + gfortran-compiled LAPACK/BLAS
#MACOS="yes"           #gfortran+ gcc + Apple Accelerate (needs wrapped ZDOTC) 
#INTEL="yes"           #Intel OneAPI ifx, ifc compilers, + MKL library
#INTEL_CLASSIC="yes"   #Intel ifort, icc compilers, + MKL library
#NVIDIA="yes"          # NVIDIA HPC compilers + libraries from HPC SDK


if [[ -n "$MACOS" ]] ; then
    FC=gfortran
    CC=gcc
    LAPACK="-framework Accelerate"
    ZDOTC='-ff2c -fno-second-underscore'       
fi

if [[ -n "$INTEL" ]] ; then
    FC=ifx
    CC=icx
    LAPACK="-qmkl=sequential"  
    #LAPACK="-qmkl=parallel"   #use with OpenMP
    NO_FORTRAN_MAIN="-nofor-main"
fi

if [[ -n "$INTEL_CLASSIC" ]] ; then
    FC=ifort
    CC=icc
    LAPACK="-qmkl=sequential"  
    #LAPACK="-qmkl=parallel"   #use with OpenMP
    NO_FORTRAN_MAIN="-nofor-main"
fi

if [[ -n "$NVIDIA" ]] ; then
    FC=nvfortran
    CC=nvc
    LAPACK="-llapack -lblas"  
    #LAPACK="-qmkl=parallel"   #use with OpenMP
    NO_FORTRAN_MAIN="-Mnomain"
fi

if [[ -n "$GNU" ]] ; then
    FC=gfortran
    CC=gcc
    LAPACK="-llapack -lblas"
    NO_FORTRAN_MAIN=
fi

if [ -z "$FC" ] ; then
    echo "*** No FORTRAN compiler has yet been specified ***" ;
    echo "*** Edit setup.sh to fix this ***"
    exit
fi

if [[ -n "$DEBUG" ]] ; then
    if [[ "$FC" == "gfortran" ]] ; then
	OPT2="-fcheck=all"
    elif [[ "$FC" == "ifx" || "$FC" == "ifort" ]] ; then
	OPT2="-check all"
    else
	OPT2=
    fi
fi

LIB=$LIBNAME.a
#LIB=$LIBNAME_$FC.a

#remove pre-existing versions in case you are rebuilding
if [ -f $LIB ] ; then
    rm $LIB
fi

cd lib
$FC $OPT -c vector_m.f90 $OPT2
$FC $OPT -c matrix_algebra_m.f90 $OPT2 
$FC $OPT -c tmatrix_m.f90 $OPT2
$FC $OPT -c lanczos_m.f90 $OPT2
$FC $OPT -c *.f90  -I . $OPT2
ar -q $LIB *.o
rm *.o
rm *.mod
mv $LIB ..
cd ..
cd C_interface
$FC $OPT -c blas_settings_wrapper.f90 $OPT2
$FC $OPT -c rs_lanczos_wrapper.f90 $OPT2 
$FC $OPT -c ch_lanczos_wrapper.f90 $OPT2
ar -q ../$LIB *.o
rm *.o
cd ..
ranlib $LIB

if [ -z "$LAPACK" ] ; then
    echo "to build test programs, specify path to LAPACK/BLAS"
    exit
fi

if [ -f test_zdotc ] ; then
    rm test_zdotc
fi
echo " building test program test_zdotc"
$FC test_programs/test_zdotc.f90 $LIB $LAPACK -o test_zdotc

echo " building fortran test program test_lanczos_f"

if [ -f test_lanczos_f ] ; then
    rm test_lanczos_f
fi
$FC $OPT test_programs/test_lanczos.f90 $LIB $LAPACK  $OPT2 -o test_lanczos_f 

if  [ -z "$CC" ] ; then
    echo "to build the C test program, edit setup.sh to specify your C compiler"
    exit
fi

echo " building C test program test_lanczos_c "
$CC -c test_programs/test_lanczos.c -I include

if [ -f test_lanczos_c ] ; then
    rm test_lanczos_c
fi

# Linking of combined C and Fortran code is done by the FC compiler
# May need $NO_FORTRAN_MAIN flag when the main program is a C program

$FC test_lanczos.o $LIB $LAPACK $NO_FORTRAN_MAIN  -o test_lanczos_c
rm test_lanczos.o
