# Lanczos_Iterator v3.3
##A "ROBUST" Lanczos diagonalization package (in Fortran95, with C99 interface)

(c) F. D. M. Haldane (2014-2024),  Princeton University,    v3.3   2024-02-20  

Purpose: given an efficent implementation "matvec" of (sparse) matrix-vector multiplication

```
             v_out(i) = sum_j matrix(i,j) * v_in(j)
```
where matrix is a large nxn Hermitian matrix (real symmetric or complex Hermitian):

  * (1) determine the maximum eigenvalue precision ("resolution") possible, given
      that matvec is implemented in (64 bit) floating point as opposed to exact arithmetic.

  * (2) obtain the lowest 1 << n_v << n eigenvalues AND corresponding eigenvectors
      (with their "variances" sqrt(<evec(i)|(M-eval(i))^2|evec(i)>), to
      the **MAXIMUM POSSIBLE PRECISION** consistent with the matvec-procedure resolution.
      No user input about desired precision is needed, and eigenvalue
      degeneracies are handled correctly.   The standard Lanczos algorithm is supplemented with a
      novel "Lanczos polishing" technique to reduce the
      variance (residual norm) of eigenvectors down to the floating-point resolution limit.

##  Applications:

*    Obtain maximally-accurate low-lying eigenvalues  **and the corresponding eigenvectors** of _huge_ sparse Hermitian
     matrices (dimensions of order tens of millions or more) that describe the Hamiltonians of interacting quantum many-body systems of
     interest in condensed matter physics  (fractional quantum Hall systems, fractional Chern insulators, Hubbard models, spin chains, etc.),
     with full resolution of degeneracies, etc.    The method is practical when the matrices are sparse, so the matrix-times-vector multiplication
     operations can be efficiently implemented. (Full-matrix diagonalization methods are limited to matrix dimensions of
     at most 100,000 with todays (2024) most powerful computers.)

*    Whether there are any other areas needing this tool is an open question.





##Usage:

Note: Both Fortran and C interfaces are described: "_TRUE_" and "_FALSE_" will mean
 values "1" and "0" \(C), "true" and "false" (C with #include <stdbool.h>) and
 ".true." and ".false." (Fortran)

 Use driver  "ch_lanczos" to obtain eigenvalues/vectors 1:nevec of a complex Hermitian matrix (64-bit reals)
 using an "all-at-once" method

  *  C  (C99):
  
```
    #include <complex.h>
    #include "lanczos.h"
    int n, nevec, maxstep, report, seed;
    void matvec(int n, double complex *vec);
       // double complex *matrix;   /* a pointer to the full matrix */
       // void matvec(int n, double complex *vec) {
       //    double complex *vec_out = calloc(n, sizeof(double complex));
       //    for (int i = 0; i < n; i++) {
       //        for (int j = 0; j < n; j++)
       //            vec_out[i] += matrix[i][j] * vec_in[j];
       //    }
       //    memcpy(vec, vec_out, n * sizeof(double complex));
       //    free (vec_out);
       // }
    scanf("%d %d %d %d %d", &n, &nevec, &maxstep, &report, &seed);
    double *eval = (double *) malloc(nevec * sizeof(double))
    double complex *evec = (double complex *) malloc(nevec * n *
                               sizeof(double complex));
    double *variance = (double *) malloc(nevec * sizeof(double));
    double resolution;
    ch_lanczos_c(n, nevec, *eval, *evec, *variance, &resolution,
                 *matvec, maxstep, report, seed, NULL, 0);
```

  * Fortran (Fortran95):

```
      implicit none
      integer :: n, nevec, maxstep, report, seed
      interface
         subroutine matvec(n, vec_in, vec_out, add)
	   use matrix_m ! module provides "complex (kind(1.0d0)):: matrix(n,n)"
	   integer, intent(in) :: n
	   complex (kind(1.0d0)), intent(in) :: vec_in(n)
	   complex (kind(1.0d0)), intent(inout) :: vec_out(n)
	   logical, intent(in) :: add
           ! if (add) vec_out = matmul(matrix,vec_in)
           ! if (.not. add) vec_out = vec_out + matmul(matrix,vec_in)
           ! return
         end subroutine matvec
      end interface
      real (kind(1.0d0)) :: resolution
      real (kind(1.0d0)), allocatable :: eval(:), variance(:)
      complex (kind(1.0d0)), allocatable :: evec(:,:)
      read *, n, nevec, maxstep, report, seed
      allocate (eval(nevec), evec(n,nevec), variance(nevec))
      call ch_lanczos(n, nevec, eval, evec, variance, resolution,&
                      matvec, maxstep, report, seed)

```
Obviously, a sparse-matrix implementation of matvec that is more efficient
than the reference implementations shown above (which use the full matrix)
should be used.

The method requires that all eigenvectors of eigenvalues lower than a target
eigenvalue have already been found.  The code also supports more flexible
calling methods that can find eigenvectors  n1,..,n2 <= nevec when eigenvectors
1,...,n1-1 have previously been found, including "one-at-a-time" mode n1 = n2.
The above is for complex Hermitian matrices defined by MATVEC; a real-symmetric
variant "rs_lanczos" is also supplied.

The input parameters are:
    
        n        length of eigenvector (matrix is  n X n)
        nevec    lowest nevec eigenvalues + eigenvectors are requested
        maxstep  maximum allowed size of Krylov subspace: set maxstep = 0
                 for no limit.
        report   set report = 0 for no progress reports; report > 0
                 outputs a summary of the Lanczos process each time
                 an eigenvector is produced.  Also e.g. report=50 gives
                 a progress report after every 50'th Lanczos step.
        seed     an integer used to generate a random starting vector


The "more flexible" extra calling arguments modify the above as follows: 

in C, the call to  ch_lanczos_c(....,seed, NULL, 0)  becomes

```
     int range[2] = { n1, n2 };
     bool nocheck;
     ch_lanczos_c(....,seed, range, nocheck); 
```

in Fortran  the call to ch_lanczos(....,seed) becomes

```
     integer :: range(2) = (/ n1, n2 /)
     logical  :: nocheck
     call ch_lanczos(...,seed, range, nocheck)         

```
Here eigenvalues in the range [max(1,n1):min(nevec,n2)] are being requested.
It is assumed that eigenvalues and eigenvectors in the range [1:n1-1] are present frm a prior
calculation.  These will be checked (by recomputing their eigenvalues and variances, then
comparing with values in eval and variance) unless nocheck = _TRUE_.

This procedure is particularly useful when full **EIGENVECTOR** accuracy is needed.
Accurate sequential computation of eigenvectors, in order of increasing eigenvalue,
is central to the algorithm. A new "polishing" algorithm to finish refining the
eigenvector after "standard" Lanczos yields no futher improvement is implemented.
All results for eigenvalues and variance (residual norm) are intrinsically validated
with matvec itself, without reference to the Lanczos procedure that obtained them.

==============================================================================

### Installation and building the (static) library:

To use the code:  build a (static) library lanczos.a with your fortran compiler
(or run the shell script setup.sh which does it for you)   You will need to know
how to link to the BLAS and LAPACK libraries.   Here we assume this is done with
"-llapack -lblas", and your compilers are FC and CC (plus a special option OPT),
where

                               (FC,        CC,    OPT                )
            gcc   (GNU):       (gfortran,  gcc,   ""                 )
            Intel (OneAPI):    (ifx,       icx,   "-nofor-main"      )
            Intel (Classic):   (ifort,     icc,   "-nofor-main"      )
            NVIDIA HPC:        (nvfortran, nvc,   "-Mnomain"         )
            Flang (LLVM):      (flang,     clang, "-fno-fortran-main")
            Flang (Classic):   (flang,     clang, "-fno-fortran-main")


### Building Lanczos using the shell script "setup.sh"

The easiest way to build the Lanczos library (plus some test programs) is:

```
tar -xvzf lanczos_iterator-v3.3.tgz
cd lanczos_iterator-v3.3
```

Then edit the shell script  setup.sh to identify your Fortran compiler, library call for
LAPACK + BLASS, and (if you will use C ) your C compiler.   Some common choices labeled
"GNU", "Intel OneAPI", "Intel Classic", "MacOS", "NVIDIA" provide preselected values when
uncommented. Uncomment the line #DEBUG="yes" to build the library
with debugging options.

Then just run

```
./setup.sh
```

This will produce the static library "liblanczos.a" and test programs "test_zdotc", "test_lanczos_f" (Fortran) and
"test_lanczos_c" \(C).    The code for the test programs is in subdirectory "test_programs".
The code in test_lanczos.f90 and test_lanczos.c may be helpful in showing how to use the
drivers  **ch_lanczos**  and **rs_lanczos**.s

### Building without the shell script:

To build the Lanczos library without using "setup.sh",  substituting appropriate values  (FC, CC, OPT) for
the Fortran compiler, C compiler, and a special option OPT (sometimes needed when
the main program is in C rather than Fortran), the workflow is:

```
cd lib
FC -c vector_m.f90 matrix_algebra_m.f90 lanczos_m.f90 tmatrix_m.f90
FC *.f90 -c
ar -q ../liblanczos.a *.o
cd ../C_interface
FC *.f90 -c
ar -q ../liblanczos.a *.o
cd ..
ranlib liblanczos.a
```

To build the test programs, it will be assumed that the LAPACK/BLAS libraries are accessed
with "`-llapack -lblas`".   If this is not the case, find out how to call them on your system
and use the correct call.

```
FC test_programs/test_zdotc-f90 -llanczos -lblas
FC test_programs/test_lanczos.f90 -llanczos -llapack -lblas
CC test_programs/test_lanczos.c -c -I include
FC test_lanczos.o OPT -llanczos -llapack -lblas
```
(This assumes that liblanczos.a
is either in the current directory, or in $LIBRARY_PATH.)


Note how after the C main program in test_lanczos.c is compiled, it is linked to the Lanczos library using the _Fortran_
compiler, which may need the option OPT to tell it that the main program is not Fortran.

### Using the test programs:


#### Check your BLAS (Basis Linear Algebra Subprograms) library

It is preferable to use a processor-optimized version of the BLAS 
library (a companion of LAPACK) for linear algebra operations on  large (dimension  _n_) vectors. These are
often vendor-supplied (Intel's MKL, Apple's Accelerate, NVIDIA's NVPL, AMD's ACML, _etc._).  However they may have
a compatibility issue with certain Fortran compilers (in particular gfortran)  involving complex  Fortran  functions,
in particular ZDOTC. If so, this can be fixed by using an option to "wrap" ZDOTC; alternatively there is  an option to
use native Fortran vector procedures instaed of BLAS.

To check compatibility of your ZDOTC implementation with your compiler:

```
./test_zdotc
```

This offers testing with or without the wrapper.

If this shows zdotc works without a "wrapper", you are fine.   If it only
works with the wrapper, in you will need to call "set_wrap_zdotc(_TRUE_)".
You can instead replace use of the BLAS procedures with native Fortran95 vector
procedures by calling "set_use_blas(FALSE).  The only reason to use BLAS
rather than native procedures is if they have been optimized for your processor 
and are faster than native ones.

#### Testing Lanczos

Next test Lanczos in C with

```
./test_lanczos_c     
```

of in Fortran with

```
./test_lanczos_f
```
These should behave identically.  You will be asked for input parameters **n**, **nevec**, **maxstep**,
and **report**  (see above); maxstep should either be zero (no limit on the size of the Krylov subspace)
or large (it is just a protection against the process getting out of control, which should not happen).
To see the whole Lanczos process, set **report** = 1, or set it to 0  for no progress reports.
You will then be asked to choose between real-symmetric of complex-Hermition matrices, whether or not to
use BLAS, and if using BLAS in the complex Hermitian case, whether or not to "wrap" ZDOTC.

 * The test-program codes in directory "test_programs"  provide examples on how to invoke rs_lanczos and ch_lanczos in both Fortran and C.

=================================================================================

### Memory requirements

The implementation of Lanczos-Iteractor is in Fortran95 with a Fortran2003/C99 C
interface, with 64-bit "double" precision floating-point numbers, optionally using the BLAS
(Basic Linear Algebra Subprograms) for linear algebra of dimension-n vectors,
and LAPACK for finding the lowest eigenvector and eigenvalue of tridiagonal
Krylov matrices.   It is practical for very large sparse matrices, represented
by a list of non-zero matrix elements, so an efficient matvec is possible.

For very large n (millions), the main memory requirement is (n_v + 3)*n double
real or double complex floating-point numbers (Fortran only), with an additional
2*n when the C language interface is used (the 2*n are needed by the matvec
implementation v_out = matvec(n,v_in) in C).

### Distribution contents

```
Contents: (in subdirectory lib)

rs_lanczos.f90               simple driver subprogram for real-symmetric case
ch_lanczos.f90               simple driver subprogram for complex-Hermitian case
lanczos_iterator.f90         The basic Lanczos procedure
lanczos_m.f90                provides Hilbert-space tools for Lanczos
tmatrix_m.f90                provides Krylov-subspace (T-matrix) tools
vector_m.f90                 provides linear algebra operations in Hilbert space
matrix_algebra_m.f90         interface to LAPACK subprograms DSTEBZ and SPBTRD
zdotc_wrapper.f90            a wrapper for BLAS function zdotc
matrix_subspace_detector.f90 Can be used to check if the matrix has disconnected
                             subspaces in the MATVEC basis (in which case it 
                             could be put in block diagonal form by merely
                             reordering the basis), and produces a mask that 
                             can be used to project on individual subspaces.
Subdirectory C_interface

rs_lanczos_wrapper.f90       A Fortran2003 wrapper allowing rs_lanczos to be called from C
ch_lanczos_wrapper.f90       A Fortran2003  wrapper allowing ch_lanczos to be called from C
blas_settings_wrapper.f90    allows access to BLAS settings from C
README_C                     C Documentation

subdirectory test-programs

test_lanczos.f90             A test program in Fortran95
test_lanczos.c               A test program in C99
zdotc_test.f90               small program to test if your LAPACK/BLAS zdotc
                             is compatible with your compiler (if not, BLAS
                             double complex function zdotc needs a wrapper.)

subdirectory include

lanczos.h                    header file for C

Top level directory:

README                       Installation,use of  ch/rs_lanczos  drivers,C/C++ and Fortran
CHANGELOG                    Log of changes     
Implementation_Details       Full description of method
README_ZDOTC                 About problems with gfortran and BLAS subroutine zdotc 
setup.sh                     setup script

```                   
===============================================================================

