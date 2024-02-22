#ifndef LANCZOS_H
#define LANCZOS_H
#ifdef __cplusplus
extern "C" {
#endif
  
/* matvec must be supplied by the user.  For the real symmetric case (rs_lanczos_c),
   its C code  should be equivalent to 
   
   double *matrix;  //a pointer to the full matrix       
   void matvec(int n, double *vec) {
     double *vec_out = (double *) calloc(n, sizeof(double));
     for (int i = 0; i < n; i++) {
       for (int j = 0; j < n; j++) {
	 vec_out[i] += matrix[i][j] * vec[j];
       }
     }
   memcpy(vec, vec_out, n * sizeof(double));
   free (vec_out);
   }
   
   For the complex Hermitian case (ch_lanczos_c), replace "double" by "double _Complex".
   
   Preferably, an equivalent efficient implementation of matvec for sparse matrices
   should be used in matvec.
   
*/
  void rs_lanczos_c(int n, int nevec, double *eval, double *evec, double *variance,
		    double *resolution,
		    void (*matvec) (int n, double *v),
		    int maxstep,  int report, int seed, int *range, _Bool nocheck);
  
  void ch_lanczos_c(int n, int nevec, double *eval, double _Complex *evec,
		    double *variance,double *resolution,
		    void (*matvec) (int n, double _Complex *v),
		    int maxstep,  int report, int seed, int *range, _Bool nocheck);
  
  void set_use_blas(_Bool use_blas);
  void set_wrap_zdotc(_Bool wrap_zdotc);


  double *create_real_storage_c(int n);
  void destroy_real_storage_c(double *ptr);
  double _Complex *create_complex_storage_c(int n);
  void destroy_complex_storage_c(double _Complex *ptr);
  
  /* executable code for these functions is included below: */

  double *create_real_storage_c(int n) {
    return (double *) malloc(n * sizeof(double));
  }
  void destroy_real_storage_c(double *ptr) {
    free (ptr);
  }
  double _Complex *create_complex_storage_c(int n) {
    return (double _Complex *) malloc(n * sizeof(double _Complex));
  }
  void destroy_complex_storage_c(double _Complex *ptr) {
    free (ptr);
  }  

  
  /* BLAS (Basic Linear Algebra Subprograms) and BLAS zdotc implementation issues:
     
     Vendor-optimized BLAS (e.g. Intel MKL, Apple framework Accelerate, NVIDIA NVPL)
     is preferable for carrying out vector algebra on very large vectors.
     However there are two conventions for complex return values from  double complex functions
     such as ZDOTC, and the optimized implementations may not be compatible with your
     FORTRAN compiler (in particular gfortran) unless a wrapper is used.

     To use the wrapper, call "set_wrap_zdotc(1)"
     To restore the default (no wrapper) call "set_wrap_zdotc(0)"
     
     Alternatively, BLAS use can be "switched off" by calling "set_use_blas(0)"
     Native Fortran95 vector algebra will then be used instead of BLAS.
     (To restore BLAS use: call "set_use_blas(1)".  
  */  

#endif   
#ifdef __cplusplus
}
#endif

