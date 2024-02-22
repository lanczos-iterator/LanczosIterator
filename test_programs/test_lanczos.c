#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <complex.h>
#include "lanczos.h"

/* simple sample ch_matvec_c for testing */
void my_ch_matvec(int n, double complex *vec) {
  for (int i = 0; i < n; i++) {
    vec[i] =  vec[i] * i;
  }
}

/* simple sample rs_matvec_c for testing */
void my_rs_matvec(int n, double  *vec) {
  for (int i = 0; i < n; i++) {
    vec[i] =  vec[i] * i;
  }
}

int main() {
  int n, nevec, maxstep, report, seed;
  int x, n1, matrix_type;
  double *eval, *variance, resolution;
  double *evec_rs = NULL;
  double complex *evec_ch = NULL;
  bool use_blas, wrap_zdotc;
  bool nocheck = false;
  int *range = (int *) malloc(2 * sizeof(int));
    
  printf("Lanczos test program in C\n");
  printf("  n = dimension of matrix\n");
  printf("  nevec = number of eigenvalues requested\n");
  printf("  maxstep is number of Lanczos steps allowed (0 for no limit)\n");
  printf("  give report = r > 0 to see progress every r Lanczos steps\n"
	 "  (plus summary of Lanczos process for each eigenvector)\n\n");
  printf("give n, nevec, maxstep, report\n");
  scanf("%i %i %i %i", &n, &nevec, &maxstep, &report);
  seed = 12345678;

  eval = (double *) malloc(n * sizeof(double));
  variance = (double *) malloc(n * sizeof(double));

  printf(" This program can use BLAS = Basic Linear Algebra Suprograms\n");
  printf("enter 1 to use BLAS, 0 for no BLAS\n");
  scanf("%i", &x);
  switch (x) {
  case 0:
  case 1:
    use_blas = (bool) x;
    set_use_blas(use_blas);
    break;
  default:
    printf("invalid choice\n");
    return -1;
  }
  
  printf("\nThe test matrix is  diagonal with eigenvalues (0,1,2,3,....)\n");
  printf("(The Lanczos algorithm is not aware that the matrix is diagonal\n\n");
 
  printf("give 1 for real symmetric, 2 for complex Hermitian\n");
  scanf("%i", &matrix_type);

  switch (matrix_type) {
  case 1:
    evec_rs = (double *) malloc( n * nevec * sizeof(double));
    /* examples of the two computational modes: */
    /* compute the first n1 eigenvectors in all-in-one-go mode (*range = NULL) */
    /* when *range = NULL, value of nocheck does not matter */
    n1 = 1 + nevec/2;
    seed = (17*seed)%7654321;
    rs_lanczos_c(n, n1,  eval,  evec_rs,  variance,
		 &resolution, &my_rs_matvec, maxstep, report, seed, NULL, 0);

    /* compute the remaining eigenvectors in one-by-one mode (range[0]=range[1])*/
    /* when (range[0],range[1]) is specified,the validity of the previously-found
     * eigenvectors will by default (nocheck= false)  be checked by computing their
     *  eigenvalue and variance and comparing with the values in eval and variance. */

    nocheck = false;   /* check that the previous eigenvector data is correct*/
    for (int i = n1 + 1; i <= nevec; i++) {
      seed = (17*seed)%7654321;
      range[0] = i; range[1] = i;  /* compute only eigenvalue i */
      rs_lanczos_c(n, i,  eval,  evec_rs,  variance,
		   &resolution, &my_rs_matvec, maxstep, report, seed, range, nocheck);
      nocheck = true;  /* no point in rechecking again  after the restart*/       
    }
    break;
  case 2:
    if (use_blas) {
      printf("enter 0 to use BLAS without a ZDOTC wrapper, 1 for wrapper\n");
      scanf("%i", &x);
      switch (x) { 
      case 0:
      case 1:
	wrap_zdotc = (bool) x;
	set_wrap_zdotc(wrap_zdotc);
	break;
      default:
	printf("invalid choice\n");
	return -1;
      }
    } 
    evec_ch = (double complex *) malloc( n * nevec * sizeof(double complex));
    /* examples of the two computational modes: */
    /* compute the first n1 eigenvectors in all-in-one-go mode  (*range = NULL) */
    /*when *range = NULL, value of nocheck does not matter */
    n1 = 1 + nevec/2;
    seed = (17*seed)%7654321;
    ch_lanczos_c(n, n1,  eval,  evec_ch,  variance,
		 &resolution, &my_ch_matvec, maxstep, report, seed, NULL, 0);
    
    /* compute the remaining eigenvectors in one-by-one mode (range[0]=range[1])*/
    /* when (range[0],range[1]) is specified,the validity of the previously-found
     * eigenvectors will by default (nocheck= false)  be checked by computing their
     *  eigenvalue and variance and comparing with the values in eval and variance. */
    
    nocheck = false;   /* check that the previous eigenvector data is correct*/
    for (int i = n1 + 1; i <= nevec; i++) {
      seed = (17*seed)%7654321;
      range[0] = i; range[1] = i;  /* compute only eigenvalue i */
      ch_lanczos_c(n, i,  eval,  evec_ch,  variance,
		   &resolution, &my_ch_matvec, maxstep, report, seed, range, nocheck);
      nocheck = true;  /* no point in rechecking again  after the restart*/
    }
    break;
    default:
      printf("invalid choice\n");
      return -1;
  }
  for (int i = 0; i < nevec; i++) {
    printf("%10d %25.16f variance = %10.3e\n", i+1, eval[i], variance[i]);
  }
  printf("  resolution: %10.3e\n",resolution);
  printf("  eigenvalue range: [%d:%d]\n", range[0], range[1]);
  if (evec_rs) free(evec_rs);
  if (evec_ch) free(evec_ch);
  
  free (eval);
  free (variance);
  return 0;
}
  
  
