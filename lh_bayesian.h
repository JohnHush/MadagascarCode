/* This file is automatically generated. DO NOT EDIT! */

#ifndef _lh_bayesian_h
#define _lh_bayesian_h


void lh_exp_array( float *a	/* Arrary A */,
		   float *b	/* Array b=exp(a) */,
		   int n	/* Dimension */);
/*< Compute the exponential of the array a >*/


void lh_inv_array( float *a	/* Array A */,
		   float *b	/* Array b = 1/a */,
		   int n	/* DImension */);
/*< Compute the inverse elements in a >*/


void lh_vector_add_vector( float *a     /* Vector a */,
                           float *b     /* Vector b */,
                           float *c     /* Vector c=a+b */,
                           int n);
/*< c=a+b >*/


void lh_vector_sub_vector( float *a	/* Vector a */,
			   float *b	/* Vector b */,
			   float *c	/* Vector c=a-b */,
			   int n);
/*< c=a-b >*/


void lh_matrix_sub_matrix(float **a     /* Matrix A */,
                          float **b     /* Matrix B */,
                          float **c     /* Matrix C=A-B */,
                          int m         /* row of Matrix A&B&C */,
                          int n         /* column of the Matrixes */);
/*< C=A-B >*/


void lh_matrix_transpose( float **a	/* Matrix A */,
			  float **at	/* Matrix A' */,
			  int m		/* row of A */,
			  int n		/* column of A */);
/*< Transpose the matrix A to A' >*/


void lh_matrix_add_matrix(float **a	/* Matrix A */,
			  float **b	/* Matrix B */,
			  float **c	/* Matrix C=A+B */,
			  int m		/* row of Matrix A&B&C */,
			  int n		/* column of the Matrixes */);
/*< Matrix A + Matrix B= Matrix C, [m,n] >*/


void lh_matrix_mu_vector( float **a	/* Matrix A */,
			  float *b	/* Vector b */,
			  float *c	/* Vector c=Ab */,
			  int m		/* row of Matrix A */,
			  int n		/* column of Matrix A */);
/*< Using BLAS to compute matrix multiply vector, c= Ab >*/


void lh_matrix_mu_matrix( float **a	/* Matrix A */, 
			  float **b	/* Matrix B */,
			  float **c	/* Matrix C=AB */,
			  int m		/* row of Matrix A */,
			  int n		/* column of Matrix A & row of Matrix B */,
			  int k		/* column of Matrix B */);
/*< Using BLAS to compute matrix multiply matrix , C=AB>*/


void lh_matrix_inverse( float **a	/* Matrix A */,
			float **inv_a	/* Inverse of Matrix A */,
			int n		/* Dimension of the Matrix A */);
/*< Using Lapack Package to compute the inverse of the matrix A, LU factorization first! >*/


float lh_float_expectation( float *a	/* The array to be Computed */,
			    int n	/* The dimension of the Array*/);
/*< To compute the expectation of the array a[i] >*/


float lh_float_variance( float *a 	/* The array to be computed */,
			  int n		/* The dimension of the Array */);
/*< To compute the variance of the array a[i] >*/


float lh_float_covariance( float *x	/* The array x[i] */,
			   float *y	/* The array y[i] */,
			   int n	/* The dimension */);
/*< To compute the covariance of X and Y >*/


void lh_convolution_cut( float *x	/* Array X[i] */,
			 int nx		/* Dimension of X */,
			 float *wave	/* The wavelet array */,
			 int nwave	/* Length of wavelet */,
			 float *z	/* Output array */,
			 int phase 	/* The Time delay of wavelet */);
/*< Using routine lh_convolution plus a cut process to synthetic a convolved trace keep with the size of X >*/


void lh_convolution( float *x	/* Array X[i] */,
		     int nx	/* Dimension of X */,
		     float *y	/* Array Y[i] */,
		     int ny	/* Dimension of Y */,
		     float *z	/* Output array */);
/*< Convolution of Two array, Z = cov<x,y> >*/

#endif
