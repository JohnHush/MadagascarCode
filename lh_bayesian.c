#include <rsf.h>
#include "lh_bayesian.h"

void lh_exp_array( float *a	/* Arrary A */,
		   float *b	/* Array b=exp(a) */,
		   int n	/* Dimension */)
/*< Compute the exponential of the array a >*/
{
	int i;
	for( i=0 ; i<n ; i++ )
		b[i] = exp( a[i]);
}

void lh_inv_array( float *a	/* Array A */,
		   float *b	/* Array b = 1/a */,
		   int n	/* DImension */)
/*< Compute the inverse elements in a >*/
{
	int i;
	for( i=0 ; i<n ; i++ )
	{
		if( a[i]==0 )
		{
			printf( "lh_inv_array, value too mmall!!\n" );
			exit(0);
		}
		b[i] = 1./a[i];
	}
}

void lh_vector_add_vector( float *a     /* Vector a */,
                           float *b     /* Vector b */,
                           float *c     /* Vector c=a+b */,
                           int n)
/*< c=a+b >*/
{
        int i;
        for( i=0 ; i<n ; i++ )
                c[i] = a[i]+b[i];
}

void lh_vector_sub_vector( float *a	/* Vector a */,
			   float *b	/* Vector b */,
			   float *c	/* Vector c=a-b */,
			   int n)
/*< c=a-b >*/
{
	int i;
	for( i=0 ; i<n ; i++ )
		c[i] = a[i]-b[i];
}

void lh_matrix_sub_matrix(float **a     /* Matrix A */,
                          float **b     /* Matrix B */,
                          float **c     /* Matrix C=A-B */,
                          int m         /* row of Matrix A&B&C */,
                          int n         /* column of the Matrixes */)
/*< C=A-B >*/
{
	int i,j;
	for( i=0 ; i<m ; i++ )
	for( j=0 ; j<n ; j++ )
		c[i][j] = a[i][j]-b[i][j];
}

void lh_matrix_transpose( float **a	/* Matrix A */,
			  float **at	/* Matrix A' */,
			  int m		/* row of A */,
			  int n		/* column of A */)
/*< Transpose the matrix A to A' >*/
{
	int i,j;
	for( i=0 ; i<m ; i++ )
	for( j=0 ; j<n ; j++ )
	{
		at[j][i] = a[i][j];
	}
}

void lh_matrix_add_matrix(float **a	/* Matrix A */,
			  float **b	/* Matrix B */,
			  float **c	/* Matrix C=A+B */,
			  int m		/* row of Matrix A&B&C */,
			  int n		/* column of the Matrixes */)
/*< Matrix A + Matrix B= Matrix C, [m,n] >*/
{
	int i,j;
	for( i=0 ; i<m ; i++ )
	for( j=0 ; j<n ; j++ )
		c[i][j] = a[i][j]+b[i][j];
}

void lh_matrix_mu_vector( float **a	/* Matrix A */,
			  float *b	/* Vector b */,
			  float *c	/* Vector c=Ab */,
			  int m		/* row of Matrix A */,
			  int n		/* column of Matrix A */)
/*< Using BLAS to compute matrix multiply vector, c= Ab >*/
{
	float alpha=1., beta=0.;
	char ctr = 'N';
	int i,j, inc=1;
	float *blas_a , *blas_b , *blas_c;

	blas_a = sf_floatalloc( m*n );
	blas_b = sf_floatalloc( n );
	blas_c = sf_floatalloc( m );

	for( i=0 ; i<m ; i++ )
	for( j=0 ; j<n ; j++ )
		blas_a[i+j*m] = a[i][j];

	for( i=0 ; i<n ; i++ )
		blas_b[i] = b[i];
	sgemv_( &ctr , &m , &n , &alpha, blas_a, &m, blas_b, &inc, &beta, blas_c, &inc );
	for( i=0 ; i<m ; i++ )
		c[i] = blas_c[i];

	free1float(blas_a);
	free1float(blas_b);
	free1float(blas_c);
}

void lh_matrix_mu_matrix( float **a	/* Matrix A */, 
			  float **b	/* Matrix B */,
			  float **c	/* Matrix C=AB */,
			  int m		/* row of Matrix A */,
			  int n		/* column of Matrix A & row of Matrix B */,
			  int k		/* column of Matrix B */)
/*< Using BLAS to compute matrix multiply matrix , C=AB>*/
{
	float alpha=1., beta=0;
	char ctr = 'N';
	int i,j;
	float *blas_a, *blas_b , *blas_c;
	
	blas_a = sf_floatalloc( m*n );
	blas_b = sf_floatalloc( n*k );
	blas_c = sf_floatalloc( m*k );

	for( i=0 ; i<m ; i++ )
	for( j=0 ; j<n ; j++ )
		blas_a[i+j*m] = a[i][j];

	for( i=0 ; i<n ; i++ )
	for( j=0 ; j<k ; j++ )
		blas_b[i+j*n] = b[i][j];

	sgemm_( &ctr , &ctr , &m , &k , &n , &alpha , blas_a , &m , blas_b , &n , &beta , blas_c , &m );

	for( i=0 ; i<m ; i++ )
	for( j=0 ; j<k ; j++ )
		c[i][j] = blas_c[i+j*m];

	free1float(blas_a);
	free1float(blas_b);
	free1float(blas_c);
}

void lh_matrix_inverse( float **a	/* Matrix A */,
			float **inv_a	/* Inverse of Matrix A */,
			int n		/* Dimension of the Matrix A */)
/*< Using Lapack Package to compute the inverse of the matrix A, LU factorization first! >*/
{
	float *lapack_a , *WORK;
	int i, j, LWORK=n, info, *ipiv;

	ipiv     = sf_intalloc ( n );
	WORK     = sf_floatalloc( LWORK );
	lapack_a = sf_floatalloc( n*n );
	
	for( i=0 ; i<n ; i++ )
	for( j=0 ; j<n ; j++ )
		lapack_a[i+j*n] = a[i][j];

	sgetrf_( &n , &n , lapack_a , &n , ipiv , &info );
	sgetri_( &n , lapack_a , &n , ipiv , WORK, &LWORK , &info );

	for( i=0 ; i<n ; i++ )
	for( j=0 ; j<n ; j++ )
		inv_a[i][j] = lapack_a[i+j*n];

	free1float(ipiv);
	free1float(WORK);
	free1float(lapack_a);
}

float lh_float_expectation( float *a	/* The array to be Computed */,
			    int n	/* The dimension of the Array*/)
/*< To compute the expectation of the array a[i] >*/
{
	float sum=0.;
	int i;
	for( i=0 ; i<n ; i++ )
		sum+=a[i]/n;

	return sum;
}

float lh_float_variance( float *a 	/* The array to be computed */,
			  int n		/* The dimension of the Array */)
/*< To compute the variance of the array a[i] >*/
{
	float exp_a;
	float sum=0.;
	int i;

	exp_a = lh_float_expectation( a , n );

	for( i=0 ; i<n ; i++ )
		sum+=(a[i]-exp_a)*(a[i]-exp_a)/n;

	return sum;
}

float lh_float_covariance( float *x	/* The array x[i] */,
			   float *y	/* The array y[i] */,
			   int n	/* The dimension */)
/*< To compute the covariance of X and Y >*/
{
	float exp_x, exp_y;
	float sum=0.;
	int i;

	exp_x = lh_float_expectation( x , n );
	exp_y = lh_float_expectation( y , n );

	for( i=0 ; i<n ; i++ )
		sum+=(x[i]-exp_x)*(y[i]-exp_y)/n;

	return sum;
}

void lh_convolution_cut( float *x	/* Array X[i] */,
			 int nx		/* Dimension of X */,
			 float *wave	/* The wavelet array */,
			 int nwave	/* Length of wavelet */,
			 float *z	/* Output array */,
			 int phase 	/* The Time delay of wavelet */)
/*< Using routine lh_convolution plus a cut process to synthetic a convolved trace keep with the size of X >*/
{
	float tmp[nx+nwave-1];
        int i;
        lh_convolution( x , nx , wave , nwave , tmp);
        for( i=0 ; i<nx ; i++ )
                z[i] = tmp[i+phase];
}

void lh_convolution( float *x	/* Array X[i] */,
		     int nx	/* Dimension of X */,
		     float *y	/* Array Y[i] */,
		     int ny	/* Dimension of Y */,
		     float *z	/* Output array */)
/*< Convolution of Two array, Z = cov<x,y> >*/
{
        int i,j;
        float tmp1[nx];
        float tmp2[ny];

        for(i=0;i<nx+ny-1;i++)
        {
                z[i]=0.0;
        }

        if(nx>ny)
        {
                for(i=0;i<ny;i++)
                {
                        tmp2[i]=y[ny-1-i];
                }
                for(i=0;i<ny;i++)
                {
                        for(j=0;j<=i;j++)
                        {
                                z[i]=z[i]+x[j]*tmp2[ny-1-i+j];
                        }
                }
                for(i=ny;i<nx;i++)
                {
                        for(j=0;j<ny;j++)
                        {
                                z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
                        }
                }
                for(i=nx;i<nx+ny-1;i++)
                {
                        for(j=0;j<ny-i+nx-1;j++)
                        {
                                z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
                        }
                }
        }
        else if(nx==ny)
        {
                for(i=0;i<ny;i++)
                {
                        tmp2[i]=y[ny-1-i];
                }
                for(i=0;i<ny;i++)
                {
                        for(j=0;j<=i;j++)
                        {
                                z[i]=z[i]+x[j]*tmp2[ny-1-i+j];
                        }
                }
                for(i=ny;i<2*ny-1;i++)
                {
                        for(j=0;j<2*ny-i-1;j++)
                        {
                                z[i]=z[i]+tmp2[j]*x[i-ny+1+j];
                        }
                }
        }
        else
        {
                for(i=0;i<nx;i++)
                {
                        tmp1[i]=x[nx-1-i];
                }
                for(i=0;i<nx;i++)
                {
                        for(j=0;j<=i;j++)
                        {
                                z[i]=z[i]+y[j]*tmp1[nx-1-i+j];
                        }
                }
                for(i=nx;i<ny;i++)
                {
                        for(j=0;j<nx;j++)
                        {
                                z[i]=z[i]+tmp1[j]*y[i-nx+1+j];
                        }
                }
                for(i=ny;i<nx+ny-1;i++)
                {
                        for(j=0;j<nx-i+ny-1;j++)
                        {
                                z[i]=z[i]+tmp1[j]*y[i-nx+1+j];
                        }
                }
        }
}


