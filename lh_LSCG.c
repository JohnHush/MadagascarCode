#include <rsf.h>
#include "lh_LSCG.h"
#include "lh_bayesian.h"
#include "lh_mp.h"

static void matrix_mu_vector( float **a , float *b , int nz , int nx , float *c );
static void Tmatrix_mu_vector( float **a , float *b , int nz , int nx , float *c );

static void matrix_mu_vector( float **a , float *b , int nz , int nx , float *c )
{
        int i,j;
        for( i=0 ; i<nz ; i++ )
        {
                c[i]=0.0;
                for( j=0 ; j<nx ; j++ )
                {
                        c[i]+=a[i][j]*b[j];
                }
        }
}

static void Tmatrix_mu_vector( float **a , float *b , int nz , int nx , float *c )
{
        int i,j;
        for( j=0 ; j<nx ; j++ )
        {
                c[j]=0.0;
                for( i=0 ; i<nz ; i++ )
                {
                        c[j]+=a[i][j]*b[i];
                }
        }
}

void lh_direct_LS( float **G	/* The Forward Operator */,
		   int Nr	/* Row Of G */,
		   int Nc	/* Column of G */,
		   float **Cd	/* The Covariance Matrix of Data Cd[Nr][Nr] */,
		   float **Cm 	/* The Covariance Matrix of model Cm[Nc][Nc] */,
		   float *d	/* Data */,
		   float *m_prior /* m prior */,
		   float *m_post  /* m posterior */ )
/*< Solve the Least Square Optimization Problem directly Using Matrix Inverse >*/
{
    float **Gt, **A, **B1, **B2, **B3, **B, *C1, *C, **AB;

    Gt = alloc2float( Nr , Nc );
    A  = alloc2float( Nr , Nc );
    B1 = alloc2float( Nc , Nr );
    B2 = alloc2float( Nr , Nr );
    B3 = alloc2float( Nr , Nr );
    B  = alloc2float( Nr , Nr );
    C1 = alloc1float( Nr );
    C  = alloc1float( Nr );
    AB = alloc2float( Nr , Nc );

    zero2float( Gt , Nr , Nc );
    zero2float( A  , Nr , Nc );
    zero2float( B1 , Nc , Nr );
    zero2float( B2 , Nr , Nr );
    zero2float( B3 , Nr , Nr );
    zero2float( B  , Nr , Nr );
    zero1float( C1 , Nr );
    zero1float( C  , Nr );
    zero2float( AB , Nr , Nc );

    lh_matrix_transpose( G , Gt , Nr , Nc );
    lh_matrix_mu_matrix( Cm , Gt , A , Nc , Nc , Nr );
    lh_matrix_mu_matrix( G , Cm , B1 , Nr , Nc , Nc );
    lh_matrix_mu_matrix( B1 , Gt , B2 , Nr , Nc , Nr );
    lh_matrix_add_matrix( B2 , Cd , B3 , Nr , Nr );
    lh_matrix_inverse( B3 , B , Nr );
    lh_matrix_mu_vector( G , m_prior , C1 , Nr , Nc );
    lh_vector_sub_vector( d , C1 , C , Nr );
    lh_matrix_mu_matrix( A , B , AB , Nc , Nr , Nr );
    lh_matrix_mu_vector( AB , C , m_post , Nc , Nr );
    lh_vector_add_vector( m_prior , m_post , m_post , Nc );

    free2float( Gt );
    free2float( A  );
    free2float( B1 );
    free2float( B2 );
    free2float( B3 );
    free2float( B  );
    free1float( C1 );
    free1float( C );
    free2float( AB );
}

void LSCG_revised( float **a		/* The Matrix A in the left Hand A[NZ][NX] NZ row number NX column number */,
		   float *x		/* The solution vector */,
		   float *b		/* The data vector in the right Hand */,
		   int nz		/* The Row number of A */,
		   int nx		/* The Column Number of A */,
		   float *lamda		/* The Pre-whiting Coefficient */,
		   float eps		/* The tolerance error */,
		   int MaxStep		/* The Maximum Iterative step*/)
/*< LSCG Method To Solve Ax=b >*/
{
        float *h,*g,*p,*ax,*ath,*ap,b_sum=0.0,alpha,beta;
        float sum_g,sum_ap,sum_p,sum_g_new,sum_error, *deviation;
        int i,istep;

        h=malloc( nz*sizeof( float ) );
        g=malloc( nx*sizeof( float ) );
        p=malloc( nx*sizeof( float ) );
        ax=malloc( nz*sizeof( float ) );
        ath=malloc( nx*sizeof( float ) );
        ap=malloc( nz*sizeof( float ) );
        deviation=malloc( ( MaxStep+10 )*sizeof( float ) );
// ************STEP1****COMPUTE the Energy of the Input Solution ************************ //
        for( i=0;i<nz;i++ )
                b_sum+=b[i]*b[i];
        b_sum=sqrt(b_sum);
// ************STEP1i********* END******************************************************* //
// ************STEP2******ZERO the Arrays *********************************************** //
        for( i=0 ; i<nx ; i++ )
        {
                g[i]=0.0;
                p[i]=0.0;
                ath[i]=0.0;
        }
        for( i=0 ; i<nz ; i++ )
        {
                h[i]=0.0;
                ax[i]=0.0;
                ap[i]=0.0;
        }
        for( i=0 ; i<MaxStep+10 ; i++ )
        {
                deviation[i]=0.0;
        }
// ************STEP2******END*********************************************************** //
// ************STEP3******Compute Initial Error **************************************** //
        matrix_mu_vector( a , x , nz , nx , ax );
        for( i=0;i<nz;i++ )
        {
                h[i]=b[i]-ax[i];
        }
        Tmatrix_mu_vector( a , h , nz , nx , ath );
        for( i=0;i<nx;i++ )
        {
                g[i]=ath[i]-lamda[i]*x[i];
                p[i]=g[i];
        }
// ***********STEP3******END********************************************************** //
// ***********STEP4********Main Iteration********************************************* //

        for( istep=0 ; istep<MaxStep ; istep++ )
        {
                //printf("LSCG_step=%d\n",istep);
        // **********STEP4-1******Compute the length of the Iteration***************** //
                sum_error=0.0;sum_g=0.0;sum_ap=0.0;sum_p=0.0;sum_g_new=0.0;alpha=0.0;beta=0.0;
                matrix_mu_vector( a , p , nz , nx , ap );
                for( i=0 ; i<nx ; i++ )
                {
                        sum_g+=g[i]*g[i];
                        sum_p+=p[i]*lamda[i]*p[i];
                }
                for( i=0 ; i<nz ; i++ )
                {
                        sum_ap+=ap[i]*ap[i];
                }
                alpha=sum_g/(sum_ap+sum_p);
        // **********STEP4-1******END************************************************ //
        // **********STEP4-2***Compute the New solution and the Error**************** //
                for( i=0 ; i<nx ; i++ )
                {
                        x[i]=x[i]+alpha*p[i];
                }
                for( i=0 ; i<nz ; i++ )
                {
                        h[i]=h[i]-alpha*ap[i];
                }
                Tmatrix_mu_vector( a , h , nz , nx , ath );
        // **********STEP4-2********END********************************************* //
        // **********STEP4-3*****COmpute the direction and the length of the Iter*** //
                for( i=0 ; i<nx ; i++ )
                {
                        g[i]=ath[i]-lamda[i]*x[i];
                }
                for( i=0 ; i<nx ; i++ )
                {
                        sum_g_new+=g[i]*g[i];
                }
                beta=sum_g_new/sum_g;
                for( i=0 ; i<nx ; i++ )
                {
                        p[i]=g[i]+beta*p[i];
                }
        // **********STEP4-3*********END******************************************* //
        // **********STEP4-4****BReak Out when the error is acceptable************* //

                for( i=0 ; i<nz ; i++ )
                {
                        sum_error+=h[i]*h[i];
                }
                sum_error=sqrt( sum_error );
                if( sum_error<eps*b_sum )
                {
                        printf("the error is small enough\n");
                        break;
                }
                deviation[istep]=sum_error;
        /*      if( istep>=3 )
                {
                        if( deviation[istep]>(deviation[istep-1]+deviation[istep-2]+deviation[istep-3])/3.0 )
                        {
                                printf("Break From the Iteration while the Functional of Error is growing\n");
                                printf("The Number of Iteration = %d\n" , istep );
                                break;
                        }
                }
*/
                // *********STEP4-4*********END******************************************* //

        }
// *******STEP4*********END******************************************************* //
        free( h );free( g );free( p );free( ax );free( ath );free( ap ),free(deviation);

}

void lh_CG( float **a	/* The Matrix A[N][N] should be positive definite */,
	   float *b	/* THe Right-hand vecotr of the Equation Ax=b */,
	   float *x	/* The Solution */,
	   int n	/* The dimension */,
	   float eps	/* The error torrelence */)
/*< Solve Equation Ax=b, A is positive definite matrix ... Using Conjugate Gradient Method >*/
{
	float *r;
	float *p;
	float *ax;
	float *ap;
	int i,j;

	float alpha,beta,sum1,sum2;

	r=malloc( n*sizeof( float ) );
	p=malloc( n*sizeof( float ) );
	ax=malloc( n*sizeof( float ) );
	ap=malloc( n*sizeof( float ) );

	for( i=0;i<n;i++ )
	{
		r[i]=0.0;
		p[i]=0.0;
		ap[i]=0.0;
		ax[i]=0.0;
	}

	matrix_mu_vector( a , x , n , n , ax );

	for( i=0;i<n;i++ )
	{
		r[i]=b[i]-ax[i];
		p[i]=b[i]-ax[i];
        }

        for( i=0;i<n;i++ )
        {
		sum1 = lh_float_iproduct( r , p , n );
		matrix_mu_vector( a , p , n , n , ap );
		sum2 = lh_float_iproduct( ap , p , n );

		alpha=sum1/sum2;

		for( j=0;j<n;j++ )
			x[j]=x[j]+alpha*p[j];

		matrix_mu_vector( a , x , n , n , ax );

		lh_vector_sub_vector( b , ax , r , n );

		sum1 = lh_float_iproduct( r , ap , n );

		beta=-sum1/sum2;

		if( lh_float_iproduct( r , r , n )/lh_float_iproduct( ax , ax , n )<eps )
			break;

		for( j=0;j<n;j++ )
			p[j]=r[j]+beta*p[j];
        }
	free( p );
	free( r );
	free( ap );
	free( ax );
}
