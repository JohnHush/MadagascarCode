/* AVA Inversion In Elastic Media */

/*********************************************************************************
* Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
*
* File Name: Mava_inversion_casestudy.c
*
* Program Context: Bayesian Inversion scheme in ELASTIC Media
*
* Author: Long Teng( Originally designed in 2013 summer).
*
* Revisor:  Heng Luo ( Revised to C code in 2014 Spring).
*
* Version: 1.1.0
*
* Date: 2014/03/26
*
* History:
*         version 1.0.0: 2014/02/13	( Revised by Heng Luo )
*         version 1.0.1: 2014/02/18	( Revised by Heng Luo )
*         version 1.0.2: 2014/02/26	( Revised by Heng Luo )
*         version 1.0.3: 2014/03/06	( Revised by Heng Luo )
*         version 1.0.4: 2014/03/10	( Revised by Long Teng & Heng Luo )
*		    Add Wavelet and welllog interface & changeable parameters
*	      version 1.0.5: 2014/03/18     ( Revised by Long Teng)
*           1.Add normization programme 
*		    2.Amend the D matrix to (nt*npara x nt*npara)size
*			3.Output inline NO.445(near the well)	for invserion QC
*         version 1.0.6: 2014/03/26	( Revised by Heng Luo )
*         version 1.0.7: 2014/09/30	( Revised by Long TENG )
*         version 1.1.0: 2014/11/06 ( Revised by Heng Luo, Totally New Version )
* 
* Reference.: Buland and Omre, 2003, Baysian Linearized AVO inversion
*             Innanen, 2011, Anelastic AVF/AVA inversion
*
* Equation1: R(sita) = 0.5*delta_rou/rou  + 0.5*(1+tan(sita)*tan(sita))*delta_v/v;
* Equation2: Aki-Richards Equation
*
**********************************************************************************/

#include <rsf.h>
#include "su_alloc.h"
#include "lh_readwrite.h"
#include "lh_wavelet.h"
#include "lh_bayesian.h"
#include "lh_LSCG.h" 

int main( int argc , char* argv[] )
{
	sf_init( argc , argv );

    sf_axis at, ax, ay, at_wave;
    /*
     * @ Setup Axises for DATASETS
     */
    sf_file FD1, FD2, FD3, FD4, FW1, FW2, FW3, FW4, FPR1, FPR2, FPR3, FPO1, FPO2, FPO3, FTR1, FTR2, FTR3;
    /*
     * @ Setup FILE INTERFACES for DATASETS
     */
    int ix, nx, nx_s, NX, iy, ny, ny_s, NY, it, nt, nt_s, NT, nwave, jt;
    /*
     * DIMENSTION PARAMETERS, IN TIME, INLINE & CROSSLINE, RESPECTIVELY
     */
    float sita1, sita2, sita3, sita4;
    /*
     * @ Angle Values of NearOffset, Mid1, Mid2 and TeleOffset, res.
     */
    float crange, corder, var_e;
    /*
     * @ Bayesian Inversion Parameters
     */
    int wsft;
    /*
     * @ Wavelet Shift
     */
    int dim, num_dim[SF_MAX_DIM];
    /*
     * @ Dimensinon detecting
     */
    float dt, ot;
    /*
     * @ TIME DIMENSINO PARAMETERS
     */
    float *NearWave, *Mid1Wave, *Mid2Wave, *TeleWave;
    /*
     * Four Input Wavelets
     */
    float *d, *M_prior, **Cm, **Cn, **A, **W, **D, **WA, **G, *M_post, ****MPOST, *M_true;
    /*
     * @ d=DATA, M_prior= A Prior Model,
     * @ Cm is the covaraince matrix of the model
     * @ Cn is the covariance matrix of the data
     * @ G=WAD
     */
    int ipara, jpara, npara=3, nsita=4;
    float vpvsratio;
    /*
     * @ auxiliary variables
     */

    if(!sf_getint( "ninline_start" , &nx_s ))        sf_error( "Missing inline Direction Starting Value !\n" );
    /* Input Parameter: the Starting value in inline direction in Inversion */
    if(!sf_getint( "ninline_inv" , &nx ))            sf_error( "Missing Number of CDPS in inline direction!\n" );
    /* Input Parameter: the number of CDPS in inline direction */
    if(!sf_getint( "ncrossline_start" , &ny_s ))     sf_error( "Missing crossline Direction Starting Value!\n" );
    /* Input Paramter: the starting Value in crossline direction in Inversion */
    if(!sf_getint( "ncrossline_inv" , &ny ))         sf_error( "Missing Number of inlines in Crossline direction!\n" );
    /* Input Parameter: the number of INLINES in crossline direction */
    if(!sf_getfloat( "nearsita" , &sita1 ))          sf_error( "Missing the Near Offset angle Value!\n" );
    /* Input Parameter: the Angle Value of Near Offset Trace Gather */
    if(!sf_getfloat( "mid1sita" , &sita2 ))          sf_error( "Missing the Mid1 Offset angle Value!\n" );
    /* Input Parameter: the Angle Value of Mid1 Offset Trace Gather */
    if(!sf_getfloat( "mid2sita" , &sita3 ))          sf_error( "Missing the Mid2 Offset angle Value!\n" );
    /* Input Parameter: the Angle Value of Mid2 Offset Trace Gather */
    if(!sf_getfloat( "telesita" , &sita4 ))          sf_error( "Missing the Tele Offset angle Value!\n" );
    /* Input Parameter: the Angle Value of Tele Offset Trace Gather */
    if(!sf_getint( "wave_shift" , &wsft ))           sf_error( "Missing the Shift of Wavelet!\n" );
    /* Input Parameter: Wave Shift in Inversion */
    if(!sf_getfloat( "correlation_range" , &crange ))sf_error( "Missing the Correlation Range in Inversion!\n" );
    /* Input Parameter: COrrelation range in Inversion */
    if(!sf_getfloat( "correlation_order" , &corder ))sf_error( "Missing the Correlation Order in Inversion!\n" );
    /* Input Parameter: Correlation Order in Inversion */
    if(!sf_getfloat( "variance_noise" , &var_e ))    sf_error( "Missing the Variance of Noise in Inversion!\n" );
    /* Input Parameter: Variance of Noise in Inversion */


    FD1 = sf_input( "NearOffsetGather" );
    FD2 = sf_input( "Mid1OffsetGather" );
    FD3 = sf_input( "Mid2OffsetGather" );
    FD4 = sf_input( "TeleOffsetGather" );
    /*
     * @ Input DATASETS
     */
    FW1 = sf_input( "NearWavelet" );
    FW2 = sf_input( "Mid1Wavelet" );
    FW3 = sf_input( "Mid2Wavelet" );
    FW4 = sf_input( "TeleWavelet" );
    /*
     * @ Input WAVELETs
     */
    FPR1= sf_input( "vp_prior" );
    FPR2= sf_input( "vs_prior" );
    FPR3= sf_input( "ro_prior" );

    FTR1= sf_input( "vp_true" );
    FTR2= sf_input( "vs_true" );
    FTR3= sf_input( "ro_true" );

    FPO1= sf_output( "vp_post" );
    FPO2= sf_output( "vs_post" );
    FPO3= sf_output( "ro_post" );
    /*
     * @ A Prior Models, True MOdels and Posterior Models
     */
	
	dim = sf_filedims ( FD1 , num_dim );

	if( dim==1 )
	{
		at = sf_iaxa( FD1 , 1 );
		NT = sf_n( at );
		dt = sf_d( at );
		ot = sf_o( at );

		NX = 1;
		NY = 1;
	}
	else if( dim==2 )
	{
		at = sf_iaxa( FD1 , 1 );
		NT = sf_n( at );
		dt = sf_d( at );
		ot = sf_o( at );

		ax = sf_iaxa( FD1 , 2 );
		NX = sf_n( ax );
		
		NY = 1;
	}
	else if( dim==3 )
	{
		at = sf_iaxa( FD1 , 1 );
		NT = sf_n( at );
		dt = sf_d( at );
		ot = sf_o( at );

		ax = sf_iaxa( FD1 , 2 );
		NX = sf_n( ax );
	
        ay = sf_iaxa( FD1 , 3 );    
		NY = sf_n( ay );

	}
	else
	{
		sf_error( "Can't handle the problem with dimension higher than 3!\n " );
	}

	at_wave = sf_iaxa( FW1 , 1 );
	nwave   = sf_n( at_wave );

    if( sf_n(sf_iaxa( FPR1 , 1 ))!=NT+1 )   sf_error( "DATASETS UNMATCHED!\n" );
    if( sf_n(sf_iaxa( FPR1 , 2 ))!=NX )     sf_error( "DATASETS UNMATCHED!\n" );
    if( nx_s<0 || nx_s>=NX )                sf_error( "The ninline start is inappropriate!\n" );
    if( ny_s<0 || ny_s>=NY )                sf_error( "The ncrossline start is inappropriate!\n" );
    if( nx_s+nx>NX )                        sf_error( "The nx inversion is inappropriate!\n" );
    if( ny_s+ny>NY )                        sf_error( "The ny inversion is inappropriate!\n" );
    /*
     * READ & CHECK
     */

    NearWave = sf_floatalloc( nwave );
    Mid1Wave = sf_floatalloc( nwave );
    Mid2Wave = sf_floatalloc( nwave );
    TeleWave = sf_floatalloc( nwave );

    sf_floatread( NearWave , nwave , FW1 );
    sf_floatread( Mid1Wave , nwave , FW2 );
    sf_floatread( Mid2Wave , nwave , FW3 );
    sf_floatread( TeleWave , nwave , FW4 );

    MPOST = sf_floatalloc4( NT+1 , NX , NY , npara );

    zero4float( MPOST , NT+1 , NX , NY , npara );

    for( iy=0 ; iy<ny ; iy++ )
    for( ix=0 ; ix<nx ; ix++ )
    {
        nt  = NT;
        nt_s=0;

        d       = sf_floatalloc ( nsita*nt );
        M_prior = sf_floatalloc ( npara*(nt+1) );
        M_post  = sf_floatalloc ( npara*(nt+1) );
        M_true  = sf_floatalloc ( npara*(nt+1) );
        Cm      = sf_floatalloc2( npara*(nt+1) , npara*(nt+1) );
        Cn      = sf_floatalloc2( 4*nt , 4*nt );
        A       = sf_floatalloc2( npara*nt , nsita*nt );
        W       = sf_floatalloc2( nsita*nt , nsita*nt );
        D       = sf_floatalloc2( npara*(nt+1) , npara*nt );
        WA      = sf_floatalloc2( npara*nt , nsita*nt );
        G       = sf_floatalloc2( npara*(nt+1) , nsita*nt );

        sf_seek( FD1 , sizeof(float)*((ny_s+iy)*NX*NT+(nx_s+ix)*NT+nt_s) , SEEK_SET );
        sf_seek( FD2 , sizeof(float)*((ny_s+iy)*NX*NT+(nx_s+ix)*NT+nt_s) , SEEK_SET );
        sf_seek( FD3 , sizeof(float)*((ny_s+iy)*NX*NT+(nx_s+ix)*NT+nt_s) , SEEK_SET );
        sf_seek( FD4 , sizeof(float)*((ny_s+iy)*NX*NT+(nx_s+ix)*NT+nt_s) , SEEK_SET );

        sf_floatread( &d[0*nt] , nt , FD1 );
        sf_floatread( &d[1*nt] , nt , FD2 );
        sf_floatread( &d[2*nt] , nt , FD3 );
        sf_floatread( &d[3*nt] , nt , FD4 );

        sf_seek( FPR1 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );
        sf_seek( FPR2 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );
        sf_seek( FPR3 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );

        sf_floatread( &M_prior[0*(nt+1)] , nt+1 , FPR1 );
        sf_floatread( &M_prior[1*(nt+1)] , nt+1 , FPR2 );
        sf_floatread( &M_prior[2*(nt+1)] , nt+1 , FPR3 );

        sf_seek( FTR1 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );
        sf_seek( FTR2 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );
        sf_seek( FTR3 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );

        sf_floatread( &M_true[0*(nt+1)] , nt+1 , FTR1 );
        sf_floatread( &M_true[1*(nt+1)] , nt+1 , FTR2 );
        sf_floatread( &M_true[2*(nt+1)] , nt+1 , FTR3 );
        /*
         * READ DATASETS INCLUDING SEISMIC DATA, A PRIOR MODEL AND TRUE MODELS
         */

        for( it=0 ; it<nt ; it++ )
        {
            vpvsratio = (M_prior[nt+1+it]+M_prior[nt+1+it+1])/(M_prior[it]+M_prior[it+1]);

            A[0*nt+it][it]       = 0.5*(1.+tan(sita1*SF_PI/180.)*tan(sita1*SF_PI/180.));
            A[0*nt+it][nt+it]    = -4.*sin(sita1*SF_PI/180.)*sin(sita1*SF_PI/180.)*pow(vpvsratio,2.);
            A[0*nt+it][2*nt+it]  = 0.5*(1+A[0*nt+it][nt+it]);

            A[1*nt+it][it]       = 0.5*(1.+tan(sita2*SF_PI/180.)*tan(sita2*SF_PI/180.));
            A[1*nt+it][nt+it]    = -4.*sin(sita2*SF_PI/180.)*sin(sita2*SF_PI/180.)*pow(vpvsratio,2.);
            A[1*nt+it][2*nt+it]  = 0.5*(1+A[1*nt+it][nt+it]);

            A[2*nt+it][it]       = 0.5*(1.+tan(sita3*SF_PI/180.)*tan(sita3*SF_PI/180.));
            A[2*nt+it][nt+it]    = -4.*sin(sita3*SF_PI/180.)*sin(sita3*SF_PI/180.)*pow(vpvsratio,2.);
            A[2*nt+it][2*nt+it]  = 0.5*(1+A[2*nt+it][nt+it]);

            A[3*nt+it][it]       = 0.5*(1.+tan(sita4*SF_PI/180.)*tan(sita4*SF_PI/180.));
            A[3*nt+it][nt+it]    = -4.*sin(sita4*SF_PI/180.)*sin(sita4*SF_PI/180.)*pow(vpvsratio,2.);
            A[3*nt+it][2*nt+it]  = 0.5*(1+A[3*nt+it][nt+it]);
        }

        for( it=0 ; it<nt ; it++ )
        for( jt=0 ; jt<nt ; jt++ )
        {
            if( (it-jt)<=nwave-1-wsft && (it-jt)>=-wsft )
            {
                W[0*nt+it][0*nt+jt] = NearWave[it-jt+wsft];
                W[1*nt+it][1*nt+jt] = Mid1Wave[it-jt+wsft];
                W[2*nt+it][2*nt+jt] = Mid2Wave[it-jt+wsft];
                W[3*nt+it][3*nt+jt] = TeleWave[it-jt+wsft];
            }
            else
            {
                W[0*nt+it][0*nt+jt] = 0.;
                W[1*nt+it][1*nt+jt] = 0.;
                W[2*nt+it][2*nt+jt] = 0.;
                W[3*nt+it][3*nt+jt] = 0.;
            }
        }
        for( ipara=0 ; ipara<npara ; ipara++ )
        for( it=0    ; it<nt       ; it++ )
        {
            D[ipara*nt+it][ipara*(nt+1)+it]   = -1;
            D[ipara*nt+it][ipara*(nt+1)+it+1] =  1;
        }
        /*
         * @Build Matrix A,W,D
         * @G=WAD
         */
        lh_matrix_mu_matrix( W , A , WA , nsita*nt , nsita*nt , npara*nt );
        lh_matrix_mu_matrix( WA , D , G , nsita*nt , npara*nt , npara*(nt+1));

        for( it=0 ; it<(nt+1)*npara ; it++ )
        {
            M_prior[it] = log( M_prior[it] );
            M_true[it]  = log( M_true[it]  );
        }

        for( it=0 ; it<nsita*nt ; it++ )
        for( jt=0 ; jt<nsita*nt ; jt++ )
        {
            if( it==jt )
                Cn[it][jt] = var_e;
            else
                Cn[it][jt] = 0.;
        }
        for( it=0 ; it<nt+1 ; it++ )
        for( jt=0 ; jt<nt+1 ; jt++ )
        for( ipara=0 ; ipara<npara ; ipara++ )
        for( jpara=0 ; jpara<npara ; jpara++ )
        {
            if( ipara==jpara )
                Cm[ipara*(nt+1)+it][jpara*(nt+1)+jt] = lh_float_variance( &M_true[ipara*(nt+1)] , nt+1 )
                                                        *exp(-1.*pow(fabs((it-jt)*dt/crange),corder));
            else
                Cm[ipara*(nt+1)+it][jpara*(nt+1)+jt] = lh_float_covariance( &M_true[ipara*(nt+1)] , &M_true[jpara*(nt+1)] , nt+1 )
                                                        *exp(-1.*pow(fabs((it-jt)*dt/crange),corder));

        }
        /*
         * @ Building Cm and Cn
         */
        lh_direct_LS( G , nsita*nt , npara*(nt+1) , Cn , Cm , d, M_prior, M_post );
        /*
         * @Inversion Under Bayesian Framework, the formulae comes from Tarantola(2005) directly
         */
        for( it=0 ; it<nt+1 ; it++ )
        {
            MPOST[0][iy+ny_s][ix+nx_s][nt_s+it] = exp(M_post[0*(nt+1)+it]);
            MPOST[1][iy+ny_s][ix+nx_s][nt_s+it] = exp(M_post[1*(nt+1)+it]);
            MPOST[2][iy+ny_s][ix+nx_s][nt_s+it] = exp(M_post[2*(nt+1)+it]);
        }
        free1float( d );
        free1float( M_prior );
        free1float( M_post );
        free1float( M_true );
        free2float( Cm );
        free2float( Cn );
        free2float( A );
        free2float( W );
        free2float( D );
        free2float( WA );
        free2float( G);
    }
    if( dim==1 )
    {
        sf_oaxa( FPO1 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO2 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO3 , sf_maxa( NT+1 , 0., dt ) , 1 );
    }
    if( dim==2 )
    {
        sf_oaxa( FPO1 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO2 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO3 , sf_maxa( NT+1 , 0., dt ) , 1 );

        sf_oaxa( FPO1 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO2 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO3 , sf_maxa( NX , 0., 1 ) , 2 );
    }
    if( dim==3 )
    {
        sf_oaxa( FPO1 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO2 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO3 , sf_maxa( NT+1 , 0., dt ) , 1 );

        sf_oaxa( FPO1 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO2 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO3 , sf_maxa( NX , 0., 1 ) , 2 );

        sf_oaxa( FPO1 , sf_maxa( NY , 0., 1 ) , 3 );
        sf_oaxa( FPO2 , sf_maxa( NY , 0., 1 ) , 3 );
        sf_oaxa( FPO3 , sf_maxa( NY , 0., 1 ) , 3 );
    }
    sf_floatwrite( &MPOST[0][0][0][0] , (NT+1)*NX*NY , FPO1 );
    sf_floatwrite( &MPOST[1][0][0][0] , (NT+1)*NX*NY , FPO2 );
    sf_floatwrite( &MPOST[2][0][0][0] , (NT+1)*NX*NY , FPO3 );

	exit( 0 );
}
