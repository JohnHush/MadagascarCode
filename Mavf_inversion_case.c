/* Quality Factor of P wave and S wave inversion Program */

/*********************************************************************************
* Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
*
* File Name: Mavf_inversion_casestudy.c
*
* Program Context: Bayesian Inversion scheme in ANACOUSTIC/ANELASTIC Media
*
* Others: 
*
* Author: Long Teng( Originally designed in 2013 summer).
*
* Revisor:  Heng Luo ( Revised to C code in 2014 Spring).
*
* Version: 
*
* Date: 2014/03/27
*
* History:
* 
* Reference.:   1.  Buland and Omre, 2003, Baysian Linearized AVO inversion
*               2.  Innanen, 2011, Anelastic AVF/AVA inversion
*               3.  Heng Luo and Huazhong Wang, 2013, A New Linearized AVAF Expression in Anelastic Media 
*
**********************************************************************************/

#include <rsf.h>
#include "lh_readwrite.h"
#include "lh_wavelet.h"
#include "lh_bayesian.h"
#include "tl_bayesian.h"
#include "lh_LSCG.h"
#include "su_alloc.h"

int main( int argc , char* argv[] )
{
	sf_init( argc , argv );

    sf_axis at, ax, ay, aw, at_wave, aw_wave;
    /*
     * @Define AXIS for DATASETS
     */
    sf_file FD1, FD2, FD3, FD4, FW1, FW2, FW3, FW4, FPR1, FPR2, FPR3, FPR4, FPR5, FPO1, FPO2, FPO3, FPO4, FPO5, FLA1, FLA2;
    /*
     * @Define File Interfaces for Datasets
     */
    int ix, nx, nx_s, NX, iy, ny, ny_s, NY, iw, nw, nw_s, NW, it, jt, nt, nt_s, NT, nwave;
    /*
     * @Define the Dimension of the Input Data & Inversion
     * @it  : index of time axis
     * @nx  : Dimension of Inversion in Inline Direction
     * @nx_s: Starting place of Inline Direction
     * @NX  : Dimension of Input Data
     * @NY,ny,ny_s: Parameters in Crossline Direction
     * @NT,nt,nt_s: Parameters in time direction
     * @nwave: length of the wavelet
     */
    float sita1, sita2, sita3, sita4;
    /*
     * @ Angle Values of NearOffset, Mid1, Mid2 and TeleOffset, respectively!
     */
    float w0, crange, corder, var_e;
    int wsft;
    /*
     * @ Parameters of Inversion, w0 is the reference frequency of the Q model
     * @ crange and corder are the parameters of Gaussian Shaped Window function
     * @ var_e is the level of noise estimated
     * @ wsft is the wave shift of the wavelet to the original point
     */
    int dim, num_dim[SF_MAX_DIM];
    /*
     * @ Detect the Dimension of the DATASETS Using the Interior Function
     */
    float dt, ot, dw, ow;
    /*
     * @ PARAMETER of the DATASETS
     * @ dt is the time interval, ot is the starting time value
     */
    float **NearWave, **Mid1Wave, **Mid2Wave, **TeleWave;
    /*
     * @ Input Wavelet file in different frequency and Offset
     */
    float **upper, **lower;
    /*
     * @ Define Two layers to control inversion interval
     */
    float *d, *M_prior, **Cm, **Cn, **A, **W, **D, **WA, **G, *M_post, ****MPOST;
    /*
     * @
     * @
     * @Cm is the covariance matrix of the model
     * @Cn is the covariance matrix of the noise
     */
    int ipara, jpara;
    float vpvsratio;
    /*
     * @auxiliary variable in calculating the coefficient in AKI-RICHARD formulae
     */

    if(!sf_getint( "nw_start" , &nw_s ))             sf_error( "Missing Frequency Start!\n" );
    /* Input Parameter: the starting frequency in inversion */
    if(!sf_getint( "nw_inv" , &nw ))                 sf_error( "Missing the Number of Frequency used!\n" );
    /* Input Parameter: the number of frequency in inversion */
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
    if(!sf_getfloat( "w0" , &w0 ))                   sf_error( "Missing w0!\n" );
    /* Input Parameter: Reference Frequency in Inversion */
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
     * @ Input DATASETs
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
    FPR4= sf_input( "qp_prior" );
    FPR5= sf_input( "qs_prior" );
    /*
     * @ Input a Prior Models
     */
    FLA1= sf_input( "upper_layer" );
    FLA2= sf_input( "lower_layer" );
    /*
     * @ Input Upper layer and Lower layer
     */
    FPO1= sf_output( "vp_post" );
    FPO2= sf_output( "vs_post" );
    FPO3= sf_output( "ro_post" );
    FPO4= sf_output( "qp_post" );
    FPO5= sf_output( "qs_post" );
    /*
     * @ Output Post MODELs by Bayesian Inversion
     */

    dim = sf_filedims( FD1 , num_dim );

    if( dim==3 )
    {
        at = sf_iaxa( FD1 , 1 );
        NT = sf_n( at );
        dt = sf_d( at );
        ot = sf_o( at );

        aw = sf_iaxa( FD1 , 2 );
        NW = sf_n( aw );
        dw = sf_d( aw );
        ow = sf_o( aw );

        ax = sf_iaxa( FD1 , 3 );
        NX = sf_n( ax );

        NY = 1;
    }
    else if( dim==4 )
    {
        at = sf_iaxa( FD1 , 1 );
        NT = sf_n( at );
        dt = sf_d( at );
        ot = sf_o( at );

        aw = sf_iaxa( FD1 , 2 );
        NW = sf_n( aw );
        dw = sf_d( aw );
        ow = sf_o( aw );

        ax = sf_iaxa( FD1 , 3 );
        NX = sf_n( ax );

        ay = sf_iaxa( FD1 , 4 );
        NY = sf_n( ay );
    }
    else
        sf_error( "DATA DIMENSION FAULT!\n" );
    /*
     * @ Only Handle the problem owns at least an inline section
     * @ dimension of the data MUST be 3 or 4!
     */

    at_wave = sf_iaxa( FW1 , 1 );
    aw_wave = sf_iaxa( FW1 , 2 );
    nwave   = sf_n( at_wave );

    if( sf_n(aw_wave)!=NW )                 sf_error( "The NW dimension of wavelet should equal to the one of the DATASETS!\n" );
    if( sf_n(sf_iaxa( FPR1 , 1 ))!=NT+1 )   sf_error( "DATASETS UNMATCHED!\n" );
    if( sf_n(sf_iaxa( FPR1 , 2 ))!=NX )     sf_error( "DATASETS UNMATCHED!\n" );
    if( nw_s<0 || nw_s>=NW )                sf_error( "The nw_start is inappropriate!\n" );
    if( nx_s<0 || nx_s>=NX )                sf_error( "The ninline start is inappropriate!\n" );
    if( ny_s<0 || ny_s>=NY )                sf_error( "The ncrossline start is inappropriate!\n" );
    if( nw_s+nw>NW )                        sf_error( "The nw inversion is inappropriate!\n" );
    if( nx_s+nx>NX )                        sf_error( "The nx inversion is inappropriate!\n" );
    if( ny_s+ny>NY )                        sf_error( "The ny inversion is inappropriate!\n" );

    /*
     * @Read Parameter of Wavelet & Check Parameters 
     * @if occurs inappropriate parameter, the program breaks immediately
     */

    NearWave = sf_floatalloc2( nwave , nw );
    Mid1Wave = sf_floatalloc2( nwave , nw );
    Mid2Wave = sf_floatalloc2( nwave , nw );
    TeleWave = sf_floatalloc2( nwave , nw );

    sf_seek( FW1 , sizeof(float)*nwave*nw_s , SEEK_SET );
    sf_seek( FW2 , sizeof(float)*nwave*nw_s , SEEK_SET );
    sf_seek( FW3 , sizeof(float)*nwave*nw_s , SEEK_SET );
    sf_seek( FW4 , sizeof(float)*nwave*nw_s , SEEK_SET );

    sf_floatread( &NearWave[0][0] , nwave*nw , FW1 );
    sf_floatread( &Mid1Wave[0][0] , nwave*nw , FW2 );
    sf_floatread( &Mid2Wave[0][0] , nwave*nw , FW3 );
    sf_floatread( &TeleWave[0][0] , nwave*nw , FW4 );
    /*
     * @Clip the WAVELET DATASETS and Read Needed Wavelet
     */
    tl_normalize2d( NearWave , NearWave , nwave , nw );
    tl_normalize2d( Mid1Wave , Mid1Wave , nwave , nw );
    tl_normalize2d( Mid2Wave , Mid2Wave , nwave , nw );
    tl_normalize2d( TeleWave , TeleWave , nwave , nw );

    upper = sf_floatalloc2( NX , NY );
    lower = sf_floatalloc2( NX , NY );
    /*
     * @ layers dimension set to be NY*NX
     */
    sf_floatread( &upper[0][0] , NX*NY , FLA1 );
    sf_floatread( &lower[0][0] , NX*NY , FLA2 );

    MPOST   = sf_floatalloc4( NT+1 , NX , NY , 5 );

    zero4float( MPOST , NT+1 , NX , NY , 5 );
    /*
     * @ The MAP solution of all the PARAMETERS
     */

    for( iy=0 ; iy<ny ; iy++ )
    for( ix=0 ; ix<nx ; ix++ )
    {
        nt = (int)((lower[iy][ix]-upper[iy][ix])/dt);
        nt_s=(int)((upper[iy][ix]-ot)/dt);
        printf( "ix=%d\tnt=%d\tnt_s=%d\n" , ix , nt , nt_s );

        d       = sf_floatalloc( 4*nw*nt );
        M_prior = sf_floatalloc ( 5*(nt+1) );
        M_post  = sf_floatalloc ( 5*(nt+1) );
        Cm      = sf_floatalloc2( 5*(nt+1) , 5*(nt+1) );
        Cn      = sf_floatalloc2( 4*nw*nt , 4*nw*nt );
        A       = sf_floatalloc2( 5*nt , 4*nw*nt );
        W       = sf_floatalloc2( 4*nw*nt , 4*nw*nt );
        D       = sf_floatalloc2( 5*(nt+1) , 5*nt );
        WA      = sf_floatalloc2( 5*nt , 4*nw*nt );
        G       = sf_floatalloc2( 5*(nt+1) , 4*nw*nt );
        /*
         * @4 means 4 angle gather
         * @5 means 5 parameter,
         * @4 &5 could be set up in the parameter list if it's needed!
         * @To keep simplicity, the two parameters are given explicitly
         */
        
        for( iw=0 ; iw<nw ; iw++ )
        {
            sf_seek( FD1 , sizeof(float)*((ny_s+iy)*NX*NW*NT+(nx_s+ix)*NW*NT+(nw_s+iw)*NT+nt_s) , SEEK_SET );
            sf_seek( FD2 , sizeof(float)*((ny_s+iy)*NX*NW*NT+(nx_s+ix)*NW*NT+(nw_s+iw)*NT+nt_s) , SEEK_SET );
            sf_seek( FD3 , sizeof(float)*((ny_s+iy)*NX*NW*NT+(nx_s+ix)*NW*NT+(nw_s+iw)*NT+nt_s) , SEEK_SET );
            sf_seek( FD4 , sizeof(float)*((ny_s+iy)*NX*NW*NT+(nx_s+ix)*NW*NT+(nw_s+iw)*NT+nt_s) , SEEK_SET );

            sf_floatread( &d[(0*nw+iw)*nt] , nt , FD1 );
            sf_floatread( &d[(1*nw+iw)*nt] , nt , FD2 );
            sf_floatread( &d[(2*nw+iw)*nt] , nt , FD3 );
            sf_floatread( &d[(3*nw+iw)*nt] , nt , FD4 );
        }

        tl_normalize( d , d , 4*nw*nt );

        sf_seek( FPR1 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );
        sf_seek( FPR2 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );
        sf_seek( FPR3 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );
        sf_seek( FPR4 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );
        sf_seek( FPR5 , sizeof(float)*((ny_s+iy)*NX*(NT+1)+(nx_s+ix)*(NT+1)+nt_s) , SEEK_SET );

        sf_floatread( &M_prior[0*(nt+1)] , nt+1 , FPR1 );
        sf_floatread( &M_prior[1*(nt+1)] , nt+1 , FPR2 );
        sf_floatread( &M_prior[2*(nt+1)] , nt+1 , FPR3 );
        sf_floatread( &M_prior[3*(nt+1)] , nt+1 , FPR4 );
        sf_floatread( &M_prior[4*(nt+1)] , nt+1 , FPR5 );
        /*
         * @Read Angle freuquency Gather use FILE SEEK function
         * @Read a Prior MODEL
         */

        for( iw=0 ; iw<nw ; iw++ )
        for( it=0 ; it<nt ; it++ )
        {
            vpvsratio = (M_prior[nt+1+it]+M_prior[nt+1+it+1])/(M_prior[it]+M_prior[it+1]);

            A[0*nw*nt+iw*nt+it][it]      = 0.5*(1.+tan(sita1*SF_PI/180.)*tan(sita1*SF_PI/180.));
            A[0*nw*nt+iw*nt+it][nt+it]   = -4.*sin(sita1*SF_PI/180.)*sin(sita1*SF_PI/180.)*pow(vpvsratio,2.);
            A[0*nw*nt+iw*nt+it][2*nt+it] = 0.5*(1.-4.*sin(sita1*SF_PI/180.)*sin(sita1*SF_PI/180.)*pow(vpvsratio,2.));
            A[0*nw*nt+iw*nt+it][3*nt+it] = A[0*nw*nt+iw*nt+it][it]*((1./SF_PI)*log((dw*(nw_s+iw)+ow)/w0));
            A[0*nw*nt+iw*nt+it][4*nt+it] = A[0*nw*nt+iw*nt+it][nt+it]*((1./SF_PI)*log((dw*(nw_s+iw)+ow)/w0));

            A[1*nw*nt+iw*nt+it][it]      = 0.5*(1.+tan(sita2*SF_PI/180.)*tan(sita2*SF_PI/180.));
            A[1*nw*nt+iw*nt+it][nt+it]   = -4.*sin(sita2*SF_PI/180.)*sin(sita2*SF_PI/180.)*pow(vpvsratio,2.);
            A[1*nw*nt+iw*nt+it][2*nt+it] = 0.5*(1.-4.*sin(sita2*SF_PI/180.)*sin(sita2*SF_PI/180.)*pow(vpvsratio,2.));
            A[1*nw*nt+iw*nt+it][3*nt+it] = A[1*nw*nt+iw*nt+it][it]*((1./SF_PI)*log((dw*(nw_s+iw)+ow)/w0));
            A[1*nw*nt+iw*nt+it][4*nt+it] = A[1*nw*nt+iw*nt+it][nt+it]*((1./SF_PI)*log((dw*(nw_s+iw)+ow)/w0));

            A[2*nw*nt+iw*nt+it][it]      = 0.5*(1.+tan(sita3*SF_PI/180.)*tan(sita3*SF_PI/180.));
            A[2*nw*nt+iw*nt+it][nt+it]   = -4.*sin(sita3*SF_PI/180.)*sin(sita3*SF_PI/180.)*pow(vpvsratio,2.);
            A[2*nw*nt+iw*nt+it][2*nt+it] = 0.5*(1.-4.*sin(sita3*SF_PI/180.)*sin(sita3*SF_PI/180.)*pow(vpvsratio,2.));
            A[2*nw*nt+iw*nt+it][3*nt+it] = A[2*nw*nt+iw*nt+it][it]*((1./SF_PI)*log((dw*(nw_s+iw)+ow)/w0));
            A[2*nw*nt+iw*nt+it][4*nt+it] = A[2*nw*nt+iw*nt+it][nt+it]*((1./SF_PI)*log((dw*(nw_s+iw)+ow)/w0));

            A[3*nw*nt+iw*nt+it][it]      = 0.5*(1.+tan(sita4*SF_PI/180.)*tan(sita4*SF_PI/180.));
            A[3*nw*nt+iw*nt+it][nt+it]   = -4.*sin(sita4*SF_PI/180.)*sin(sita4*SF_PI/180.)*pow(vpvsratio,2.);
            A[3*nw*nt+iw*nt+it][2*nt+it] = 0.5*(1.-4.*sin(sita4*SF_PI/180.)*sin(sita4*SF_PI/180.)*pow(vpvsratio,2.));
            A[3*nw*nt+iw*nt+it][3*nt+it] = A[3*nw*nt+iw*nt+it][it]*((1./SF_PI)*log((dw*(nw_s+iw)+ow)/w0));
            A[3*nw*nt+iw*nt+it][4*nt+it] = A[3*nw*nt+iw*nt+it][nt+it]*((1./SF_PI)*log((dw*(nw_s+iw)+ow)/w0));
        }

        for( iw=0; iw<nw; iw++ )
        for( it=0; it<nt; it++ )
        for( jt=0; jt<nt; jt++ )
        {
            if ( (it-jt)<=nwave-1-wsft && (it-jt)>=-wsft )
            {
                W[0*nw*nt+iw*nt+it][0*nw*nt+iw*nt+jt]= NearWave[iw][it-jt+wsft];
                W[1*nw*nt+iw*nt+it][1*nw*nt+iw*nt+jt]= Mid1Wave[iw][it-jt+wsft];
                W[2*nw*nt+iw*nt+it][2*nw*nt+iw*nt+jt]= Mid2Wave[iw][it-jt+wsft];
                W[3*nw*nt+iw*nt+it][3*nw*nt+iw*nt+jt]= TeleWave[iw][it-jt+wsft];
            }
            else
            {
                W[0*nw*nt+iw*nt+it][0*nw*nt+iw*nt+jt]=0.;
                W[1*nw*nt+iw*nt+it][1*nw*nt+iw*nt+jt]=0.;
                W[2*nw*nt+iw*nt+it][2*nw*nt+iw*nt+jt]=0.;
                W[3*nw*nt+iw*nt+it][3*nw*nt+iw*nt+jt]=0.;
            }
        }

        for( ipara=0 ; ipara<5 ; ipara++ )
        for( it=0    ; it<nt   ; it++ )
        {
            D[ipara*nt+it][ipara*(nt+1)+it]   = -1;
            D[ipara*nt+it][ipara*(nt+1)+it+1] = 1;
        }
        /*
         * @Build three essential Matrix, A,W,D
         * @A is the coefficient Matrix in multiple angle(4) and multiple frequency form
         * @W is the wavelet matrix
         * @D is the difference matrix
         * @the forward operator is G=WAD
         */
        lh_matrix_mu_matrix( W , A , WA , 4*nw*nt , 4*nw*nt , 5*nt );
        lh_matrix_mu_matrix( WA , D , G , 4*nw*nt , 5*nt , 5*(nt+1));

        /*
         * @ Build forward operator G, using W,A,D
         */
        for( it=0 ; it<3*(nt+1) ; it++ )
        {
            M_prior[it] = log(M_prior[it]);
        }
        for( it=3*(nt+1) ; it<5*(nt+1) ; it++ )
        {
            M_prior[it] = 1./M_prior[it];
        }
        /*
         * @ Revise the A Prior Model Into the expectation value of True Value
         */

        for( it=0 ; it<nt+1 ; it++ )
        for( jt=0 ; jt<nt+1 ; jt++ )
        for( ipara=0 ; ipara<5 ; ipara++ )
        for( jpara=0 ; jpara<5 ; jpara++ )
        {
            if( ipara==jpara )
                Cm[ipara*(nt+1)+it][jpara*(nt+1)+jt] = lh_float_variance( &M_prior[ipara*(nt+1)] , nt+1 )
                                                        *exp(-1.*pow(fabs((it-jt)*dt/crange),corder));
            else
                Cm[ipara*(nt+1)+it][jpara*(nt+1)+jt] = lh_float_covariance( &M_prior[ipara*(nt+1)] , &M_prior[jpara*(nt+1)] , nt+1 )
                                                        *exp(-1.*pow(fabs((it-jt)*dt/crange),corder));

        }
        for( it=0 ; it<4*nw*nt ; it++ )
        for( jt=0 ; jt<4*nw*nt ; jt++ )
        {
            if( it==jt )
                Cn[it][jt] = var_e;
            else
                Cn[it][jt] = 0.;
        }
        /*
         * @ Transfer the ELASTIC parameters into LOG value
         * @ Transfer the ANELASTIC parameters into RECIPROCAL value
         * @ Setup the variance of the MODEL and the NOISE,
         * @ Cm Usually comes from welldata, here is given from MODELS
         */
        lh_direct_LS( G , 4*nw*nt , 5*(nt+1) , Cn , Cm , d, M_prior, M_post );
        /*
         * @ Inversion Under Bayesian Framework, the formulae comes from Tarantola(2005) directly
         */
        if( ix==0 )
        {
            lh_write_1d_float_bin( M_post , 5*(nt+1) , "M_post.bin" );
            lh_write_1d_float_bin( M_prior, 5*(nt+1) , "M_prior.bin");
        }

        for( it=0 ; it<nt+1 ; it++ )
        {
            MPOST[0][iy+ny_s][ix+nx_s][nt_s+it] = exp(M_post[0*(nt+1)+it]);
            MPOST[1][iy+ny_s][ix+nx_s][nt_s+it] = exp(M_post[1*(nt+1)+it]);
            MPOST[2][iy+ny_s][ix+nx_s][nt_s+it] = exp(M_post[2*(nt+1)+it]);
            MPOST[3][iy+ny_s][ix+nx_s][nt_s+it] = 1.0/M_post[3*(nt+1)+it];
            MPOST[4][iy+ny_s][ix+nx_s][nt_s+it] = 1.0/M_post[4*(nt+1)+it];
        }

        free1float( d );
        free1float( M_prior );
        free1float( M_post );
        free2float( Cm );
        free2float( Cn );
        free2float( A );
        free2float( W );
        free2float( D );
        free2float( WA );
        free2float( G );
    }
    if( dim==3 )
    {
        sf_oaxa( FPO1 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO2 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO3 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO4 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO5 , sf_maxa( NT+1 , 0., dt ) , 1 );

        sf_oaxa( FPO1 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO2 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO3 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO4 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO5 , sf_maxa( NX , 0., 1 ) , 2 );
    }
    if( dim==4 )
    {
        sf_oaxa( FPO1 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO2 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO3 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO4 , sf_maxa( NT+1 , 0., dt ) , 1 );
        sf_oaxa( FPO5 , sf_maxa( NT+1 , 0., dt ) , 1 );

        sf_oaxa( FPO1 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO2 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO3 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO4 , sf_maxa( NX , 0., 1 ) , 2 );
        sf_oaxa( FPO5 , sf_maxa( NX , 0., 1 ) , 2 );

        sf_oaxa( FPO1 , sf_maxa( NY , 0., 1 ) , 3 );
        sf_oaxa( FPO2 , sf_maxa( NY , 0., 1 ) , 3 );
        sf_oaxa( FPO3 , sf_maxa( NY , 0., 1 ) , 3 );
        sf_oaxa( FPO4 , sf_maxa( NY , 0., 1 ) , 3 );
        sf_oaxa( FPO5 , sf_maxa( NY , 0., 1 ) , 3 );
    }
    sf_floatwrite( &MPOST[0][0][0][0] , (NT+1)*NX*NY , FPO1 );
    sf_floatwrite( &MPOST[1][0][0][0] , (NT+1)*NX*NY , FPO2 );
    sf_floatwrite( &MPOST[2][0][0][0] , (NT+1)*NX*NY , FPO3 );
    sf_floatwrite( &MPOST[3][0][0][0] , (NT+1)*NX*NY , FPO4 );
    sf_floatwrite( &MPOST[4][0][0][0] , (NT+1)*NX*NY , FPO5 );

    free( MPOST );

    exit( 0 );
}
