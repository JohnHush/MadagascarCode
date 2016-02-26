/* Generate attenuated traces using 1D--5D reflectivity & Q values with the same wavelet */

/*********************************************************************************
 * * Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
 * *
 * * File Name: MqConvolution.c
 * *
 * * Program Context: 
 * *
 * * Author: John Hush
 * *
 * * Version: 1.0
 * *
 * * Date: 2015/06/15
 * *
 * * History:
 * *
 * **********************************************************************************/

#include <rsf.h>
#include "su_alloc.h"
#include "lh_wavelet.h"
#include "lh_bayesian.h"
#include "lh_readwrite.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    sf_file FI1, FI2, FI3, FO1;

    sf_axis aa;

    int dim, num_dim[SF_MAX_DIM], NTRACE, NT, nw, nwave, wsft, itrace, idim;

    float dt, **ref, **Q, *wave, **gather, **TVM;

    if(!sf_getint( "wsft" , &wsft ))  sf_error( "MISSING wave_shift PARAMETER!!\n" );
    /* The wave shift should be set up! */

    FI1 = sf_input( "ref" );
    FI2 = sf_input( "Q" );
    FI3 = sf_input( "wave" );
    FO1 = sf_output( "gather" );

    NTRACE = sf_leftsize( FI1 , 1 );
    dim = sf_filedims( FI1 , num_dim );

    /* define the LEFTSIZE bigger than the first 
     * dimension, so we can treat the data as two
     * dimensional-data, 
     */

    NT = sf_n( sf_iaxa( FI1 , 1 ) );
    dt = sf_d( sf_iaxa( FI1 , 1 ) );
    nw = sf_n( sf_iaxa( FI3 , 1 ) );

    /* extract the information about the first
     * dimension, which is needed by the function
     * of calculating the attenuated wavelet
     */

    nwave = lh_powerof2( nw );

    /* compute the suitable length for the FFT
     */

    ref     = alloc2float( NT , NTRACE );
    Q       = alloc2float( NT , NTRACE );
    gather  = alloc2float( NT , NTRACE );
    wave    = alloc1float( nwave );
    TVM     = alloc2float( NT , NT );

    zero2float( ref   , NT , NTRACE );
    zero2float( Q     , NT , NTRACE );
    zero2float( gather, NT , NTRACE );
    zero1float( wave  , nwave );

    /* alloc and setup */

    sf_floatread( &ref[0][0] , NT*NTRACE , FI1 );
    sf_floatread( &Q[0][0]   , NT*NTRACE , FI2 );
    sf_floatread( &wave[0]   , nw        , FI3 );

    for( itrace=0 ; itrace<NTRACE ; itrace++ )
    {
        zero2float( TVM , NT , NT );
        lh_time_variant_matrix( wave , wsft , nwave , dt , &Q[itrace][0], NT , TVM );
        lh_matrix_mu_vector( TVM , &ref[itrace][0] , &gather[itrace][0], NT , NT );
    }
/*
float **spec_real, **spec_imag, **spec_amp;
int nwave1 = lh_powerof2( NT );

spec_real = alloc2float( nwave1 , NT );
spec_imag = alloc2float( nwave1 , NT );
spec_amp  = alloc2float( 200 , NT );

zero2float( spec_real , nwave1 , NT );
zero2float( spec_imag , nwave1 , NT );
zero2float( spec_amp  , 200   , NT );

int i,j;
for( i=0 ; i<NT ; i++ )
{
    for( j=0 ; j<NT ; j++ )
    spec_real[i][j] = TVM[j][i];

    lh_fft( &spec_real[i][0] , &spec_imag[i][0], nwave1 , 1 );

    for( j=0 ; j<200 ; j++ )
    spec_amp[i][j] = sqrt( spec_real[i][j]*spec_real[i][j]+spec_imag[i][j]*spec_imag[i][j] );
}
lh_write_2d_float_bin_row( spec_amp , NT , 200 , "spec_amp_qwave.bin" );
*/

    for( idim=1 ; idim<=dim ; idim++ )
    {
        aa = sf_iaxa( FI1 , idim );
        sf_oaxa( FO1 , aa , idim );
    }
    sf_floatwrite( &gather[0][0] , NT*NTRACE , FO1 );

    exit( 0 );
}
