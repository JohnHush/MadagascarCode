/* Generate the time-frequency spectrum of 1D-5D data by S TRANSFORM */

/*********************************************************************************
 * * Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
 * *
 * * File Name: Mlh_timefreq.c
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
#include "st_avf.h"
#include "lh_bayesian.h"
#include "lh_wavelet.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    sf_file FI1, FO;

    sf_axis aa;

    int dim, num_dim[SF_MAX_DIM], NTRACE, NT, itrace, idim, nw, nfft, iw, it, iw_index;

    float dt, ow, dw, *s_real, *s_imag, *tmp, **data, ***spec_amp;

    if(!sf_getint( "nw" , &nw ))    sf_error( "missing nw!!\n" );
    /* the number of frequency in TIME-FREQUENCY transformation */
    if(!sf_getfloat( "dw" , &dw ))    sf_error( "missing dw!!\n" );
    /* the INTERVAL of frequency in TIME-FREQUENCY transformation */
    if(!sf_getfloat( "ow" , &ow ))    sf_error( "missing ow!!\n" );
    /* the STARTING frequency in TIME-FREQUENCY transformation */

    FI1 = sf_input( "in" );
    FO = sf_output( "out" );

    NTRACE = sf_leftsize( FI1 , 1 );
    dim = sf_filedims( FI1 , num_dim );

    /* define the LEFTSIZE bigger than the first 
     * dimension, so we can treat the data as two
     * dimensional-data, 
     */

    NT = sf_n( sf_iaxa( FI1 , 1 ) );
    dt = sf_d( sf_iaxa( FI1 , 1 ) );

    nfft = lh_powerof2( NT );

    /* compute the suitable length for the FFT
     */

    data   = alloc2float( NT , NTRACE );
    s_real = alloc1float( nfft );
    s_imag = alloc1float( nfft );
    tmp    = alloc1float( nfft );

    spec_amp = alloc3float( NT , nw , NTRACE );

    zero2float( data , NT , NTRACE );
    zero3float( spec_amp , NT , nw , NTRACE );

    /* alloc and setup */

    sf_floatread( &data[0][0] , NT*NTRACE , FI1 );

    for( itrace=0 ; itrace<NTRACE ; itrace++ )
    {
        zero1float( tmp , nfft );
        for( it=0 ; it<NT ; it++ )
            tmp[it] = data[itrace][it];
        for( iw=0 ; iw<nw ; iw++ )
        {
            iw_index = (int)((ow+dw*iw)*nfft*dt);
            zero1float( s_real , nfft );
            zero1float( s_imag , nfft );
            gst_fre( tmp , s_real , s_imag , dt , nfft , iw_index );
            for( it=0 ; it<NT ; it++ )
            {
                spec_amp[itrace][iw][it] = sqrt( s_real[it]*s_real[it]+s_imag[it]*s_imag[it] );
            }
        }
    }

    aa = sf_iaxa( FI1 , 1 );
    sf_oaxa( FO , aa , 1 );

    aa = sf_maxa( nw , ow , dw );
    sf_oaxa( FO , aa , 2 );

    for( idim=2 ; idim<=dim ; idim++ )
    {
        aa = sf_iaxa( FI1 , idim );
        sf_oaxa( FO , aa , idim+1 );
    }
    sf_floatwrite( &spec_amp[0][0][0] , nw*NT*NTRACE , FO );

    exit( 0 );
}
