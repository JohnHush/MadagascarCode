/* AVF Forward Program Using Convolution Model in Frequency Domain. */

/**************************************************************************
* Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
*
* History
***************************************************************************/

#include <rsf.h>
#include "su_alloc.h"
#include "lh_readwrite.h"
#include "lh_wavelet.h"
#include "lh_bayesian.h"

int main( int argc , char* argv[] )
{
	sf_init( argc , argv );

    float *wavelet, *hwave, *fre;

    wavelet = alloc1float( 200 );
    hwave = alloc1float( 200 );
    fre = alloc1float( 200 );

    lh_unit_ricker_nonzero_phase( wavelet , 200 , 0.002 , 25 , 0.2 , 330 );

    lh_write_1d_float_bin( wavelet , 200 , "waveler.bin" );

    lh_signal_phase( wavelet , hwave , 1 , 200 );

    lh_write_1d_float_bin( hwave , 200 , "phase.bin" );


    lh_signal_frequency( hwave , fre , 0.002 , 200  );

    int i;
    for( i=0 ; i<200 ; i++)
        fre[i] = fre[i]*sqrt(SF_PI)*0.5;

    lh_write_1d_float_bin( fre , 200 , "freq.bin" );



	exit( 0 );
}
