#include <rsf.h>
#include <complex.h>
#include <fftw3.h>
#include "su_alloc.h"
#include "st_avf.h"
#include "lh_bayesian.h"
#include "lh_wavelet.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    sf_file FI, FO;

    int dim, num_dim[SF_MAX_DIM], NT , NTRACE, nfft, itrace, it;
    float **data_in, *s_real, *s_imag, *phase, *amp, **data_out;
    fftw_complex *fft_in, *fft_out, *fft_tmp;
    fftw_plan p1,p2;

    FI = sf_input( "in" );
    FO = sf_output( "out" );

    dim = sf_filedims( FI , num_dim );

    NT = sf_n( sf_iaxa( FI , 1 ) );
    NTRACE = sf_leftsize( FI , 1 );


    data_in = alloc2float( NT , NTRACE );
    data_out= alloc2float( NT , NTRACE );
    phase   = alloc1float( NT );
    amp     = alloc1float( NT );
    fft_in = (fftw_complex*) fftw_malloc( sizeof(fftw_complex) *NT );
    fft_out= (fftw_complex*) fftw_malloc( sizeof(fftw_complex) *NT );
    fft_tmp= (fftw_complex*) fftw_malloc( sizeof(fftw_complex) *NT );

    zero2float( data_in , NT , NTRACE );
    zero2float( data_out, NT , NTRACE );
    zero1float( phase   , NT );
    zero1float( amp     , NT );
    for( it=0 ; it<NT ; it++ )
        fft_tmp[it] =0;

    sf_floatread( &data_in[0][0] , NT*NTRACE , FI );

    for( itrace=0 ; itrace<30 ; itrace++ )
    {
        for( it=0 ; it<NT ; it++ )
            fft_in[it] = data_in[itrace][it];

        p1 = fftw_plan_dft_1d( NT , fft_in , fft_out , FFTW_FORWARD, FFTW_ESTIMATE );
        fftw_execute( p1 );
        for( it=0 ; it<NT ; it++ )
            fft_tmp[it]+= fft_out[it];
    }
    for( it=0 ; it<NT ; it++ )
        phase[it] = atan2( cimag(fft_tmp[it]) , creal(fft_tmp[it]) );

    for( itrace=0 ; itrace<NTRACE ; itrace++ )
    {
        
        for( it=0 ; it<NT ; it++ )
            fft_in[it] = data_in[itrace][it];

        p1 = fftw_plan_dft_1d( NT , fft_in , fft_out , FFTW_FORWARD, FFTW_ESTIMATE );
        fftw_execute( p1 );

        for( it=0 ; it<NT ; it++ )
            amp[it] = sqrt( creal(fft_out[it])*creal(fft_out[it]) + cimag(fft_out[it])*cimag(fft_out[it]) );

        for( it=0 ; it<NT ; it++ )
        {
            fft_in[it] = amp[it]*cexp( phase[it]*I );
        }
        p2 = fftw_plan_dft_1d( NT , fft_in , fft_out , FFTW_BACKWARD , FFTW_ESTIMATE );
        fftw_execute( p2 );

        for( it=0 ; it<NT ; it++ )
            data_out[itrace][it] = creal(fft_out[it]);
    }
    sf_floatwrite( &data_out[0][0] , NT*NTRACE , FO );

    fftw_destroy_plan( p1 );
    fftw_destroy_plan( p2 );
    fftw_free( fft_in );
    fftw_free( fft_out );

    return 0;
}
