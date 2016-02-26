/* decompose the trace into several individual wavelets by Matching pursuit method
 * ======================================================================
 * the number of wavelets is specified to be constant 
 * the program is developed by Luo Heng in Nanjing, CN, 23th, Dec, 2015*/

#include <rsf.h>
#include "su_alloc.h"
#include "lh_wavelet.h"
#include "lh_readwrite.h"
#include "lh_mp.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    sf_file FI1 , FO1, FO2;
    sf_axis at, ax, ay;
    int dim, num_dim[SF_MAX_DIM], ncdp, icdp, nt, it, iter, niter, iphase,  max_phase, max_enve;

    float dt,**gather_in, **gather_out, **layer_position, *enve, *freq, *phas, *resi, eps, ins_freq, *wavelet, corr, max_corr;

    if(!sf_getint  ( "niter" , &niter    ))  sf_error( "Missing NITER!\n" );
    if(!sf_getfloat( "eps_iphase" , &eps ))  sf_error( "Missing eps_iphase!\n" );

    FI1 = sf_input ( "in" );
    FO1 = sf_output( "out" );
    FO2 = sf_output( "layer_position" );

    dim = sf_filedims( FI1 , num_dim );
    ncdp = sf_leftsize( FI1 , 1 );

    at      = sf_iaxa( FI1 , 1 );

    nt      = sf_n( at );
    dt      = sf_d( at );

    gather_in   = alloc2float( nt , ncdp );
    gather_out  = alloc2float( nt , ncdp );
    enve        = alloc1float( nt );
    freq        = alloc1float( nt );
    phas        = alloc1float( nt );
    resi        = alloc1float( nt );
    wavelet     = alloc1float( nt );

    layer_position = alloc2float( nt , ncdp );

    zero2float( gather_in  , nt , ncdp );
    zero2float( gather_out , nt , ncdp );
    zero2float( layer_position , nt , ncdp );

    sf_floatread( &gather_in[0][0] , nt*ncdp , FI1 );

    for( icdp=0  ; icdp<ncdp   ; icdp++ )
    for( iter=0  ; iter<niter  ; iter++ )
    {
        if( iter==0 )
        {
            for( it=0 ; it<nt ; it++ )
                resi[it] = gather_in[icdp][it];
        }
        zero1float( enve , nt );
        zero1float( phas , nt );
        zero1float( freq , nt );

        lh_signal_envelope ( resi , enve , nt );
        lh_signal_phase    ( resi , phas , eps , nt );
        lh_signal_frequency( phas , freq , dt , nt );

        max_enve = lh_index_max( enve , nt );
        ins_freq = -freq[max_enve]*sqrt(SF_PI)/2.;

        if ( (int)layer_position[icdp][max_enve] == 0 )
            layer_position[icdp][max_enve] = 1.;

        for( iphase=0 ; iphase<60 ; iphase=iphase+3 )
        {
            corr = 0.;
            lh_unit_ricker_nonzero_phase( wavelet , nt , dt , ins_freq , dt*max_enve , (float)iphase );

            for( it=0 ; it<nt ; it++ )
                corr += wavelet[it]*resi[it];
            if( iphase==0 )
            {
                max_phase= 0;
                max_corr = corr;
            }
            if( fabs(corr) > fabs(max_corr) )
            {
                max_phase= iphase;
                max_corr = corr;
            }
        }

//        printf( "maxphase=%d\tifold=%d\titer=%d\tins_freq=%f\tmaxcrr=%f\n" , max_phase , ifold , iter , ins_freq , max_corr);
        lh_unit_ricker_nonzero_phase( wavelet , nt , dt , ins_freq , dt*(max_enve) , (float)max_phase );

        for( it=0 ; it<nt ; it++ )
            resi[it] = resi[it] - wavelet[it]*max_corr;

        for( it=0 ; it<nt ; it++ )
            gather_out[icdp][it] += wavelet[it]*max_corr;
    }

    sf_oaxa( FO1 , at   , 1 );
    sf_oaxa( FO2 , at   , 1 );
    if( dim==2 )
    {
        ax = sf_iaxa( FI1 , 2 );
        sf_oaxa( FO1 , ax , 2 );
        sf_oaxa( FO2 , ax , 2 );
    }
    if( dim==3 )
    {
        ax = sf_iaxa( FI1 , 2 );
        ay = sf_iaxa( FI1 , 3 );
        sf_oaxa( FO1 , ax , 2 );
        sf_oaxa( FO2 , ax , 2 );
        sf_oaxa( FO1 , ay , 3 );
        sf_oaxa( FO2 , ay , 3 );
    }
    sf_floatwrite( &gather_out[0][0] , nt*ncdp , FO1 );
    sf_floatwrite( &layer_position[0][0] , nt*ncdp , FO2 );

    exit( 0 );
}
