 /* unstretched NMO based on matching pursuit method using ricker wavelet */

#include <rsf.h>
#include "su_alloc.h"
#include "lh_wavelet.h"
#include "lh_readwrite.h"
#include "lh_mp.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    sf_file FI1 , FI2 , FO1;
    sf_axis at, afold, ax, ay;
    int dim, num_dim[SF_MAX_DIM], ncdp, icdp, nt, it, nfold, ifold, iter, niter, iphase, imax, max_phase, tindex, max_enve;
    float dt, ofold, dfold,***gather_in, ***gather_nmo, **nmo_vel, *enve, *freq, *phas, *resi, eps, ins_freq, *wavelet, corr, max_corr;

    if(!sf_getint  ( "niter" , &niter    ))  sf_error( "Missing NITER!\n" );
    if(!sf_getfloat( "eps_iphase" , &eps ))  sf_error( "Missing eps_iphase!\n" );

    FI1 = sf_input ( "in" );
    FI2 = sf_input ( "nmo_vel" );
    FO1 = sf_output( "nmo_gather" );

    dim = sf_filedims( FI1 , num_dim );
    ncdp = sf_leftsize( FI1 , 2 );

    at      = sf_iaxa( FI1 , 1 );
    afold   = sf_iaxa( FI1 , 2 );

    nt      = sf_n( at );
    dt      = sf_d( at );
    nfold   = sf_n( afold );
    dfold   = sf_d( afold );
    ofold   = sf_o( afold );

    gather_in   = alloc3float( nt , nfold , ncdp );
    gather_nmo  = alloc3float( nt , nfold , ncdp );
    nmo_vel     = alloc2float( nt , ncdp );
    enve        = alloc1float( nt );
    freq        = alloc1float( nt );
    phas        = alloc1float( nt );
    resi        = alloc1float( nt );
    wavelet     = alloc1float( nt );

    zero3float( gather_in , nt , nfold , ncdp );
    zero3float( gather_nmo, nt , nfold , ncdp );
    zero2float( nmo_vel   , nt , ncdp );

    sf_floatread( &gather_in[0][0][0] , nt*nfold*ncdp , FI1 );
    sf_floatread( &nmo_vel[0][0]      , nt*ncdp       , FI2 );

    for( icdp=0  ; icdp<ncdp   ; icdp++ )
    for( ifold=0 ; ifold<nfold ; ifold++ )
    for( iter=0  ; iter<niter  ; iter++ )
    {
        if( iter==0 )
        {
            for( it=0 ; it<nt ; it++ )
                resi[it] = gather_in[icdp][ifold][it];
        }
        zero1float( enve , nt );
        zero1float( phas , nt );
        zero1float( freq , nt );

        lh_signal_envelope ( resi , enve , nt );
        lh_signal_phase    ( resi , phas , eps , nt );
        lh_signal_frequency( phas , freq , dt , nt );

        max_enve = lh_index_max( enve , nt );
        ins_freq = -freq[max_enve]*sqrt(SF_PI)/2.;

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

        printf( "maxphase=%d\tifold=%d\titer=%d\tins_freq=%f\tmaxcrr=%f\n" , max_phase , ifold , iter , ins_freq , max_corr);
        lh_unit_ricker_nonzero_phase( wavelet , nt , dt , ins_freq , dt*(max_enve) , (float)max_phase );

        for( it=0 ; it<nt ; it++ )
            resi[it] = resi[it] - wavelet[it]*max_corr;

        if( ifold==0 && iter==niter-1 )
            lh_write_1d_float_bin( resi , nt , "resi.bin" );

        if( dt*max_enve*dt*max_enve>=((ofold+dfold*ifold)*(ofold+dfold*ifold))/(nmo_vel[icdp][max_enve]*nmo_vel[icdp][max_enve]) )
            tindex =(int) ((sqrt( dt*max_enve*dt*max_enve - 
                            ((ofold+dfold*ifold)*(ofold+dfold*ifold))/(nmo_vel[icdp][max_enve]*nmo_vel[icdp][max_enve]) ))/dt);
        else
        {
            tindex = 0;
            zero1float( wavelet , nt );
        }

        lh_shift_array( wavelet , nt , max_enve-tindex );

        for( it=0 ; it<nt ; it++ )
            gather_nmo[icdp][ifold][it] += wavelet[it]*max_corr;
    }

    sf_oaxa( FO1 , at   , 1 );
    sf_oaxa( FO1 , afold, 2 );
    if( dim==3 )
    {
        ax = sf_iaxa( FI1 , 3 );
        sf_oaxa( FO1 , ax , 3 );
    }
    if( dim==4 )
    {
        ax = sf_iaxa( FI1 , 3 );
        ay = sf_iaxa( FI1 , 4 );
        sf_oaxa( FO1 , ax , 3 );
        sf_oaxa( FO1 , ay , 4 );
    }
    sf_floatwrite( &gather_nmo[0][0][0] , nt*nfold*ncdp , FO1 );

    exit( 0 );
}
