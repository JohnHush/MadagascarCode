/* Windowed FFT to get Frequency dependent Trace gather */

/**********************************************************************************
 * Copyright(C), By SINOPEC Geophysical Research Institute, Nanjing, CN
 *
 * File Name: Mwinft.c
 *
 * Authors: Heng Luo & Long Teng
 *
 * Date: 2014/11/03
 *
 * ********************************************************************************/

#include <rsf.h>
#include "st_avf.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    float ***data, ****data_out;

    int ix, nx, iy, ny, it, nt, nfft, nw, iw;
    /*
     * @ Setup the parameters of the dimensions of the DATASETS
     * @ the Output Parameters also be setup
     */
    sf_file Fi=NULL, Fo=NULL;
    /*
     * @ Setup the FILE POINT of Input Data and Output Frequency Dependent Trace
     */
    sf_axis at, ax, ay, aw;
    /*
     * @ Setup the AXISes of DataSETS
     */
    float var;
    /*
     * @ Variance of the Window Function in FT
     */
    int dim, num_dim[SF_MAX_DIM];

    float dt, dw, nw_s, fre;

    if(!sf_getint  ( "nw" , &nw ))         sf_error( "Missing nw!\n" );
    /* Input Parameter: the Number of Slicing in FD */
    if(!sf_getfloat( "nw_s" , &nw_s))      sf_error( "Missing nw_s!\n" );
    /* Input Parameter: the starting number of frequency */
    if(!sf_getfloat( "dw" , &dw ))         sf_error( "Missing dw!\n" );
    /* Input Parameter: the interval of OUTPUT frequency */
    if(!sf_getint( "nfft" , &nfft ))       sf_error( "Missing nfft!\n" );
    /* Input Parameter: the length of FFT */
    if(!sf_getfloat( "variance" , &var ))  sf_error( "Missing variance!\n" );
    /* the variance of the window */

    Fi = sf_input ( "data_in" );
    Fo = sf_output( "data_out" );

    dim = sf_filedims( Fi , num_dim );

    if( dim==1 )
    {
        at = sf_iaxa( Fi , 1 );
        nt = sf_n( at );
        dt = sf_d( at );

        nx = 1;
        ny = 1;
    }
    else if( dim==2 )
    {
        at = sf_iaxa( Fi , 1 );
        nt = sf_n( at );
        dt = sf_d( at );

        ax = sf_iaxa( Fi , 2 );
        nx = sf_n( ax );

        ny = 1;
    }
    else if( dim==3 )
    {
        at = sf_iaxa( Fi , 1 );
        nt = sf_n( at );
        dt = sf_d( at );

        ax = sf_iaxa( Fi , 2 );
        nx = sf_n( ax );

        ay = sf_iaxa( Fi , 3 );
        ny = sf_n( ay );
    }
    else
    {
        sf_error( "DIMENSION WRONG!\n" );
    }

    data        = sf_floatalloc3( nt , nx , ny );
    data_out    = sf_floatalloc4( nt , nw , nx , ny );

    sf_floatread( &data[0][0][0] , nt*nx*ny , Fi );

    for( iy=0 ; iy<ny ; iy++ )
    for( ix=0 ; ix<nx ; ix++ )
    for( iw=0 ; iw<nw ; iw++ )
    {
        fre = nw_s+dw*iw;
        winft( &data[iy][ix][0] , &data_out[iy][ix][iw][0] , nt , nfft , dt , fre , var );
    }

    aw = sf_maxa( nw , nw_s , dw );

    sf_oaxa( Fo , at , 1 );
    sf_oaxa( Fo , aw , 2 );
    if( dim==2 )
    {
        sf_oaxa( Fo , ax , 3 );
    }
    if( dim==3 )
    {
        sf_oaxa( Fo , ax , 3 );
        sf_oaxa( Fo , ay , 4 );
    }

    sf_floatwrite( &data_out[0][0][0][0] , ny*nx*nw*nt , Fo );

    free( data );
    free( data_out );

    exit( 0 );
}
