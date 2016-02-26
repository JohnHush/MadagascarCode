/* Generate ricker wavelet */

/*********************************************************************************
 * * Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
 * *
 * * File Name: Mlh_ricker.c
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
#include "lh_wavelet.h"
#include "su_alloc.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    sf_file FO1;

    sf_axis aa;

    float dt, amp, *ricker, fm;
    
    int nt, wsft;

    if(!sf_getint( "wsft" , &wsft) )    sf_error( "MISSING wavelet shift!\n" );
    /* The wavelet shift should be setup! */
    if(!sf_getint( "nt" , &nt) )        sf_error( "MISSING wavelet shift!\n" );
    /* The length should be setup! */
    if(!sf_getfloat( "dt" , &dt ))      sf_error( "MIssing dt!\n" );
    /* dt should be set up! */
    if(!sf_getfloat( "amp" , &amp ))    sf_error( "MIssing amp!\n" );
    /* Amplitude should be set up! */
    if(!sf_getfloat( "fm" , &fm ))      sf_error( "MIssing fm!\n" );
    /* Main frequency should be set up! */

    FO1 = sf_output( "out" );

    ricker = alloc1float( nt );

    zero1float( ricker , nt );

    lh_ricker( ricker , amp , nt , dt , fm , dt*wsft );

    aa = sf_maxa( nt , 0. , dt );

    sf_oaxa( FO1 , aa , 1 );

    sf_floatwrite( &ricker[0] , nt , FO1 );

    exit( 0 );
}
