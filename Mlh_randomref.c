/* Random reflectivity Generator */

/*********************************************************************************
 * * Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
 * *
 * * File Name: Mrandomref.c
 * *
 * * Program Context: Generate pseudo-random reflectivity using rand() function 
 * *
 * * Author: John Hush
 * *
 * * Version: 1.0
 * *
 * * Date: 2015/06/09
 * *
 * * History:
 * *
 * **********************************************************************************/

#include <rsf.h>
#include "su_alloc.h"
#include <time.h>

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    sf_axis at;

    sf_file FO=NULL;

    int nt, it, seed;

    float dt, *ref, density, ot, maxref;

    if(!sf_getint( "nt" , &nt ))            sf_error( "Missing the Length of the reflectivity!!!\n" );
    /* Input the LENGTH of the reflectivity */
    if(!sf_getfloat( "dt" , &dt ))          sf_error( "Missing the Sample rate of the reflectivity!!!\n" );
    /* Input the SAMPLE RATE of the reflectivity */
    if(!sf_getfloat( "ot" , &ot ))          sf_error( "Missing the Initial time of the reflectivity!!!\n" );
    /* Input the INITIAL TIME of the reflectivity */
    if(!sf_getfloat( "density" , &density ))density=0.5; 
    /* density is the approximate density of the non-zero reflection, valued from 0.0 to 1.0 */
    if(!sf_getfloat( "maxref" , &maxref ))  maxref=0.8;
    /* value of the maximum reflectivity could be specfied from 0.0 to 1.0 corresponds to the reality */
    if(!sf_getint( "seed" , &seed ))        srand((unsigned)time(0)); 
    /* Specify the SEEDS of the random function, if not been given, then use the current time */
    else
        srand( seed );

    FO = sf_output( "out" );

    ref = alloc1float( nt );

    zero1float( ref , nt );

    for( it=0 ; it<nt ; it++ )
    {
        if( ((float)rand()/RAND_MAX)<=density )
            ref[it] = 2.*maxref*((float)rand()/RAND_MAX)-maxref;
    }

    at = sf_maxa( nt , ot , dt );
    sf_oaxa( FO , at , 1 );
    sf_floatwrite( ref , nt , FO );
}
