/* Convert the Reflectivity to logarithmic Impedance by Trace Integration */

/*********************************************************************************
 * * Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
 * *
 * * File Name: Mlh_ref2imp.c
 * *
 * * Program Context: 
 * *
 * * Author: John Hush
 * *
 * * Version: 1.0
 * *
 * * Date: 2015/06/28
 * *
 * * History:
 * *
 * **********************************************************************************/

#include <rsf.h>
#include "su_alloc.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );

    sf_axis at, aa;

    sf_file FI=NULL, FO=NULL;
    int NT, it, dim, num_dim[SF_MAX_DIM], NTRACE, itrace, idim;
    float imp0, *ref, *imp;

    if(!sf_getfloat( "imp0" , &imp0 ))    sf_error( "Missing the value of the first impedance!!!\n" );
    /* Input the value of the first impedance imp0=velocity0*density0, could be 3000*2500 */

    FI = sf_input ( "in" );
    FO = sf_output( "out" );

    NTRACE = sf_leftsize( FI , 1 );
    dim    = sf_filedims( FI , num_dim );

    at = sf_iaxa( FI , 1 );
    NT = sf_n( at );

    ref = alloc1float( NT );
    imp = alloc1float( NT+1 );

    imp[0] = log(imp0);

    sf_setn( at , NT+1 );

    sf_oaxa( FO , at , 1 );

    if( dim>1 )
    {
        for( idim=2 ; idim<=dim ; idim++ )
        {
            at = sf_iaxa( FI , idim );
            sf_oaxa( FO , at , idim );
        }
    }

    for( itrace=0 ; itrace<NTRACE ; itrace++ )
    {
        sf_seek( FI , sizeof(float)*NT*itrace , SEEK_SET );
        sf_floatread( &ref[0] , NT , FI );

        for( it=1 ; it<NT+1 ; it++ )
            imp[it] = 2.*ref[it-1]+imp[it-1];

        sf_seek( FO , sizeof(float)*(NT+1)*itrace , SEEK_SET );
        sf_floatwrite( &imp[0] , NT+1 , FO );
    }

    exit( 0 );
}
