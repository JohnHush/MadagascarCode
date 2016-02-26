/* compute RMS velocity for layered model using formulae */

#include <rsf.h>
#include "su_alloc.h"

int main( int argc , char* argv[] )
{
    sf_init( argc , argv );
    sf_file Fi1, Fi2, Fo1, Fo2;
    sf_axis at, afold, ax, ay;
    int dim, num_dim[SF_MAX_DIM];
    int ncdp, icdp, nfold, ifold, nt , nz , iz, it, jz, index, iindex, tindex;
    float dt, dz , dfold, ofold;
    float **layered_vel, **rms_vel, **grid_time, sum_time, ***nmo_gather, ***ori_gather;

    Fi1 = sf_input( "input_gather" );   /* input gather before nmo correction*/
    Fi2 = sf_input( "input_velocity" ); /* input layered velocity */
    Fo1 = sf_output( "output_gather" ); /* output gather after nmo correction */
    Fo2 = sf_output( "output_velocity" );/* output rms velocity */

    dim = sf_filedims( Fi1 , num_dim); 
    ncdp = sf_leftsize( Fi1 , 2 );

    nt    = sf_n( sf_iaxa( Fi1 , 1 ) );
    dt    = sf_d( sf_iaxa( Fi1 , 1 ) );
    nfold = sf_n( sf_iaxa( Fi1 , 2 ) );
    dfold = sf_d( sf_iaxa( Fi1 , 2 ) );
    ofold = sf_o( sf_iaxa( Fi1 , 2 ) );
    nz    = sf_n( sf_iaxa( Fi2 , 1 ) );
    dz    = sf_d( sf_iaxa( Fi2 , 1 ) );

    layered_vel = alloc2float( nz , ncdp );
    rms_vel     = alloc2float( nt , ncdp );
    grid_time   = alloc2float( nz , ncdp );
    nmo_gather  = alloc3float( nt , nfold , ncdp );
    ori_gather  = alloc3float( nt , nfold , ncdp );

    zero2float( layered_vel , nz , ncdp );
    zero2float( rms_vel     , nt , ncdp );
    zero2float( grid_time   , nz , ncdp );
    zero3float( nmo_gather  , nt , nfold , ncdp );
    zero3float( ori_gather  , nt , nfold , ncdp );

    sf_floatread( &layered_vel[0][0]   , ncdp*nz       , Fi2 );
    sf_floatread( &ori_gather[0][0][0] , nt*nfold*ncdp , Fi1 );

    for( icdp=0 ; icdp<ncdp ; icdp++ )
    {
        for( iz=0 ; iz<nz ; iz++ )
            grid_time[icdp][iz] = 2.*dz/layered_vel[icdp][iz];
        
        for( it=0 ; it<nt ; it++ )
        {
            index = 0 ;
            sum_time = grid_time[icdp][0];
            for( iz=1 ; iz<nz ; iz++ )
            {
                if( (it*dt)<sum_time )
                    break;

                sum_time += grid_time[icdp][iz];
                index += 1;
            }
            for( iindex=0 ; iindex<index ; iindex++ )
            {
                rms_vel[icdp][it] += grid_time[icdp][iindex]*layered_vel[icdp][iindex]*layered_vel[icdp][iindex];
            }
            if( index!=nz-1 )
                rms_vel[icdp][it] += layered_vel[icdp][index]*layered_vel[icdp][index]*(grid_time[icdp][index]-sum_time+it*dt);
            else
                rms_vel[icdp][it] += layered_vel[icdp][index]*layered_vel[icdp][index]*(it*dt-sum_time);

            if( it!=0 )
            {
                rms_vel[icdp][it] /= (it*dt);
                rms_vel[icdp][it] = sqrt( rms_vel[icdp][it] );
            }
            else
                rms_vel[icdp][it] = layered_vel[icdp][it];
        }
    }

    for( icdp=0 ; icdp<ncdp   ; icdp++ )
    for( it=0   ; it<nt       ; it++ )
    for( ifold=0; ifold<nfold ; ifold++ )
    {
        tindex = (int)((sqrt( dt*it*dt*it + ((dfold*ifold+ofold)*(dfold*ifold+ofold))/(rms_vel[icdp][it]*rms_vel[icdp][it]) ) 
                    - (dt*it))/dt) + it;

        if( tindex<0 || tindex>=nt )
            nmo_gather[icdp][ifold][it] = 0.;
        else
            nmo_gather[icdp][ifold][it] = ori_gather[icdp][ifold][tindex];
    }

    at    = sf_iaxa( Fi1 , 1 );
    afold = sf_iaxa( Fi1 , 2 );
    if( dim==3 )
        ax= sf_iaxa( Fi1 , 3 );
    if( dim==4 )
    {
        ax= sf_iaxa( Fi1 , 3 );
        ay= sf_iaxa( Fi1 , 4 );
    }

    sf_oaxa( Fo1 , at    , 1 );
    sf_oaxa( Fo1 , afold , 2 );
    sf_oaxa( Fo2 , at    , 1 );
    if( dim==3 )
    {
        sf_oaxa( Fo1 , ax , 3 );
        sf_oaxa( Fo2 , ax , 2 );
    }
    if( dim==4 )
    {
        sf_oaxa( Fo1 , ax , 3 );
        sf_oaxa( Fo1 , ay , 4 );
        sf_oaxa( Fo2 , ax , 2 );
        sf_oaxa( Fo2 , ay , 3 );
    }

    sf_floatwrite( &nmo_gather[0][0][0] , nt*nfold*ncdp , Fo1 );
    sf_floatwrite( &rms_vel[0][0]       , nt*ncdp , Fo2 );

    exit( 0 );
}
