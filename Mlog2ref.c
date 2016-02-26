
#include <rsf.h>
#include "su_alloc.h"
#include "lh_readwrite.h"
#include "lh_wavelet.h"
#include "lh_bayesian.h"

int main( int argc , char* argv[] )
{
	int flag;
	int nt=99 , nsita = 40 , i , it , isita;
	float dt=0.002 , dsita=1. , osita=0. , ot , sita ,radian;
	float *logdata , *vp , *vs , *rou , **ref;
	sf_file Fr=NULL;
	float a, b , c;

	sf_init( argc , argv );

	if(!sf_getint( "flag" , &flag ))	sf_error( "Must specify flag!\n" );
	/* flag=1, welllog data to acoustic reflectivity; flag=2 , to elastic reflectivity */

	Fr = sf_output( "out" );

	logdata = sf_floatalloc ( 6*nt );
	vp      = sf_floatalloc ( nt );
	vs      = sf_floatalloc ( nt );
	rou     = sf_floatalloc ( nt );
	ref     = sf_floatalloc2( nt-1 , nsita );

	lh_read_1d_float_bin( logdata , 6*nt , "/home/luoheng/log_viscoelastic.bin" );

	ot = logdata[0];

	sf_putint  ( Fr , "n1" , nt-1 );
	sf_putint  ( Fr , "n2" , nsita );
	sf_putfloat( Fr , "d1" , dt );
	sf_putfloat( Fr , "d2" , dsita );
	sf_putfloat( Fr , "o1" , ot );
	sf_putfloat( Fr , "o2" , osita );

	for( i=99 ; i<396 ; i++ )
	{
		if( i>=99 && i<198 )
			vp[i-99] = logdata[i];
		else if( i>=198 && i<297 )
			vs[i-198] = logdata[i];
		else
			rou[i-297] = logdata[i];
	}

	for( isita=0 , sita=0. ; isita<nsita ; isita++, sita+=dsita )
	{
		radian = sita*SF_PI/180;
		for( it=1 ; it<nt ; it++ )
		{
			if(flag==2){
			a = 0.5*(1.+tan(radian)*tan(radian));
			b = -4.*sin(radian)*sin(radian)*(vs[it]/2+vs[it-1]/2)*(vs[it]/2+vs[it-1]/2)/((vp[it]/2+vp[it-1]/2)*(vp[it]/2+vp[it-1]/2));
			c = 0.5*(1.+b);
			}
			else if( flag==1 ){
			a = 0.5*(1.+tan(radian)*tan(radian));
			b = 0.;
			c = 0.5;
			}
			else{
			sf_error( "The flag is wrong!\n" );
			}

			ref[isita][it-1] = a*((vp[it]-vp[it-1])/(vp[it]/2+vp[it-1]/2))+b*((vs[it]-vs[it-1])/(vs[it]/2+vs[it-1]/2))+c*((rou[it]-rou[it-1])/(rou[it]/2+rou[it-1]/2));
		}
	}

	sf_floatwrite( &ref[0][0] , nsita*(nt-1) , Fr );

	exit( 0 );
}
