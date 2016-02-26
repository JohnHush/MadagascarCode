#include <rsf.h>
#include "su_alloc.h"
#include "lh_readwrite.h"
#include "lh_wavelet.h"
#include "lh_bayesian.h"

int main( int argc , char* argv[] )
{
	int flag, nt=99, nsita=40, nw=100, i, it, isita, iw;
	float dt=0.002, dsita=1., osita=0., ow=0., dw, ot, sita, radian, fre , w0=30.;
	float **logdata, *vp, *vs, *rou, *qp, *qs, ***ref, *vp1, *vs1 , *rou1, *qp1 , *qs1;
	float a, b, c, d, e;
	float **log2 , sum;
	sf_file Fr=NULL , Fo1=NULL , Fo2=NULL, Fo3=NULL, Fo4=NULL, Fo5=NULL;
	int num=1;

	dw = 1./(dt*512);

	sf_init( argc , argv );

	if(!sf_getint( "flag" , &flag ))	sf_error( "Must specify flag!\n" );

	Fr = sf_output( "data" );
	Fo1= sf_output( "vp" );
	Fo2= sf_output( "vs" );
	Fo3= sf_output( "rou" );
	Fo4= sf_output( "qp" );
	Fo5= sf_output( "qs" );

	logdata = sf_floatalloc2( nt , 6 );
	vp      = sf_floatalloc ( nt );
	vs      = sf_floatalloc ( nt );
	rou     = sf_floatalloc ( nt );
	qp      = sf_floatalloc ( nt );
	qs      = sf_floatalloc ( nt );

	log2	= sf_floatalloc2( nt , 6 );
	vp1     = sf_floatalloc ( nt );
	vs1     = sf_floatalloc ( nt );
	rou1    = sf_floatalloc ( nt );
	qp1     = sf_floatalloc ( nt );
	qs1     = sf_floatalloc ( nt );
	ref     = sf_floatalloc3( nt-1 , nw , nsita );

	zero1float( vp1 , nt );
	zero1float( vs1 , nt );
	zero1float( rou1, nt );
	zero1float( qp1 , nt );
	zero1float( qs1 , nt );
	zero2float( log2, nt , 6 );

	lh_read_2d_float_bin( logdata , 6 , nt , "/home/Pitaloveu/DeskTop/log_viscoelastic.bin" );

	ot = logdata[0][0];

	sf_putint  ( Fr , "n1" , nt-1 );
	sf_putint  ( Fr , "n2" , nw );
	sf_putint  ( Fr , "n3" , nsita );
	sf_putfloat( Fr , "d1" , dt );
	sf_putfloat( Fr , "d2" , dw );
	sf_putfloat( Fr , "d3" , dsita );
	sf_putfloat( Fr , "o1" , ot );
	sf_putfloat( Fr , "o2" , ow );
	sf_putfloat( Fr , "o3" , osita );
	sf_putstring( Fr , "label1" , "Time" );
	sf_putstring( Fr , "label2" , "Frequency" );
	sf_putstring( Fr , "label3" , "Angle" );
	sf_putstring( Fr , "unit1" , "ms" );
	sf_putstring( Fr , "unit2" , "hz" );
	sf_putstring( Fr , "unit3" , "degree" );

	for( i=0 ; i<nt ; i++ )
	{
		vp[i]  = logdata[1][i];
		vs[i]  = logdata[2][i];
		rou[i] = logdata[3][i];
		qp[i]  = logdata[4][i];
		qs[i]  = logdata[5][i];
	}
	for( it=0 ; it<nt ; it++ )
	{
		vp1[it]  = vp[it/num];
		vs1[it]  = vs[it/num];
		rou1[it] = rou[it/num];
		qp1[it]  = qp[it/num];
		qs1[it]  = qs[it/num];
	}
	for( i=0 ; i<nt ; i++ )
	{
		log2[1][i] = vp1[i];
		log2[2][i] = vs1[i];
		log2[3][i] = rou1[i];
		log2[4][i] = qp1[i];
		log2[5][i] = qs1[i];
	}
	lh_write_2d_float_bin( log2 , 6 , nt , "/home/Pitaloveu/DeskTop/log_viscoelastic_resampling.bin" );

	for( isita=0 , sita=0. ; isita<nsita ; isita++, sita+=dsita )
	for( iw=1 ,fre=dw ; iw<nw ; iw++ , fre+=dw )
	{
		radian = sita*SF_PI/180;
		for( it=1 ; it<nt ; it++ )
		{
			if( flag==3 )
			{
				a = 0.5*(1.+tan(radian)*tan(radian));
				b = 0.;
				c = 0.5;
				d = a*((1./SF_PI)*log(fre/w0));
				e = 0.;
			}
			if( flag==5 )
			{
				a = 0.5*(1.+tan(radian)*tan(radian));
				b = -4.*sin(radian)*sin(radian)*(vs1[it]/2+vs1[it-1]/2)*(vs1[it]/2+vs1[it-1]/2)/((vp1[it]/2+vp1[it-1]/2)*(vp1[it]/2+vp1[it-1]/2));
				c = 0.5*(1.+b);
				d = a*((1./SF_PI)*log(fre/w0));
				e = b*((1./SF_PI)*log(fre/w0));
			}

			ref[isita][iw][it-1] = a*((vp1[it]-vp1[it-1])/(vp1[it]/2+vp1[it-1]/2))+
					   b*((vs1[it]-vs1[it-1])/(vs1[it]/2+vs1[it-1]/2))+
					   c*((rou1[it]-rou1[it-1])/(rou1[it]/2+rou1[it-1]/2))+
					   d*((qp1[it])-(qp1[it-1]))+
					   e*((qs1[it])-(qs1[it-1]));
		}
	}
	sf_floatwrite( &ref[0][0][0], nsita*nw*(nt-1) , Fr );
	sf_floatwrite( vp1 , nt , Fo1 );
	sf_floatwrite( vs1 , nt , Fo2 );
	sf_floatwrite( rou1, nt , Fo3 );
	sf_floatwrite( qp1 , nt , Fo4 );
	sf_floatwrite( qs1 , nt , Fo5 );

	exit( 0 );
}
