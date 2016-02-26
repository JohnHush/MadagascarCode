/* AVA Forward Program Using Convolution Model in Time Domain. */

/**************************************************************************

* Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
*
* File Name: Mava_forward.c
*
* Program Context: Construst seismic trace using Convolution model in Frequency domain
*
* Others: input a RVF data & parameters of the source wavelet is optional
*
* Author: Heng Luo
*
* Version: 1.0.2
*
* Date: 2014/02/11
* Revised: 2014/02/18
*
* History
***************************************************************************/

#include <rsf.h>
#include "su_alloc.h"
#include "lh_readwrite.h"
#include "lh_wavelet.h"
#include "lh_bayesian.h"

int main( int argc , char* argv[] )
{
	sf_init( argc , argv );

	int isita;
	int nt, nsita, nwave, wsft;
	float dt, dsita, fm, amp;
	float *wavelet, **trace;
	sf_file Fr=NULL, Ft=NULL, Fw=NULL;
	sf_axis at , asita, awave ;
	float **rvf;

    if(!sf_getint( "wave_shift" , &wsft ))  wsft=40;
    /* Get Wave Shift Parameter , default value id 40 */
	if(!sf_getfloat( "fm" , &fm ))			fm=50.;
	/* The default main frequency of the wavelet is 50HZ */
	if(!sf_getfloat("amp", &amp ))			amp=10.;
	/* The default amp of the wavelet is 10. */

	Fr = sf_input ( "input" );
	Ft = sf_output( "out" );
    Fw = sf_output( "wavelet" );

	at    = sf_iaxa( Fr , 1 );   nt    = sf_n( at );      dt    = sf_d( at );
	asita = sf_iaxa( Fr , 2 );   nsita = sf_n( asita );   dsita = sf_d( asita );

	sf_oaxa( Ft , at    , 1 );
	sf_oaxa( Ft , asita , 2 );

	nwave = 256;

    awave = sf_maxa( nwave , 0., dt );
    sf_setlabel( awave , "time" );
    sf_setunit( awave , "s" );
    sf_oaxa( Fw , awave , 1 );

	trace   = sf_floatalloc2( nt , nsita );
	wavelet = sf_floatalloc ( nwave );
	rvf	= sf_floatalloc2( nt , nsita );

	sf_floatread( &rvf[0][0] , nt*nsita , Fr );
	lh_ricker( wavelet, amp, nwave, dt, fm, wsft*dt );

	for( isita=0 ; isita<nsita ; isita++ )
		lh_convolution_cut( &rvf[isita][0] , nt , wavelet , nwave , &trace[isita][0] , wsft );

	sf_floatwrite( &trace[0][0] , nsita*nt , Ft );
    sf_floatwrite( wavelet , nwave , Fw );

	exit( 0 );

}
