/* AVF Forward Program Using Convolution Model in Frequency Domain. */

/**************************************************************************

* Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
*
* File Name: Mavf_forward.c
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
* Revised: 2014/02/17
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

	int it, isita;
	int nt, nsita, nwave;
	float dt, dsita, t0, fm, amp;
	float *wavelet, **trace;
	sf_file Fr=NULL, Ft=NULL, Fw=NULL;
	sf_axis at , asita, awave;
	kiss_fftr_cfg cfg, icfg;
	kiss_fft_cpx *ori_spec, *res_spec;
	int nfft , nw , iw;
	float *res_wavelet , dw, ***rvf;
	int itao;
	sf_axis aw;

	if(!sf_getfloat( "fm" , &fm ))			fm=35.;
	/* The default main frequency of the wavelet is 50HZ */
	if(!sf_getfloat("amp", &amp ))			amp=10.;
	/* The default amp of the wavelet is 10. */

	Fr = sf_input ( "input" );
	Ft = sf_output( "out" );
	Fw = sf_output( "wavelet" );

	at    = sf_iaxa( Fr , 1 );   nt    = sf_n( at );      dt    = sf_d( at );
	aw    = sf_iaxa( Fr , 2 );   nw    = sf_n( aw );      dw    = sf_d( aw );
	asita = sf_iaxa( Fr , 3 );   nsita = sf_n( asita );   dsita = sf_d( asita );

	sf_oaxa( Ft , at    , 1 );
	sf_oaxa( Ft , asita , 2 );

	nwave = 1./( dw*dt );
	if(nwave%2) nwave++;
	nfft = nwave/2+1;

    t0 = dt*(nwave/2);

	awave = sf_maxa( nwave , 0. , dt );
	sf_setlabel( awave , "time" );
	sf_setunit( awave , "s" );

	sf_oaxa( Fw , awave , 1 );

	trace	    = sf_floatalloc2( nt , nsita );
	wavelet	    = sf_floatalloc ( nwave );
	res_wavelet = sf_floatalloc( nwave );
	rvf	    = sf_floatalloc3( nt , nw , nsita );
	ori_spec    = (kiss_fft_cpx *)sf_complexalloc( nfft );
	res_spec    = (kiss_fft_cpx *)sf_complexalloc( nfft );
	cfg	    = kiss_fftr_alloc( nwave , 0 , NULL , NULL );
	icfg	    = kiss_fftr_alloc( nwave , 1 , NULL , NULL );

	zero2float( trace , nt , nsita );

	sf_floatread( &rvf[0][0][0] , nt*nw*nsita , Fr );
	lh_ricker( wavelet, amp, nwave, dt, fm, t0 );
	sf_floatwrite( wavelet , nwave , Fw );
	kiss_fftr( cfg , wavelet , ori_spec );
	/* Calculate the Fourier spectrum of the original wavelet */

	for( isita=0 ; isita<nsita ; isita++ )
	for( it=0    ; it<nt       ; it++ )
	{
		zero1float( res_wavelet , nwave );
		for( iw=0 ; iw<nfft ; iw++ )
		{
			res_spec[iw].r = 0.;
			res_spec[iw].i = 0.;
		}
		/* Clean the result spectrum arrays first before we calculate
                        the calculated spectrum of the reflected wavelet */
		for( iw=0 ; iw<(nfft>nw?nw:nfft) ; iw++ )
			res_spec[iw+1] = sf_crmul( ori_spec[iw+1] , rvf[isita][iw][it] );

		kiss_fftri( icfg , res_spec , res_wavelet );
		lh_shift_array( res_wavelet , nwave , (int)(t0/dt+0.5)-it );
		/* Time-shift the wavelet to simulate the propagation of the wave */

		for( itao=0 ; itao<nwave ; itao++ )
			res_wavelet[itao] = res_wavelet[itao]/nwave;

		for( itao=0 ; itao<(nt<nwave?nt:nwave) ; itao++ )
			trace[isita][itao]+=res_wavelet[itao];
	}
	sf_floatwrite( &trace[0][0] , nsita*nt , Ft );
	exit( 0 );
}
