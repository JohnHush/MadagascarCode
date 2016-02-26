/* Matching Pursuit Program in Q media to extract RVF */

/*********************************************************************************
* Copyright(C), SINOPEC Geophysical Research Institute, Nanjing, CN
*
* File Name: Mmp_rvf.c
*
* Program Context: Matching Pursuit Method based on linearized AVF equation
*
* Others: 
*
* Author:   Heng Luo ( Originally designed in 2013 winter ).
*
* Revisor:  Heng Luo ( Revised in Feb,Mar/2014 ).
*
* Version: 1.0.2
*
* Date: 2014/03/27
*
* History:
*         version 1.0.2: 2014/03/27     ( Revised by Heng Luo )
* 

**********************************************************************************/

#include <rsf.h>
#include "lh_mp.c"
#include "lh_readwrite.c"
#include "su_alloc.c"
#include "lh_wavelet.c"

int main( int argc , char* argv[] )
{
	sf_init( argc , argv );
	sf_file Fd=NULL, Fw=NULL, Fo=NULL;
	sf_axis at, aw, asita, ainline, acrossline, awave;

	int nt, nsita, ninline, ncrossline, nfft, nfre;
	int i, dim, num_dim[SF_MAX_DIM];
	int nsita_in, ninline_in, ncrossline_in, nsita_start, ninline_start, ncrossline_start;
	int it, isita, iinline, icrossline, iw;
	int nk, nevent, flag, nwave, wave_shift;

	float ****data, *wavelet, *****rvf, ****error, ****data_clip, *****rvf_mod;
	float dt, dw, k_start, k_end, intercept, fm, amp, t0=0.04;

/********VPVSROU*********/

sf_file Fvp=NULL, Fvs=NULL, Frou=NULL;

Fvp = sf_input ( "vp" );
Fvs = sf_input ( "vs" );
Frou= sf_input ( "rou" );

float **vp, **vs, **rou;
float sita, radian;
float a,b,c;
float ****ref_old,****ref_new;

vp = alloc2float( 501 , 601 );
vs = alloc2float( 501 , 601 );
rou= alloc2float( 501 , 601 );

/************************/

	Fd = sf_input ( "data" );
	Fw = sf_input ( "wave" );
	Fo = sf_output( "rvf"  );

	dim = sf_filedims ( Fd , num_dim );
	/* dim=1, data[nt]*/
	/* dim=2, data[nsita][nt]*/
	/* dim=3, data[ninline][nsita][nt]*/
	/* dim=4, data[ncrossline][ninline][nsita][nt]*/

	if(!sf_getint( "flag" , &flag ))			sf_error( "Parameter flag needs to be speficied!\n" );
	/* flag=0 means that we need     input wavelet */
	/* flag=1 means that we need not input  wavelet, Genrated by the subrountine_ricker() */

	if( flag==0 ){
	if(!sf_getint( "wave_shift" , &wave_shift ))		sf_error( "Parameter wave_shift needs to be speficied!\n" );
	}
	if( flag==1 ){
	if(!sf_getfloat( "fm" , &fm ))				sf_error( "Parameter fm needs to be speficied!\n" );
	/* Parameter fm should be fixed! */
	if(!sf_getfloat( "amp", &amp))				sf_error( "Parameter amp needs to be speficied!\n" );
	/* Parameter amp should be fixed! */
	}
	if(!sf_getint  ( "nsita_start" , &nsita_start))		sf_error( "Parameter nsita_start      needs to be speficied!\n" );
	if(!sf_getint  ( "nsita" , &nsita))			sf_error( "Parameter nsita            needs to be speficied!\n" );
	if(!sf_getint  ( "ninline_start" , &ninline_start))	sf_error( "Parameter ninline_start    needs to be speficied!\n" );
	if(!sf_getint  ( "ninline" , &ninline))			sf_error( "Parameter ninline          needs to be speficied!\n" );
	if(!sf_getint  ( "ncrossline_start",&ncrossline_start))	sf_error( "Parameter ncrossline_start needs to be speficied!\n" );
	if(!sf_getint  ( "ncrossline" , &ncrossline))		sf_error( "Parameter ncrossline       needs to be speficied!\n" );

	if(!sf_getfloat( "k_start" , &k_start))			k_start		= -0.2;
	/* Optional parameter K-start is fixed to -0.2 */
	if(!sf_getfloat( "k_end" , &k_end))			k_end		= 0.2;
	/* Optional parameter K-end is fixed to 0.2 */
	if(!sf_getfloat( "intercept" , &intercept))		intercept	= 0.08;
	/* Optional parameter intercept is fixed to 0.08 */
	if(!sf_getint  ( "nk" , &nk))				nk		= 1000;
	/* Optional parameter nk is fixed to 1000 */
	if(!sf_getint  ( "nevent" , &nevent))			nevent		= 100;
	/* Optional parameter nevent is fixed to 100 */

/****0000**************Extract the axises from RSF data******************************/
	if( dim==1 )
	{
		at 		= sf_iaxa( Fd , 1 );
		nt 		= sf_n( at );
		dt		= sf_d( at );

		nsita_in	= 1;
		ninline_in	= 1;
		ncrossline_in	= 1;
	}
	else if( dim==2 )
	{
		at		= sf_iaxa( Fd , 1 );
		nt		= sf_n( at );
		dt		= sf_d( at );
	
		asita		= sf_iaxa( Fd , 2 );
		nsita_in	= sf_n( asita );

		ninline_in	= 1;
		ncrossline_in	= 1;
	}
	else if( dim==3 )
	{
		at		= sf_iaxa( Fd , 1 );
		nt		= sf_n( at );
		dt		= sf_d( at );
	
		asita		= sf_iaxa( Fd , 2 );
		nsita_in	= sf_n( asita );

		ainline		= sf_iaxa( Fd , 3 );
		ninline_in	= sf_n( ainline );

		ncrossline_in	= 1;
	}
	else if( dim==4 )
	{
		at		= sf_iaxa( Fd , 1 );
		nt		= sf_n( at );
		dt		= sf_d( at );
	
		asita		= sf_iaxa( Fd , 2 );
		nsita_in	= sf_n( asita );

		ainline		= sf_iaxa( Fd , 3 );
		ninline_in	= sf_n( ainline );

		acrossline	= sf_iaxa( Fd , 4 );
		ncrossline_in	= sf_n( acrossline );
	}
	else
	{
		sf_error( "Can't handle the problem with dimension higher than 4!\n " );
	}

	if( flag==0 )
	{
		awave 		= sf_iaxa( Fw , 1 );
		nwave 		= sf_n( awave );
		if( sf_d( awave )!=dt )
			sf_error( "Sampling of wavelet and the one of the data doesn't conincident!\n" );
	}
/****0000******************END*****************************************************************************************/
/****1111*******Check the Parameter*************************************************************************************/
	if( nsita_start>=nsita_in || nsita_start<0 )			sf_error( "Parameter nsita_start wrong!\n" );
	if( ninline_start>=ninline_in || ninline_start<0 )		sf_error( "Parameter ninline_start wrong!\n" );
	if( ncrossline_start>=ncrossline_in || ncrossline_start<0 )	sf_error( "Parameter ncrossline_start wrong!\n" );

	if( (nsita+nsita_start)>nsita_in )				nsita = nsita_in-nsita_start;
	if( (ninline+ninline_start)>ninline_in )			ninline = ninline_in-ninline_start;
	if( (ncrossline+ncrossline_start)>ncrossline_in)		ncrossline = ncrossline_in-ncrossline_start;
/****1111******END******************************************************************************************************/
/****2222**********Determine the parameters for FFT*********************************************************************/
	for( nfft=1 , i=1 ; i<16 ; i++ )
	{
		nfft = 2*nfft;
		if( nt<=nfft )
			break;
	}
	nfre = nfft/2+1;
	dw = 1./(dt*nfft);
/****2222*******************END*****************************************************************************************/
/****3333**************Output axises to the RVF.RSF data****************************************************************/

	aw = sf_maxa( nfre , 0. , dw );

	if( dim==1 )
	{
		sf_oaxa( Fo , at , 1 );
		sf_oaxa( Fo , aw , 2 );
	}
	else if( dim==2 )
	{
		sf_oaxa( Fo , at    , 1 );
		sf_oaxa( Fo , aw    , 2 );
		sf_setn( asita , nsita );
		sf_seto( asita , sf_o( asita )+nsita_start*sf_d( asita ) );
		sf_oaxa( Fo , asita , 3 );
	}
	else if( dim==3 )
	{
		sf_oaxa( Fo , at    , 1 );
		sf_oaxa( Fo , aw    , 2 );

		sf_setn( asita , nsita );
		sf_seto( asita , sf_o( asita )+nsita_start*sf_d( asita ) );
		sf_oaxa( Fo , asita , 3 );

		sf_setn( ainline , ninline );
		sf_seto( ainline , sf_o( ainline )+ninline_start*sf_d( ainline ) );
		sf_oaxa( Fo , ainline , 4 );
	}
	else if( dim==4 )
	{
		sf_oaxa( Fo , at    , 1 );
		sf_oaxa( Fo , aw    , 2 );

		sf_setn( asita , nsita );
		sf_seto( asita , sf_o( asita )+nsita_start*sf_d( asita ) );
		sf_oaxa( Fo , asita , 3 );

		sf_setn( ainline , ninline );
		sf_seto( ainline , sf_o( ainline )+ninline_start*sf_d( ainline ) );
		sf_oaxa( Fo , ainline , 4 );

		sf_setn( acrossline , ncrossline );
		sf_seto( acrossline , sf_o( acrossline )+ncrossline_start*sf_d( acrossline ) );
		sf_oaxa( Fo , acrossline , 5 );
	}
	else
	{
		sf_error( "Wrong dim!\n" );
	}


/****3333****************************END********************************************************************/

/********VpVSROU*********/
ref_old = sf_floatalloc4( 500 , nsita , 601 , ncrossline );
ref_new = sf_floatalloc4( 501 , nsita , 601 , ncrossline );

zero4float( ref_old , 500 , nsita , 601 , ncrossline );
zero4float( ref_new , 501 , nsita , 601 , ncrossline );
/************************/

	data		= sf_floatalloc4( nt , nsita_in , ninline_in , ncrossline_in );
	data_clip	= sf_floatalloc4( nt , nsita    , ninline    , ncrossline );
	rvf		= sf_floatalloc5( nt , nfre , nsita , ninline , ncrossline );
	rvf_mod		= sf_floatalloc5( nt , nsita , ninline , ncrossline , nfre );
	error		= sf_floatalloc4( nt , nsita , ninline , ncrossline );
	wavelet 	= sf_floatalloc ( nfft );

	zero4float( data      , nt , nsita_in , ninline_in , ncrossline_in );
	zero4float( data_clip , nt , nsita    , ninline    , ncrossline );
	zero5float( rvf , nt , nfre , nsita , ninline , ncrossline );
	zero4float( error , nt , nsita , ninline , ncrossline );
	zero1float( wavelet , nfft );

	if( flag==0 )
	{
		if( nwave<nfft )
		{
			sf_floatread( wavelet , nwave , Fw );
			lh_shift_array( wavelet , nfft , wave_shift );
		}
		else
		{
			sf_floatread( wavelet , nfft , Fw );
			lh_shift_array( wavelet , nfft , wave_shift );
		}
	}

	if( flag==1 )
	{
		lh_ricker( wavelet , amp , nfft , dt , fm , t0 );
		lh_shift_array( wavelet , nfft , (int)(t0/dt+0.5) );
	}

	sf_floatread( &data[0][0][0][0] , nt*nsita_in*ninline_in*ncrossline_in , Fd );

	for( icrossline=0 ; icrossline<ncrossline ; icrossline++ )
	for( iinline=0    ; iinline<ninline       ; iinline++ )
	for( isita=0      ; isita<nsita           ; isita++ )
	for( it=0         ; it<nt                 ; it++ )
		data_clip[icrossline][iinline][isita][it] = data[icrossline+ncrossline_start][iinline+ninline_start][isita+nsita_start][it];

/********VPVSROU*********/
sf_floatread( &vp[0][0] , 501*601 , Fvp );
sf_floatread( &vs[0][0] , 501*601 , Fvs );
sf_floatread( &rou[0][0], 501*601 , Frou);

for( icrossline=0 ; icrossline<ncrossline ; icrossline++ )
for( iinline=0    ; iinline<ninline       ; iinline++ )
for( isita=0, sita=sf_o( asita )+nsita_start*sf_d( asita )      ; isita<nsita           ; isita++ , sita+=sf_d(asita) )
for( it=1         ; it<nt                 ; it++ )
{
	radian = sita*SF_PI/180;

	if( vp[iinline][it]!=0 && vp[iinline][it-1]!=0 && vs[iinline][it]!=0 && vs[iinline][it-1]!=0 && rou[iinline][it]!=0 && rou[iinline][it-1]!=0 ){
	a = 0.5*(1.+tan(radian)*tan(radian));
	b = -4.*sin(radian)*sin(radian)*(vs[iinline][it]/2+vs[iinline][it-1]/2)*(vs[iinline][it]/2+vs[iinline][it-1]/2)/((vp[iinline][it]/2+vp[iinline][it-1]/2)*(vp[iinline][it]/2+vp[iinline][it-1]/2));
	c = 0.5*(1.+b);

	ref_old[icrossline][iinline][isita][it-1] = a*((vp[iinline][it]-vp[iinline][it-1])/(vp[iinline][it]/2+vp[iinline][it-1]/2))+b*((vs[iinline][it]-vs[iinline][it-1])/(vs[iinline][it]/2+vs[iinline][it-1]/2))+c*((rou[iinline][it]-rou[iinline][it-1])/(rou[iinline][it]/2+rou[iinline][it-1]/2));
	}
}
for( icrossline=0 ; icrossline<ncrossline ; icrossline++ )
for( iinline=0    ; iinline<ninline       ; iinline++ )
for( isita=0      ; isita<nsita           ; isita++ )
for( it=0         ; it<nt                 ; it++ )
{
	if( it!=nt-1 )
		ref_new[icrossline][iinline][isita][it] = ref_old[icrossline][iinline][isita][it];
	
}
lh_write_1d_float_bin( &ref_new[0][0][0][0] , ncrossline*ninline*nsita*nt , "ref_new.bin" );

/************************/

	//lh_mp_rvf_inputwave( data_clip , ncrossline , ninline , nsita , nt , dt , k_start , k_end , nk , intercept , nevent , nfft , rvf , error , wavelet );

/********VPVSROU*********/
lh_mp_rvf_inputwave_ref( data_clip , ref_new , ncrossline , ninline , nsita , nt , dt , 30. , k_start , k_end , nk , intercept  , nfft , rvf , error , wavelet );
/************************/


	sf_floatwrite( &rvf[0][0][0][0][0] , ncrossline*ninline*nsita*nfre*nt , Fo );
	for( icrossline=0 ; icrossline<ncrossline ; icrossline++ )
	for( iinline=0    ; iinline<ninline       ; iinline++ )
	for( isita=0      ; isita<nsita           ; isita++ )
	for( iw=0	  ; iw<nfre		  ; iw++ )
	for( it=0         ; it<nt                 ; it++ )
		rvf_mod[iw][icrossline][iinline][isita][it] = rvf[icrossline][iinline][isita][iw][it];

	lh_write_1d_float_bin( &rvf_mod[0][0][0][0][0] , ncrossline*ninline*nsita*nfre*nt , "rvf_mod.bin" );

	exit( 0 );
}
