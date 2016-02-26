//#include <rsf.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "lh_wavelet.h"
#include "su_alloc.h"

#ifndef SF_PI

#define SF_PI 3.1415926

#endif

void lh_generate_w1w2_for_qinversion( float *ori_wavelet /* Original Wavelet */ ,
			      	    int nwave	       /* Length of Wavelet , should equal to the power(2,k))*/,
				    float dt	       /* Time Sampling Interval Of wavelet */,
				    float *Q_value     /* The inverse of estimated Q */,
				    int nposition      /* Output w1 position, must be less than the dimension of Q_value */,
				    float *w1	       /* Output w1 */ ,
				    float *w2	       /* Output w2 */ )
/*< Generate w1&w2 in Linearized Q inversion >*/
{
	float *real_tmp, *imag_tmp, dw, w, sigma0=0., *real_tmp2, *imag_tmp2;
	int i,iw, it;

	real_tmp = alloc1float( nwave );
	imag_tmp = alloc1float( nwave );
	real_tmp2= alloc1float( nwave );
	imag_tmp2= alloc1float( nwave );

	zero1float( real_tmp , nwave );
	zero1float( imag_tmp , nwave );
	zero1float( real_tmp2, nwave );
	zero1float( imag_tmp2, nwave );

	for( i=0 ; i<nwave ; i++ )
	{
		real_tmp[i] = ori_wavelet[i];
		real_tmp2[i]= ori_wavelet[i];
	}

	dw = 1.0/(dt*nwave);

	for( i=0 ; i<nposition ; i++ )
		sigma0+=Q_value[i];

	lh_fft( real_tmp , imag_tmp , nwave , 1 );
	lh_fft( real_tmp2, imag_tmp2, nwave , 1 );

	for( iw=0 ; iw<=nwave/2 ; iw++ )
	{
		w            = dw*iw;
		real_tmp[iw] = real_tmp[iw]*exp(-SF_PI*dt*w*sigma0)*(1.+SF_PI*dt*sigma0*w);
		imag_tmp[iw] = imag_tmp[iw]*exp(-SF_PI*dt*w*sigma0)*(1.+SF_PI*dt*sigma0*w);

		real_tmp2[iw]= real_tmp2[iw]*exp(-SF_PI*dt*w*sigma0)*(-SF_PI*dt*w);
		imag_tmp2[iw]= imag_tmp2[iw]*exp(-SF_PI*dt*w*sigma0)*(-SF_PI*dt*w);
	}

	for( iw=1 ; iw<nwave/2 ; iw++ )
	{
		real_tmp[nwave-iw] = real_tmp[iw];
		imag_tmp[nwave-iw] = -imag_tmp[iw];

		real_tmp2[nwave-iw] = real_tmp2[iw];
		imag_tmp2[nwave-iw] = -imag_tmp2[iw];
	}

	lh_fft( real_tmp , imag_tmp , nwave , -1 );
	lh_fft( real_tmp2, imag_tmp2, nwave , -1 );

	for( i=0 ; i<nwave ; i++ )
	{
		w1[i] = real_tmp[i];
		w2[i] = real_tmp2[i];
	}

	free1float( real_tmp );
	free1float( imag_tmp );
	free1float( real_tmp2);
	free1float( imag_tmp2);
}

void lh_time_variant_wavelet( float *ori_wavelet /* Original Wavelet */,
			      int nwave		 /* Length of Wavelet , should equal to the power(2,k))*/,
			      float dt		 /* Time Interval of Wavelet */,
			      float *Q_value	 /* The inverse of Q value */,
			      int nposition	 /* Output wavelet position , must be less than nt */,
			      float *att_wavelet /* Output wavelet */ )
/*< Generate time variant wavelet while Q is variable with time >*/
{
	float *real_tmp, *imag_tmp, dw, w;
	int i,iw, it;

	real_tmp = alloc1float( nwave );
	imag_tmp = alloc1float( nwave );

	zero1float( real_tmp , nwave );
	zero1float( imag_tmp , nwave );

	for( i=0 ; i<nwave ; i++ )
		real_tmp[i] = ori_wavelet[i];

	dw = 1.0/(dt*nwave);

	for( it=0 ; it<nposition ; it++ )
	{
		lh_fft( real_tmp , imag_tmp , nwave , 1 );

		for( iw=0 ; iw<=nwave/2 ; iw++ )
		{
			w            = dw*iw;
			real_tmp[iw] = real_tmp[iw]*exp(-SF_PI*dt*w*Q_value[it]);
			imag_tmp[iw] = imag_tmp[iw]*exp(-SF_PI*dt*w*Q_value[it]);
		}

		for( iw=1 ; iw<nwave/2 ; iw++ )
		{
			real_tmp[nwave-iw] = real_tmp[iw];
			imag_tmp[nwave-iw] = -imag_tmp[iw];
		}
		lh_fft( real_tmp , imag_tmp , nwave , -1 );

		zero1float( imag_tmp , nwave );
	}

	for( i=0 ; i<nwave ; i++ )
		att_wavelet[i] = real_tmp[i];

	free1float( real_tmp );
	free1float( imag_tmp );
}

void lh_time_variant_matrix( float *ori_wavelet	/* Original Wavelet */,
			     int wave_shift	/* Time lapse of the Wavelet */,
			     int nwave		/* Length of wavelet, should equal power(2,k) */,
			     float dt		/* Time Interval */,
			     float *Q_value	/* The Inverse Of Q value */,
			     int nt		/* The dimension of the Matrix */,
			     float **TVM	/* The generated Matrix */)
/*< Generate A Time Variant Matrix Using an Original Wavelet and a Q value vector >*/
{
	float *real_tmp, *imag_tmp, dw, w, *sigma;
	int i,iw, it, j;

	real_tmp = alloc1float( nwave );
	imag_tmp = alloc1float( nwave );
	sigma	 = alloc1float( nt );

    zero1float( sigma , nt );

	for( i=0 ; i<nt ; i++ )
	for( j=0 ; j<i  ; j++ )
		sigma[i]+=Q_value[j];

	dw = 1.0/(dt*nwave);

	for( it=0 ; it<nt ; it++ )
	{
		zero1float( real_tmp , nwave );
		zero1float( imag_tmp , nwave );

		for( i=0 ; i<nwave ; i++ )
			real_tmp[i] = ori_wavelet[i];

		lh_fft( real_tmp , imag_tmp , nwave , 1 );

		for( iw=0 ; iw<=nwave/2 ; iw++ )
		{
			w            = dw*iw;
			real_tmp[iw] = real_tmp[iw]*exp(-SF_PI*dt*w*sigma[it]);
			imag_tmp[iw] = imag_tmp[iw]*exp(-SF_PI*dt*w*sigma[it]);
		}

		for( iw=1 ; iw<nwave/2 ; iw++ )
		{
			real_tmp[nwave-iw] = real_tmp[iw];
			imag_tmp[nwave-iw] = -imag_tmp[iw];
		}
		lh_fft( real_tmp , imag_tmp , nwave , -1 );
		
		for( i=0 ; i<nt ; i++ )
		if( wave_shift+i-it>=0 && wave_shift+i-it<nwave )
			TVM[i][it] = real_tmp[wave_shift+i-it];
	}
	free1float( real_tmp );
	free1float( imag_tmp );
	free1float( sigma );
}

void lh_wyh_Q_model( float *ori_wavelet		/* Original Wavelet */,
		     float *res_wavelet		/* Resulted Wavelet */,
		     float q			/* the constant Q in the layer */,
		     float t			/* two-way time in the layer */,
		     int nwave			/* length of the wavelet */,
		     float dt			/* sampling of the wavelet */,
		     float w0			/* the reference frequency */)
/*< Q model from Wang YangHua,2002,2006 >*/
{
	float *real_tmp, *imag_tmp, dw;
	float complex *Att, gamma, w;
	int i,iw;

	real_tmp = alloc1float( nwave );
	imag_tmp = alloc1float( nwave );

	zero1float( real_tmp , nwave );
	zero1float( imag_tmp , nwave );

	Att = (float complex *)malloc( nwave*sizeof(float complex) );

	for( i=0 ; i<nwave ; i++ )
	{
		real_tmp[i] = ori_wavelet[i];
		Att[i] = 0.+0.*I;
	}
	lh_fft( real_tmp , imag_tmp , nwave , 1 );
	gamma = -1.0/(SF_PI*q);
	dw = 1.0/(dt*nwave);

	for( iw=0 ; iw<=nwave/2 ; iw++ )
	{
		w            = dw*iw;
		Att[iw]      = exp(-pow((w/w0),-gamma)*((SF_PI*w*t)/(q)) );
		real_tmp[iw] = crealf(real_tmp[iw]*Att[iw]);
		imag_tmp[iw] = crealf(imag_tmp[iw]*Att[iw]);
	}
	for( iw=1 ; iw<nwave/2 ; iw++ )
	{
		real_tmp[nwave-iw] = real_tmp[iw];
		imag_tmp[nwave-iw] = -imag_tmp[iw];
	}
	lh_fft( real_tmp , imag_tmp , nwave , -1 );
	for( i=0 ; i<nwave; i++ )
		res_wavelet[i] = real_tmp[i];
}

void lh_ricker(	float *f 	/* array for wavelet */,
		float amp	/* Amplitude of the wavelet */,
		int n 		/* length of wavelet */,
		float dt	/* sampling of the wavelet */,
		float fm	/* main frequency of the wavelet */,
		float t0	/* the delayed time */)
/*< Ricker wavelet generator >*/
{
        float t;
        int i;
        for( i=0 , t=-t0 ; i<n ; i++ , t+=dt )
                f[i] = amp*(1.0-2.0*(SF_PI*fm*t*SF_PI*fm*t))*exp(-SF_PI*fm*t*SF_PI*fm*t);

}

void lh_unit_ricker_nonzero_phase(
        float *wavelet  /* output wavelet */,
        int nt          /* intput length */,
        float dt        /* input time interval */,
        float fm        /* input main frequency */,
        float t0        /* input starting time (unit=second)*/,
        float phase     /* input rotation phase (unit=degree) */
        )
/*< Generate ricker wavelet with different phase by integrating Hilbert transform >*/
{
    float *ricker, *hricker;
    float amp=1., power=0.;
    int i;

    ricker = alloc1float( nt );
    hricker= alloc1float( nt );

    zero1float( ricker , nt );
    zero1float( hricker, nt );

    lh_ricker( ricker , amp , nt , dt , fm , t0 );
    lh_hilbert( ricker , hricker , nt );

    for( i=0 ; i<nt ; i++ )
    {
        wavelet[i] = crealf( cexpf( (phase*SF_PI/180.)*I )*( ricker[i]+hricker[i]*I ) );
        power += wavelet[i]*wavelet[i];
    }
    for( i=0 ; i<nt ; i++ )
        wavelet[i] /= sqrt(power);

    free1float( ricker );
    free1float( hricker );
}

void lh_signal_envelope(
        float *a        /* original signal */,
        float *enve     /* envelope of the original signal */,
        int n           /* length of the signal */
        )
/*< compute the envelope of a specified signal >*/
{
    float *ha;
    int i;

    ha = alloc1float( n );
    zero1float( ha , n );

    lh_hilbert( a , ha , n );

    for( i=0 ; i<n ; i++ )
        enve[i] = sqrt( a[i]*a[i] + ha[i]*ha[i] );

    free1float( ha );
}

void lh_signal_phase(
        float *a        /* input signal */,
        float *phase    /* output unwrapped phase */,
        float eps       /* parameter controls phase unwrapping */,
        int n           /* signal length */
        )
/*< Compute the unwrapping phase of a signal based on Hilbert transform >*/
{
    float *ha;
    int i,j;

    ha = alloc1float( n );

    zero1float( ha , n );

    lh_hilbert( a , ha , n );

    for( i=0 ; i<n ; i++ )
        phase[i] = atan2( ha[i] , a[i] );

    for( i=1 ; i<n ; i++ )
    {
        if( ( phase[i]-phase[i-1] )>(2.*SF_PI-eps) )
        {
            for( j=i ; j<n ; j++ )
                phase[j] -= 2*SF_PI;
        }
        if( ( phase[i]-phase[i-1] )<(-2.*SF_PI+eps) )
        {
            for( j=i ; j<n ; j++ )
                phase[j] += 2*SF_PI;
        }
    }
    free1float( ha );
}

void lh_signal_frequency(
        float *phase    /* input instantaneous phase computed from LH_SIGNAL_PHASE function */,
        float *freq     /* output instataneous frequency */,
        float dt        /* time interval of the signal */,
        int n           /* length of input */
        )
/*< Compute the instataneous frequency by Hilbert transform >*/
{
    int i;

    for( i=1 ; i<n ; i++ )
        freq[i] = (phase[i]-phase[i-1])/(dt*2*SF_PI);

    freq[0] = freq[1];
}

void lh_shift_array( float *x	/* array to be shifted */,
		     int n	/* size of the array */,
		     int step	/* step>0, shift to the left; step<0, shift to the right */)
/*< Shift the array >*/
{
        int i;
        float *tmp;
        tmp = (float *)malloc(n*sizeof(float));

        step = step%n;
        if( step<0 )
                step = step+n;

        for( i=0;i<n;i++ )
        {
                if( i<n-step )
                        tmp[i] = x[i+step];
                else
                        tmp[i] = x[i-n+step];
        }
        for( i=0;i<n;i++ )
                x[i] = tmp[i];
        free( tmp );
}

int lh_powerof2( int x )
/*< compute the proper number 2^k<x<OUTPUT=2^(k+1) >*/
{
    int i,j;
    for( j=1,i=0 ; i<16 ; i++ )
    {
        j=2*j;
        if( j>=x )
            break;
    }
    return j;
}

void lh_hilbert(    float *a    /* input signal*/,
                    float *b    /* output signal*/,
                    int n       /* signal length*/  )
/*< hilbert transform of a signal >*/
{
    int i , nfft;

    nfft=1;
    for( i=1 ; i<16 ; i++ )
    {
        nfft=nfft*2;
        if( (n)<=nfft )
            break;
    }

    float *rm, *im , *hrm, *him ;

    rm = alloc1float( nfft );
    im = alloc1float( nfft );
    hrm = alloc1float( nfft );
    him = alloc1float( nfft );

    zero1float( rm , nfft );
    zero1float( im , nfft );
    zero1float( hrm , nfft );
    zero1float( him , nfft );

    for( i=0 ; i<nfft ; i++ )
    {
        if( i<n )
            rm[i] = a[i];
    }
    lh_fft( rm , im , nfft , 1 );

    for( i=0 ; i<nfft ; i++ )
    {
        if( i>=0 && i<(nfft/2)-1 )
        {
            hrm[i] = -im[i];
            him[i] = rm[i];
        }
        else if ( i==(nfft/2)-1 )
        {
            hrm[i] = 0.;
            him[i] = 0.;
        }
        else
        {
            hrm[i] = im[i];
            him[i] = -rm[i];
        }
    }
    lh_fft( hrm , him , nfft , -1 );

    for( i=0 ; i<n ; i++ )
        b[i] = hrm[i];

    free1float( rm );
    free1float( im );
    free1float( hrm );
    free1float( him );

}

void lh_fft( float x[]	/* Real Part of the complex Array */,
        float y[]	/* Image part of the complex Array */,
        int n		/* The DImension */,
        int sign	/* Sign=1, Forward FFT, Sign=-1, IFFT */)
/*< fft >*/
{
        int i,j,k,l,m,n1,n2;
        float c,c1,e,s,s1,t,tr,ti;
        for(j=1,i=1;i<16;i++)
        {
                m=i;
                j=2*j;
                if(j==n)
                        break;
        }
        n1=n-1;
        for(j=0,i=0;i<n1;i++)
        {
                if(i<j)
                {
                        tr=x[j];
                        ti=y[j];
                        x[j]=x[i];
                        y[j]=y[i];
                        x[i]=tr;
                        y[i]=ti;
                }
                k=n/2;
                while(k<(j+1))
                {
                        j=j-k;
                        k=k/2;
                }
                j=j+k;
        }
        n1=1;
        for(l=1;l<=m;l++)
        {
                n1=2*n1;
                n2=n1/2;
                e=3.14159265359/n2;
                c=1.0;
                s=0.0;
                c1=cos(e);
                s1=-sign*sin(e);
                for(j=0;j<n2;j++)
                {
                        for(i=j;i<n;i+=n1)
                        {
                                k=i+n2;
                                tr=c*x[k]-s*y[k];
                                ti=c*y[k]+s*x[k];
                                x[k]=x[i]-tr;
                                y[k]=y[i]-ti;
                                x[i]=x[i]+tr;
                                y[i]=y[i]+ti;
                        }
                        t=c;
                        c=c*c1-s*s1;
                        s=t*s1+s*c1;
                }
        }
        if(sign==-1)
        {
                for(i=0;i<n;i++)
                {
                        x[i]/=n;
                        y[i]/=n;
                }
        }
}

