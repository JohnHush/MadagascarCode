/* This file is automatically generated. DO NOT EDIT! */

#ifndef _lh_wavelet_h
#define _lh_wavelet_h


void lh_generate_w1w2_for_qinversion( float *ori_wavelet /* Original Wavelet */ ,
			      	    int nwave	       /* Length of Wavelet , should equal to the power(2,k))*/,
				    float dt	       /* Time Sampling Interval Of wavelet */,
				    float *Q_value     /* The inverse of estimated Q */,
				    int nposition      /* Output w1 position, must be less than the dimension of Q_value */,
				    float *w1	       /* Output w1 */ ,
				    float *w2	       /* Output w2 */ );
/*< Generate w1&w2 in Linearized Q inversion >*/


void lh_time_variant_wavelet( float *ori_wavelet /* Original Wavelet */,
			      int nwave		 /* Length of Wavelet , should equal to the power(2,k))*/,
			      float dt		 /* Time Interval of Wavelet */,
			      float *Q_value	 /* The inverse of Q value */,
			      int nposition	 /* Output wavelet position , must be less than nt */,
			      float *att_wavelet /* Output wavelet */ );
/*< Generate time variant wavelet while Q is variable with time >*/


void lh_time_variant_matrix( float *ori_wavelet	/* Original Wavelet */,
			     int wave_shift	/* Time lapse of the Wavelet */,
			     int nwave		/* Length of wavelet, should equal power(2,k) */,
			     float dt		/* Time Interval */,
			     float *Q_value	/* The Inverse Of Q value */,
			     int nt		/* The dimension of the Matrix */,
			     float **TVM	/* The generated Matrix */);
/*< Generate A Time Variant Matrix Using an Original Wavelet and a Q value vector >*/


void lh_wyh_Q_model( float *ori_wavelet		/* Original Wavelet */,
		     float *res_wavelet		/* Resulted Wavelet */,
		     float q			/* the constant Q in the layer */,
		     float t			/* two-way time in the layer */,
		     int nwave			/* length of the wavelet */,
		     float dt			/* sampling of the wavelet */,
		     float w0			/* the reference frequency */);
/*< Q model from Wang YangHua,2002,2006 >*/


void lh_ricker(	float *f 	/* array for wavelet */,
		float amp	/* Amplitude of the wavelet */,
		int n 		/* length of wavelet */,
		float dt	/* sampling of the wavelet */,
		float fm	/* main frequency of the wavelet */,
		float t0	/* the delayed time */);
/*< Ricker wavelet generator >*/


void lh_unit_ricker_nonzero_phase(
        float *wavelet  /* output wavelet */,
        int nt          /* intput length */,
        float dt        /* input time interval */,
        float fm        /* input main frequency */,
        float t0        /* input starting time (unit=second)*/,
        float phase     /* input rotation phase (unit=degree) */
        );
/*< Generate ricker wavelet with different phase by integrating Hilbert transform >*/


void lh_signal_envelope(
        float *a        /* original signal */,
        float *enve     /* envelope of the original signal */,
        int n           /* length of the signal */
        );
/*< compute the envelope of a specified signal >*/


void lh_signal_phase(
        float *a        /* input signal */,
        float *phase    /* output unwrapped phase */,
        float eps       /* parameter controls phase unwrapping */,
        int n           /* signal length */
        );
/*< Compute the unwrapping phase of a signal based on Hilbert transform >*/


void lh_signal_frequency(
        float *phase    /* input instantaneous phase computed from LH_SIGNAL_PHASE function */,
        float *freq     /* output instataneous frequency */,
        float dt        /* time interval of the signal */,
        int n           /* length of input */
        );
/*< Compute the instataneous frequency by Hilbert transform >*/


void lh_shift_array( float *x	/* array to be shifted */,
		     int n	/* size of the array */,
		     int step	/* step>0, shift to the left; step<0, shift to the right */);
/*< Shift the array >*/


int lh_powerof2( int x );
/*< compute the proper number 2^k<x<OUTPUT=2^(k+1) >*/


void lh_hilbert(    float *a    /* input signal*/,
                    float *b    /* output signal*/,
                    int n       /* signal length*/  );
/*< hilbert transform of a signal >*/


void lh_fft( float x[]	/* Real Part of the complex Array */,
        float y[]	/* Image part of the complex Array */,
        int n		/* The DImension */,
        int sign	/* Sign=1, Forward FFT, Sign=-1, IFFT */);
/*< fft >*/

#endif
