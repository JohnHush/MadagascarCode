/* This file is automatically generated. DO NOT EDIT! */

#ifndef _st_avf_h
#define _st_avf_h


void winft( float *data_in      /* Input wavelet */,
            float *data_out     /* Output wavelet */,
            int nt              /* length of wavelet */,
            int nfft            /* Length of FFT */,
            float dt            /* Time Interval of TIME DOMAIN */,
            float fre           /* Window Parameter, central frequency */,
            float var           /* Window Parameter, Variance of the Window */);
/*< Windowed FFT & IFFT >*/


void wavelet_filtered(	float *wavelet		/* The Original Wavelet */,
			float *wave_filtered	/* The Filtered Wavelet */,
			int n			/* Length Of Wavelet */,
			float dt		/* Time Interval of the Wavelets */,
			float fre		/* The Window Parameter, central frequency*/,
			float variance		/* The Window Parameter, Variance of the window */);
/*< Compute the filtered wavelet using Fourier transform >*/


void cross_correlation_ft( float *x , int nx , float *y , int ny , float *cross , int nc );
/*< Compute the correlation between array x and array y >*/


void move_array( float *x , int n , int step );
/*< Move array x in step >*/


void move_spec( float *x , float *y , int n , int step );
/*< Move array x & y in step >*/


void gst_avf( float *gt , float **s_real , float **s_imag , float tao , int n , float min_fre , int nw );
/*< Compute time-frequency using S-transform >*/


void gst_fre( float *gt , float *s_real , float *s_imag , float tao , int n , int nw );
/*< Compute time-frequency using S-Transform subordinate programme >*/


void gst_fre_plus( float *gt , float *s_real , float *s_imag , float tao , int k , int nw );
/*< Compute time-frequency using S-Transform subordinate programme >*/


void fft( float x[], float y[], int n,int sign);
/*< Fourier Transform >*/


void gist_filtered(	float *gt	/* The filtered Trace */,
			float **s_real	/* The Input real part of ST Spectrum */,
			float **s_imag	/* The Input Imag part of ST Spectrum */,
			int n		/* Length of Spectrum */,
			float dt	/* Time Interval of the data */,
			float fre	/* The Window Parameter, Central frequency */,
			float variance	/* The Window Parameter, Variance of the Window */);
/*< Compute frequency-depended seismic data using inverse S-Transform >*/

#endif
