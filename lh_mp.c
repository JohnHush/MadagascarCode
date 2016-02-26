#include <rsf.h>
#include "lh_bayesian.h"
#include "lh_wavelet.h"
#include "su_alloc.h"
#include "lh_readwrite.h"
#include "lh_mp.h"

struct wave_project
{
        float *wavelet;
        float k;
        float intercept;
        float energy;
        float *correlation;
};

int lh_index_max(   float *a    /* input array */,
                    int n       /* length */)
/*< return the index with max value >*/
{
    int i, index=0;
    float val;

    val = a[0];
    for( i=1 ; i<n ; i++ )
        if( a[i]>val )
        {
            val = a[i];
            index = i;
        }
    return index;
}

int lh_index_maxabs(float *a    /* input array */,
                    int n       /* length */)
/*< return the index with max value >*/
{
    int i, index=0;
    float val;

    val = fabs(a[0]);
    for( i=1 ; i<n ; i++ )
        if( fabs(a[i])>val )
        {
            val = fabs(a[i]);
            index = i;
        }
    return index;
}

float lh_float_iproduct( float *a       /* The array A[n] */,
			float *b   /* The array B[n] */,
			int n      /* The dimension */)
/*< Return the inner product of <a,b> >*/
{
	int i;
	float sum=0.;

	for( i=0 ; i<n ; i++ )
		sum+=a[i]*b[i];

	return sum;
}

int lh_max_abs(	float *a	/* The array */,
		int n		/* The dimension */ )
/*< Return the axes of number in array with maximum absolute value >*/
{
	int i;
	int max_ax=0;
	float max=fabs(a[0]);

	for( i=0 ; i<n ; i++ )
	{
		if( max<fabs(a[i]) )
		{
			max = fabs(a[i]);
			max_ax = i;
		}
	}
	return max_ax;
}

float lh_butterworth_filter( float frequency		/* Frequency coefficient */,
			     float cut_frequency	/* Cut frequency */,
			     int step			/* The step of the filter */)
/*< Butterworth Filter >*/
{
        float hw;
        hw = 1.0/(1.0+pow((frequency/cut_frequency),2.0*step));
        return hw;
}

void lh_cross_correlation_ft( float *x		/* Array x[i] */,
			      int nx		/* Dimention of x */,
			      float *y		/* Array y[i] */,
			      int ny		/* Dimension of y */,
			      float *cross	/* Cross correlation cross<x,y> */,
			      int nc		/* Dimension of cross */)
/*< Cross correlation of x and y, CROSS<X,Y> >*/
{
	int i , fft_value;
	int n_max =nx ;

	if( n_max<ny )
		n_max = ny;

	fft_value=1;
	for( i=1 ; i<16 ; i++ )
	{
		fft_value=fft_value*2;
		if( (n_max)<=fft_value )
			break;
	}
	float *rm_x , *im_x;
	float *rm_y , *im_y;
	float *rm_cross , *im_cross;

	rm_x     = (float *)calloc(fft_value , sizeof(float));
	im_x     = (float *)calloc(fft_value , sizeof(float));
	rm_y     = (float *)calloc(fft_value , sizeof(float));
	im_y     = (float *)calloc(fft_value , sizeof(float));
	rm_cross     = (float *)calloc(fft_value , sizeof(float));
	im_cross     = (float *)calloc(fft_value , sizeof(float));

	for( i=0 ; i<nc ; i++ )
		cross[i] = 0;
	for( i=0 ; i<nx ; i++ )
		rm_x[i] = x[i];
	for( i=0 ; i<ny ; i++ )
		rm_y[i] = y[i];

	lh_fft( rm_x , im_x , fft_value , 1 );
	lh_fft( rm_y , im_y , fft_value , 1 );
	for( i=0 ; i<fft_value ; i++ )
		im_y[i] = -im_y[i];
	for( i=0 ; i<fft_value ; i++ )
	{
		rm_cross[i] = crealf( (rm_x[i] + im_x[i]*I)*(rm_y[i] + im_y[i]*I) );
		im_cross[i] = cimagf( (rm_x[i] + im_x[i]*I)*(rm_y[i] + im_y[i]*I) );
	}

	lh_fft( rm_cross , im_cross , fft_value , -1 );
	for( i=0 ; i<nc ; i++ )
	{
		cross[i] = rm_cross[i];
	}
	free( rm_x );
	free( im_x );
	free( rm_y );
	free( im_y );
	free( rm_cross );
	free( im_cross );
}

void lh_wave_estimation( float ****profile		/* Input profile */,
			 int ny				/* Number of Offset in crossline direction */,
			 int nx				/* Number of Offset in inline direction */,
			 int nsita			/* Number of angle in gather */,
			 int nt				/* Number of Time Point */,
			 float *wavelet_estimated	/* Output wavelet with length=nt */,
			 int fft_length			/* The length of the fft and the return wavelet */)
/*< Using Profile to estimate the wavelet based on white-noise assumption >*/
{
	float *trace_rm , *trace_im , *trace_rm_sum , *trace_im_sum;
	int j, iy , ix , isita;

	trace_rm         = alloc1float( fft_length );
	trace_im         = alloc1float( fft_length );
	trace_rm_sum     = alloc1float( fft_length );
	trace_im_sum     = alloc1float( fft_length );

	zero1float( trace_rm , fft_length );
	zero1float( trace_im , fft_length );
	zero1float( trace_rm_sum , fft_length );
	zero1float( trace_im_sum , fft_length );

	for( iy=0    ; iy<ny       ; iy++ )
	for( ix=0    ; ix<nx       ; ix++ )
	for( isita=0 ; isita<nsita ; isita++ )
	{
		zero1float( trace_rm , fft_length );
		zero1float( trace_im , fft_length );
	
		for( j=0 ; j<nt ; j++ )
		{
			trace_rm[j] = profile[iy][ix][isita][j];
		}
		lh_fft( trace_rm , trace_im , fft_length , 1 );
		for( j=0 ; j<fft_length ; j++ )
		{
			trace_rm_sum[j] = trace_rm_sum[j]+trace_rm[j];
			trace_im_sum[j] = trace_im_sum[j]+trace_im[j];
		}
	}
	for( j=0 ; j<fft_length ; j++ )
	{
		trace_rm_sum[j] = trace_rm_sum[j]/(nx*ny*nsita);
		trace_im_sum[j] = trace_im_sum[j]/(nx*ny*nsita);
	}
	for( j=0 ; j<fft_length ; j++ )
		wavelet_estimated[j] = sqrt( trace_rm_sum[j]*trace_rm_sum[j] + trace_im_sum[j]*trace_im_sum[j] );

// ********Filter**The***AMP**SPECTRUM******************* //
	zero1float( trace_im , fft_length );
	lh_fft( wavelet_estimated , trace_im , fft_length , 1 );

	for( j=0 ; j<=fft_length/2 ; j++ )
	{
		wavelet_estimated[j] = wavelet_estimated[j]*lh_butterworth_filter( (float)j , 50. , 2 );
		trace_im[j]          =          trace_im[j]*lh_butterworth_filter( (float)j , 50. , 2 );
		
		if( j!=0 &&j!=fft_length/2 )
		{
			wavelet_estimated[fft_length-j] = wavelet_estimated[j];
			trace_im[fft_length-j] = -trace_im[j];
		}
	}
	lh_fft( wavelet_estimated , trace_im , fft_length , -1 );
// ********Filter***END********************************** //

	zero1float( trace_im , fft_length );
	lh_fft( wavelet_estimated , trace_im , fft_length , -1 );

	free1float( trace_rm );
	free1float( trace_im );
	free1float( trace_rm_sum );
	free1float( trace_im_sum );

}

void lh_mp_rvf_inputwave( float ****data	/* 4D Input data[NY][NX][NSITA][NT]*/,
			  int ny		/* Input parameter, NY */,
			  int nx		/* Input parameter, NX */,
			  int nsita		/* Input parameter, NSITA */,
			  int nt		/* Input parameter, NT */,
		  	  float dt		/* Input parameter, DT */,
			  float k_start		/* Scan parameter, start slope */,
			  float k_end		/* Scan parameter, end slope */,
			  int nk		/* Scan parameter, Number of scan slope */,
			  float intercept	/* Scan parameter, value of intercept */,
			  int nevent		/* Scan parameter, Number of events */,
			  int nfft		/* Derived from NT,nfft>=nt  */,
			  float *****rvf	/* Output data, rvf[NY][NX][NSITA][NW][NT] */,
			  float ****error	/* Output data, error[NY][NX][NSITA][NT] */,
			  float *wavelet	/* Input data, wavelet[NFFT] */)
/*< Reflectivity versus Frequency extraction from seismic trace >*/
{
	int ik , iw , it , ievent , ix , iy , isita;
	int max_time=0 , max_ik=0 , nfre=nfft/2+1;
	float *MP_RVF_trace , *MP_RVF_wavelet;
	float k_value , fre_value , max_k=k_start , max_correlation=0. , max_ref_ratio;
	float dk=(k_end-k_start)/nk , dfre;
	float *MP_RVF_rvf , *MP_RVF_wave_rm , *MP_RVF_wave_im;
	struct wave_project *wave;

	dfre = 1.0/(dt*nfft);
	wave = (struct wave_project *)malloc( nk*sizeof( struct wave_project ));
	for( ik=0 ; ik<nk ; ik++ )
	{
		wave[ik].wavelet     = alloc1float( nfft );
		wave[ik].correlation = alloc1float( nfft );
		wave[ik].k           = 0.;
		wave[ik].intercept   = 0.;
		wave[ik].energy      = 0.;

		zero1float( wave[ik].wavelet     , nfft );
		zero1float( wave[ik].correlation , nfft );
	}
	
	MP_RVF_trace   = alloc1float( nfft );
	MP_RVF_wavelet = alloc1float( nfft );
	MP_RVF_rvf     = alloc1float( nfft );
	MP_RVF_wave_rm = alloc1float( nfft );
	MP_RVF_wave_im = alloc1float( nfft );

	zero1float( MP_RVF_trace   , nfft );
	zero1float( MP_RVF_wavelet , nfft );
	zero1float( MP_RVF_rvf     , nfft );
	zero1float( MP_RVF_wave_rm , nfft );
	zero1float( MP_RVF_wave_im , nfft );

/*******000000wavelet estimation*****************************************************************/
	
	for( it=0 ; it<nfft ; it++ )
		MP_RVF_wavelet[it] = wavelet[it];

/*********00000END*******************************************************************************/
/******1111GENERATE THE ATOMS LIBRARY************************************************************/
	for( ik=0,k_value=k_start ; ik<nk ; ik++, k_value+= dk )
	{
		wave[ik].k         = k_value;
		wave[ik].intercept = intercept;
		zero1float( MP_RVF_rvf     , nfft );
		zero1float( MP_RVF_wave_rm , nfft );
		zero1float( MP_RVF_wave_im , nfft );

		for( iw=0,fre_value=0. ; iw<nfre ; iw++,fre_value+= dfre )
		{
			if( iw==0 )
				MP_RVF_rvf[0]             = 0;
			else if( iw>0 && iw<nfre-1 )
			{
				MP_RVF_rvf[iw]            = intercept + k_value*log(fre_value);
				MP_RVF_rvf[nfft-iw] = MP_RVF_rvf[iw];
			}
			else
				MP_RVF_rvf[iw]            = intercept + k_value*log(fre_value);
		}

		for( it=0 ; it<nfft ; it++ )
		{
			MP_RVF_wave_rm[it] = MP_RVF_wavelet[it];
		}
		lh_fft( MP_RVF_wave_rm , MP_RVF_wave_im , nfft , 1 );
		for( iw=0 ; iw<nfft ; iw++ )
		{
			MP_RVF_wave_rm[iw] = MP_RVF_wave_rm[iw]*MP_RVF_rvf[iw];
			MP_RVF_wave_im[iw] = MP_RVF_wave_im[iw]*MP_RVF_rvf[iw];
		}
		lh_fft( MP_RVF_wave_rm , MP_RVF_wave_im , nfft , -1 );
		for( it=0 ; it<nfft ; it++ )
		{
			wave[ik].wavelet[it] = MP_RVF_wave_rm[it];
		}
		for( it=0 ; it<nfft ; it++ )
		{
			wave[ik].energy+= wave[ik].wavelet[it]*wave[ik].wavelet[it];
		}
		if( wave[ik].energy>0.0000001 )
		{
			for( it=0 ; it<nfft ; it++ )
				wave[ik].wavelet[it] = wave[ik].wavelet[it]/sqrt(wave[ik].energy);
		}
	}
/*******1111END************************************************************************************/
/*******2222**Use the generated atoms library to accomplish Matching Pursuit Method****************/
	for( iy=0 ; iy<ny ; iy++ )
	for( ix=0 ; ix<nx ; ix++ )
	for( isita=0 ; isita<nsita ; isita++ )
	{
printf( "icrossline=%d\tiinline=%d\tisita=%d\n" , iy , ix , isita );
		for( it=0 ; it<nt ; it++ )
		{
			MP_RVF_trace[it] = data[iy][ix][isita][it];
		}
		for( ievent=0 ; ievent<nevent ; ievent++ )
		{
			max_time=0;
			max_ik  =0;
			max_k = k_value;
			max_correlation = 0;
			for( ik=0,k_value=k_start ; ik<nk ; ik++, k_value+= dk )
			{
				lh_cross_correlation_ft( MP_RVF_trace , nfft , wave[ik].wavelet , nfft , wave[ik].correlation , nfft );
				for(  it=0 ; it<nfft ; it++ )
				{
					if( fabs(max_correlation)<fabs(wave[ik].correlation[it]) )
					{
						max_correlation = wave[ik].correlation[it];
						max_time        = it;
						max_k           = k_value;
						max_ik          = ik;
						max_ref_ratio   = wave[ik].correlation[it]/sqrt(wave[ik].energy);
					}
				}
			
			}
			lh_shift_array( wave[max_ik].wavelet , nfft , -max_time );
			for( it=0 ; it<nfft ; it++ )
			{
				MP_RVF_trace[it] = MP_RVF_trace[it]-max_correlation*wave[max_ik].wavelet[it];
			}
			lh_shift_array( wave[max_ik].wavelet , nfft , max_time );
			for( iw=1,fre_value=dfre ; iw<nfre ; iw++,fre_value+=dfre )
			{
				if( max_time<nt )
					rvf[iy][ix][isita][iw][max_time]+= (intercept + max_k*log(fre_value))*max_ref_ratio;
			}
		}
		for( it=0 ; it<nt ; it++ )
			error[iy][ix][isita][it] = MP_RVF_trace[it];
	}
/*******2222END************************************************************************************/

	free1float( MP_RVF_trace );
	free1float( MP_RVF_wavelet );
	free1float( MP_RVF_rvf );
	free1float( MP_RVF_wave_rm );
	free1float( MP_RVF_wave_im );
	free( wave );
}

void lh_mp_rvf_inputwave_ref(	float ****data	/* 4D Input data[NY][NX][NSITA][NT]*/,
				float ****ref	/* 4D Input ref[NY][NX][NSITA][NT]*/,
			  	int ny		/* Input parameter, NY */,
			  	int nx		/* Input parameter, NX */,
				int nsita	/* Input parameter, NSITA */,
				int nt		/* Input parameter, NT */,
				float dt	/* Input parameter, DT */,
				float ave_fre	/* Input parameter, Average Frequency */,
				float k_start	/* Scan parameter, start slope */,
				float k_end	/* Scan parameter, end slope */,
				int nk		/* Scan parameter, Number of scan slope */,
				float intercept	/* Scan parameter, value of intercept */,
				int nfft	/* Derived from NT,nfft>=nt  */,
				float *****rvf	/* Output data, rvf[NY][NX][NSITA][NW][NT] */,
				float ****error	/* Output data, error[NY][NX][NSITA][NT] */,
			  	float *wavelet	/* Input data, wavelet[NFFT] */)
/*< Reflectivity versus Frequency extraction from seismic trace >*/
{
	int ik , iw , it , ievent , ix , iy , isita;
	int max_time=0 , max_ik=0 , nfre=nfft/2+1;
	float *MP_RVF_trace , *MP_RVF_wavelet;
	float k_value , fre_value , max_k=k_start , max_ref_ratio;
	float dk=(k_end-k_start)/nk , dfre;
	float *MP_RVF_rvf , *MP_RVF_wave_rm , *MP_RVF_wave_im, *MP_RVF_ref;
	struct wave_project *wave;
	float *innerproduct;

	dfre = 1.0/(dt*nfft);
	wave = (struct wave_project *)malloc( nk*sizeof( struct wave_project ));
	for( ik=0 ; ik<nk ; ik++ )
	{
		wave[ik].wavelet     = alloc1float( nfft );
		wave[ik].correlation = alloc1float( nfft );
		wave[ik].k           = 0.;
		wave[ik].intercept   = 0.;
		wave[ik].energy      = 0.;

		zero1float( wave[ik].wavelet     , nfft );
		zero1float( wave[ik].correlation , nfft );
	}
	
	MP_RVF_trace   = alloc1float( nfft );
	MP_RVF_wavelet = alloc1float( nfft );
	MP_RVF_rvf     = alloc1float( nfft );
	MP_RVF_ref     = alloc1float( nfft );
	MP_RVF_wave_rm = alloc1float( nfft );
	MP_RVF_wave_im = alloc1float( nfft );
	innerproduct   = alloc1float( nk );

	zero1float( MP_RVF_trace   , nfft );
	zero1float( MP_RVF_wavelet , nfft );
	zero1float( MP_RVF_rvf     , nfft );
	zero1float( MP_RVF_ref     , nfft );
	zero1float( MP_RVF_wave_rm , nfft );
	zero1float( MP_RVF_wave_im , nfft );

/*******000000wavelet estimation*****************************************************************/
	
	for( it=0 ; it<nfft ; it++ )
		MP_RVF_wavelet[it] = wavelet[it];

/*********00000END*******************************************************************************/
/******1111GENERATE THE ATOMS LIBRARY************************************************************/
	for( ik=0,k_value=k_start ; ik<nk ; ik++, k_value+= dk )
	{
		wave[ik].k         = k_value;
		wave[ik].intercept = intercept;
		zero1float( MP_RVF_rvf     , nfft );
		zero1float( MP_RVF_wave_rm , nfft );
		zero1float( MP_RVF_wave_im , nfft );

		for( iw=0,fre_value=0. ; iw<nfre ; iw++,fre_value+= dfre )
		{
			if( iw==0 )
				MP_RVF_rvf[0]             = 0;
			else if( iw>0 && iw<nfre-1 )
			{
				MP_RVF_rvf[iw]            = intercept + k_value*log(fre_value);
				MP_RVF_rvf[nfft-iw] 	  = MP_RVF_rvf[iw];
			}
			else
				MP_RVF_rvf[iw]            = intercept + k_value*log(fre_value);
		}

		for( it=0 ; it<nfft ; it++ )
		{
			MP_RVF_wave_rm[it] = MP_RVF_wavelet[it];
		}
		lh_fft( MP_RVF_wave_rm , MP_RVF_wave_im , nfft , 1 );
		for( iw=0 ; iw<nfft ; iw++ )
		{
			MP_RVF_wave_rm[iw] = MP_RVF_wave_rm[iw]*MP_RVF_rvf[iw];
			MP_RVF_wave_im[iw] = MP_RVF_wave_im[iw]*MP_RVF_rvf[iw];
		}
		lh_fft( MP_RVF_wave_rm , MP_RVF_wave_im , nfft , -1 );
		for( it=0 ; it<nfft ; it++ )
		{
			wave[ik].wavelet[it] = MP_RVF_wave_rm[it];
		}
		for( it=0 ; it<nfft ; it++ )
		{
			wave[ik].energy+= wave[ik].wavelet[it]*wave[ik].wavelet[it];
		}
		if( wave[ik].energy>0.0000001 )
		{
			for( it=0 ; it<nfft ; it++ )
				wave[ik].wavelet[it] = wave[ik].wavelet[it]/sqrt(wave[ik].energy);
		}
	}
/*******1111END************************************************************************************/
/*******2222**Use the generated atoms library to accomplish Matching Pursuit Method****************/

	for( iy=0 ; iy<ny ; iy++ )
	for( ix=0 ; ix<nx ; ix++ )
	for( isita=0 ; isita<nsita ; isita++ )
	{
printf( "icrossline=%d\tiinline=%d\tisita=%d\n" , iy , ix , isita );

		zero1float( MP_RVF_trace , nfft );
		zero1float( MP_RVF_ref , nfft );

		for( it=0 ; it<nt ; it++ )
		{
			MP_RVF_trace[it] = data[iy][ix][isita][it];
			MP_RVF_ref  [it] = ref [iy][ix][isita][it];
		}

		for( ievent=0 ; ievent<nfft ; ievent++ )
		{
			max_time = lh_max_abs( MP_RVF_ref , nfft );
			if( fabs(MP_RVF_ref[max_time])<0.0000001 )
				break;

			zero1float( innerproduct , nk );

			for( ik=0,k_value=k_start ; ik<nk ; ik++, k_value+= dk )
			{
				lh_shift_array( wave[ik].wavelet , nfft , -max_time );
				innerproduct[ik] = lh_float_iproduct( wave[ik].wavelet , MP_RVF_trace , nfft );
				lh_shift_array( wave[ik].wavelet , nfft , max_time );
			}

			max_ik = lh_max_abs( innerproduct , nk );
			max_k  = k_start+dk*max_ik;
			max_ref_ratio = MP_RVF_ref[max_time]/(intercept+max_k*log(ave_fre));

			lh_shift_array( wave[max_ik].wavelet , nfft , -max_time );
			for( it=0 ; it<nfft ; it++ )
			{
				MP_RVF_trace[it] = MP_RVF_trace[it]-innerproduct[max_ik]*wave[max_ik].wavelet[it];
			}
			lh_shift_array( wave[max_ik].wavelet , nfft , max_time );

			for( iw=1,fre_value=dfre ; iw<nfre ; iw++,fre_value+=dfre )
			{
				if( max_time<nt )
					rvf[iy][ix][isita][iw][max_time]+= (intercept + max_k*log(fre_value))*max_ref_ratio;
			}
			MP_RVF_ref[max_time] = 0.;
		}
		for( it=0 ; it<nt ; it++ )
			error[iy][ix][isita][it] = MP_RVF_trace[it];
	}
/*******2222END************************************************************************************/

	free1float( MP_RVF_trace );
	free1float( MP_RVF_wavelet );
	free1float( MP_RVF_rvf );
	free1float( MP_RVF_wave_rm );
	free1float( MP_RVF_wave_im );
	free( wave );
}

void lh_mp_rvf_outputwave( float ****data	/* 4D Input data[NY][NX][NSITA][NT]*/,
			  int ny		/* Input parameter, NY */,
			  int nx		/* Input parameter, NX */,
			  int nsita		/* Input parameter, NSITA */,
			  int nt		/* Input parameter, NT */,
		  	  float dt		/* Input parameter, DT */,
			  float k_start		/* Scan parameter, start slope */,
			  float k_end		/* Scan parameter, end slope */,
			  int nk		/* Scan parameter, Number of scan slope */,
			  float intercept	/* Scan parameter, value of intercept */,
			  int nevent		/* Scan parameter, Number of events */,
			  int nfft		/* Derived from NT,nfft>=nt  */,
			  float *****rvf	/* Output data, rvf[NY][NX][NSITA][NW][NT] */,
			  float ****error	/* Output data, error[NY][NX][NSITA][NT] */,
			  float *wavelet	/* Input data, wavelet[NFFT] */)
/*< Reflectivity versus Frequency extraction from seismic trace >*/
{
	int ik , iw , it , ievent , ix , iy , isita;
	int max_time=0 , max_ik=0 , nfre=nfft/2+1;
	float *MP_RVF_trace , *MP_RVF_wavelet;
	float k_value , fre_value , max_k=k_start , max_correlation=0. , max_ref_ratio;
	float dk=(k_end-k_start)/nk , dfre;
	float *MP_RVF_rvf , *MP_RVF_wave_rm , *MP_RVF_wave_im;
	struct wave_project *wave;

	dfre = 1.0/(dt*nfft);
	wave = (struct wave_project *)malloc( nk*sizeof( struct wave_project ));
	for( ik=0 ; ik<nk ; ik++ )
	{
		wave[ik].wavelet     = alloc1float( nfft );
		wave[ik].correlation = alloc1float( nfft );
		wave[ik].k           = 0.;
		wave[ik].intercept   = 0.;
		wave[ik].energy      = 0.;

		zero1float( wave[ik].wavelet     , nfft );
		zero1float( wave[ik].correlation , nfft );
	}
	
	MP_RVF_trace   = alloc1float( nfft );
	MP_RVF_wavelet = alloc1float( nfft );
	MP_RVF_rvf     = alloc1float( nfft );
	MP_RVF_wave_rm = alloc1float( nfft );
	MP_RVF_wave_im = alloc1float( nfft );

	zero1float( MP_RVF_trace   , nfft );
	zero1float( MP_RVF_wavelet , nfft );
	zero1float( MP_RVF_rvf     , nfft );
	zero1float( MP_RVF_wave_rm , nfft );
	zero1float( MP_RVF_wave_im , nfft );

/*******000000wavelet estimation*****************************************************************/

	lh_wave_estimation( data ,ny , nx , nsita , nt , MP_RVF_wavelet , nfft );
	for( it=0 ; it<nfft ; it++ )
		wavelet[it] = MP_RVF_wavelet[it];	

/*********00000END*******************************************************************************/
/******1111GENERATE THE ATOMS LIBRARY************************************************************/
	for( ik=0,k_value=k_start ; ik<nk ; ik++, k_value+= dk )
	{
		wave[ik].k         = k_value;
		wave[ik].intercept = intercept;
		zero1float( MP_RVF_rvf     , nfft );
		zero1float( MP_RVF_wave_rm , nfft );
		zero1float( MP_RVF_wave_im , nfft );

		for( iw=0,fre_value=0. ; iw<nfre ; iw++,fre_value+= dfre )
		{
			if( iw==0 )
				MP_RVF_rvf[0]             = 0;
			else if( iw>0 && iw<nfre-1 )
			{
				MP_RVF_rvf[iw]            = intercept + k_value*log(fre_value);
				MP_RVF_rvf[nfft-iw] = MP_RVF_rvf[iw];
			}
			else
				MP_RVF_rvf[iw]            = intercept + k_value*log(fre_value);
		}

		for( it=0 ; it<nfft ; it++ )
		{
			MP_RVF_wave_rm[it] = MP_RVF_wavelet[it];
		}
		lh_fft( MP_RVF_wave_rm , MP_RVF_wave_im , nfft , 1 );
		for( iw=0 ; iw<nfft ; iw++ )
		{
			MP_RVF_wave_rm[iw] = MP_RVF_wave_rm[iw]*MP_RVF_rvf[iw];
			MP_RVF_wave_im[iw] = MP_RVF_wave_im[iw]*MP_RVF_rvf[iw];
		}
		lh_fft( MP_RVF_wave_rm , MP_RVF_wave_im , nfft , -1 );
		for( it=0 ; it<nfft ; it++ )
		{
			wave[ik].wavelet[it] = MP_RVF_wave_rm[it];
		}
		for( it=0 ; it<nfft ; it++ )
		{
			wave[ik].energy+= wave[ik].wavelet[it]*wave[ik].wavelet[it];
		}
		if( wave[ik].energy>0.0000001 )
		{
			for( it=0 ; it<nfft ; it++ )
				wave[ik].wavelet[it] = wave[ik].wavelet[it]/sqrt(wave[ik].energy);
		}
	}
/*******1111END************************************************************************************/
/*******2222**Use the generated atoms library to accomplish Matching Pursuit Method****************/
	for( iy=0 ; iy<ny ; iy++ )
	for( ix=0 ; ix<nx ; ix++ )
	for( isita=0 ; isita<nsita ; isita++ )
	{
		for( it=0 ; it<nt ; it++ )
		{
			MP_RVF_trace[it] = data[iy][ix][isita][it];
		}
		for( ievent=0 ; ievent<nevent ; ievent++ )
		{
			max_time=0;
			max_ik  =0;
			max_k = k_value;
			max_correlation = 0;
			for( ik=0,k_value=k_start ; ik<nk ; ik++, k_value+= dk )
			{
				lh_cross_correlation_ft( MP_RVF_trace , nfft , wave[ik].wavelet , nfft , wave[ik].correlation , nfft );
				for(  it=0 ; it<nfft ; it++ )
				{
					if( fabs(max_correlation)<fabs(wave[ik].correlation[it]) )
					{
						max_correlation = wave[ik].correlation[it];
						max_time        = it;
						max_k           = k_value;
						max_ik          = ik;
						max_ref_ratio   = wave[ik].correlation[it]/sqrt(wave[ik].energy);
					}
				}
			
			}
			lh_shift_array( wave[max_ik].wavelet , nfft , -max_time );
			for( it=0 ; it<nfft ; it++ )
			{
				MP_RVF_trace[it] = MP_RVF_trace[it]-max_correlation*wave[max_ik].wavelet[it];
			}
			lh_shift_array( wave[max_ik].wavelet , nfft , max_time );
			for( iw=1,fre_value=dfre ; iw<nfre ; iw++,fre_value+=dfre )
			{
				if( max_time<nt )
					rvf[iy][ix][isita][iw][max_time]+= (intercept + max_k*log(fre_value))*max_ref_ratio;
			}
		}
		for( it=0 ; it<nt ; it++ )
			error[iy][ix][isita][it] = MP_RVF_trace[it];
	}
/*******2222END************************************************************************************/

	free1float( MP_RVF_trace );
	free1float( MP_RVF_wavelet );
	free1float( MP_RVF_rvf );
	free1float( MP_RVF_wave_rm );
	free1float( MP_RVF_wave_im );
	free( wave );
}
