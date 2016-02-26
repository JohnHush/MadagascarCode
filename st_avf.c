#include <rsf.h>
#include "st_avf.h"
#include "su_alloc.h"
#include "lh_wavelet.h"

void winft( float *data_in      /* Input wavelet */,
            float *data_out     /* Output wavelet */,
            int nt              /* length of wavelet */,
            int nfft            /* Length of FFT */,
            float dt            /* Time Interval of TIME DOMAIN */,
            float fre           /* Window Parameter, central frequency */,
            float var           /* Window Parameter, Variance of the Window */)
/*< Windowed FFT & IFFT >*/
{
    float *rm, *im, *ft_win;
    int i;
    float dw=1./(nfft*dt);

    rm = alloc1float( nfft );
    im = alloc1float( nfft );

    ft_win = alloc1float( 1+nfft/2);

    zero1float( rm , nfft );
    zero1float( im , nfft );

    zero1float( ft_win , 1+nfft/2 );

    for( i=0 ; i<nt ; i++ )
        rm[i] = data_in[i];

    fft( rm , im , nfft , 1 );

    for( i=0 ; i<1+nfft/2 ; i++ )
    {
        ft_win[i] = (1./(var*sqrt(2*3.14159265359)))*exp(-1.*(dw*i-fre)*(dw*i-fre)/(2.*var*var));
        rm[i]     = rm[i]*ft_win[i];
        im[i]     = im[i]*ft_win[i];
    }
    for( i=1+nfft/2 ; i<nfft ; i++ )
    {
        rm[i]     = rm[nfft-i];
        im[i]     = -im[nfft-i];
    }
    fft( rm , im , nfft , -1 );

    for( i=0 ; i<nt ; i++ )
        data_out[i] = rm[i];

    free1float( rm );
    free1float( im );
    free1float( ft_win );
}

void wavelet_filtered(	float *wavelet		/* The Original Wavelet */,
			float *wave_filtered	/* The Filtered Wavelet */,
			int n			/* Length Of Wavelet */,
			float dt		/* Time Interval of the Wavelets */,
			float fre		/* The Window Parameter, central frequency*/,
			float variance		/* The Window Parameter, Variance of the window */)
/*< Compute the filtered wavelet using Fourier transform >*/
{
	int i;

	float *rm , *im;
	float *tf_win;
	float *win;
	int frequency=(int)(fre*(n*dt));

	rm = alloc1float( n );
	im = alloc1float( n );
	
	tf_win = alloc1float( n/2+1 );
	win    = alloc1float( n );

	zero1float( rm , n );
	zero1float( im , n );

	zero1float( tf_win , n/2+1 );
	zero1float( win , n );

	if( frequency<0 || frequency >n/2 )
	{
		printf( "Wrong Frequency Para!\n" );
		exit( 0 );
	}

	for( i=0 ; i<n ; i++ )
		rm[i] = wavelet[i];

	fft( rm , im , n , 1 );

	for( i=0 ; i<n/2+1 ; i++ )
	{
		tf_win[i] = (1./(variance*sqrt(2*3.14159265359)))*exp(-1.*(i-frequency)*(i-frequency)/(2.*variance*variance));
	}
	for( i=0 ; i<n ; i++ )
	{
		if( i<n/2+1 )
			win[i] = tf_win[i];
		else
			win[i] = tf_win[n-i];
	}
	for( i=0 ; i<n ; i++ )
	{
		rm[i] = rm[i]*win[i];
		im[i] = im[i]*win[i];
	}
	fft( rm , im , n , -1 );
	for( i=0 ; i<n ; i++ )
		wave_filtered[i] = rm[i];

	free1float( rm );
	free1float( im );
	free1float( tf_win );
	free1float( win );
}

void cross_correlation_ft( float *x , int nx , float *y , int ny , float *cross , int nc )
/*< Compute the correlation between array x and array y >*/
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

	fft( rm_x , im_x , fft_value , 1 );
	fft( rm_y , im_y , fft_value , 1 );
	for( i=0 ; i<fft_value ; i++ )
		im_y[i] = -im_y[i];
	for( i=0 ; i<fft_value ; i++ )
	{
		rm_cross[i] = creal( (rm_x[i] + im_x[i]*I)*(rm_y[i] + im_y[i]*I) );
		im_cross[i] = cimag( (rm_x[i] + im_x[i]*I)*(rm_y[i] + im_y[i]*I) );
	}

	fft( rm_cross , im_cross , fft_value , -1 );
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

void move_array( float *x , int n , int step )
/*< Move array x in step >*/
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

void move_spec( float *x , float *y , int n , int step )
/*< Move array x & y in step >*/
{
        int i;
        float *tmp1; tmp1=alloc1float( n ); zero1float( tmp1 , n );
        float *tmp2; tmp2=alloc1float( n ); zero1float( tmp2 , n );

        if( step>0 )
        {
                for( i=0;i<n;i++ )
                {
                        if( i<n-step )
                        {
                                tmp1[i]=x[i+step]; tmp2[i]=y[i+step];
                        }
                        else
                        {
                                tmp1[i]=x[i-n+step]; tmp2[i]=y[i-n+step];
                        }
                }
                for( i=0;i<n;i++ )
                {
                        x[i]=tmp1[i]; y[i]=tmp2[i];
                }
        }
        else if( step==0 )
        {
                for( i=0;i<n;i++ )
                {
                        x[i]=x[i]; y[i]=y[i];
                }
        }
	else
        {
                step=-step;
                for( i=0;i<n;i++ )
                {
                        if( i<step )
                        {
                                tmp1[i]=x[n-step+i]; tmp2[i]=y[n-step+i];
                        }
                        else
                        {
                                tmp1[i]=x[i-step]; tmp2[i]=y[i-step];
                        }
                }
                for( i=0;i<n;i++ )
                {
                        x[i]=tmp1[i]; y[i]=tmp2[i];
                }
        }
        free1float( tmp1 );
        free1float( tmp2 );
}

void gst_avf( float *gt , float **s_real , float **s_imag , float tao , int n , float min_fre , int nw )
/*< Compute time-frequency using S-transform >*/
{
	int i , j, iw , iw_start , iw_end;
	float dw=1.0/(tao*n);

	iw_start = (int)(min_fre/dw);
	iw_end   = iw_start+nw-1;
//	printf("n=%d\n", n);
	if( iw_end >n/2 )
	{
		printf( "iw_end is too big!\n" );
		//return 1;
	}
	zero2float( s_real , n , n );
	zero2float( s_imag , n , n );

	for( i=0 ; i<n ; i++ )
	{
		if( i<=n/2 )
		{
			iw = n/2-i;
			if( iw<=iw_end && iw>=iw_start )
				gst_fre( gt , &s_real[i][0] , &s_imag[i][0] , tao , n , iw );
		}
		else
		{
			for( j=0 ; j<n ; j++ )
			{
				s_real[i][j] = s_real[n-i][j];
				s_imag[i][j] = -s_imag[n-i][j];
			}
		}
	}
	for( i=0 ; i<n ; i++ )
	for( j=0 ; j<n ; j++ )
		s_imag[i][j] = -s_imag[i][j];
	//return 0;
	
}

void gst_fre( float *gt , float *s_real , float *s_imag , float tao , int n , int nw )
/*< Compute time-frequency using S-Transform subordinate programme >*/
{
	int i;
	float *gt_real , *gt_imag , dw=1.0/(n*tao) , cru_fre , *tf_win , *pr_tfwin , *pi_tfwin;

	cru_fre = dw*nw;

	gt_real    = alloc1float( n );
	gt_imag    = alloc1float( n );
	tf_win     = alloc1float( n );
	pr_tfwin   = alloc1float( n );
	pi_tfwin   = alloc1float( n );

	zero1float( gt_real  , n );
	zero1float( gt_imag  , n );
	zero1float( tf_win   , n );
	zero1float( pr_tfwin , n );
	zero1float( pi_tfwin , n );

	for( i=0 ; i<n ; i++ )
	{
		gt_real[i] = gt[i];
	}
	fft( gt_real , gt_imag , n , 1 );
	move_spec( gt_real , gt_imag , n , nw );

	for( i=0 ; i<n ; i++ )
	{
		if( i<=n/2 )
			tf_win[i]   = (fabs(cru_fre)/sqrt(2.0*3.14159265359))*exp(-0.5*cru_fre*cru_fre*tao*tao*i*i);
		else
			tf_win[i]   = tf_win[n-i];
		pr_tfwin[i] = tf_win[i];
	}
	fft( pr_tfwin , pi_tfwin , n , 1 );

	for( i=0 ; i<n ; i++ )
	{
		s_real[i]=crealf( (pr_tfwin[i]+pi_tfwin[i]*I) * (gt_real[i]+gt_imag[i]*I) );
		s_imag[i]=cimagf( (pr_tfwin[i]+pi_tfwin[i]*I) * (gt_real[i]+gt_imag[i]*I) );
	}
	fft( s_real , s_imag , n , -1 );
	
	free1float( gt_real );
	free1float( gt_imag );
	free1float( tf_win );
	free1float( pr_tfwin );
	free1float( pi_tfwin );
}

void gst_fre_plus( float *gt , float *s_real , float *s_imag , float tao , int k , int nw )
/*< Compute time-frequency using S-Transform subordinate programme >*/
{
	int i,n;
    for( n=1,i=0;i<16;i++)
    {
        n=2*n;
        if(n>=k)
            break;
    }
	float *gt_real , *gt_imag , dw=1.0/(n*tao) , cru_fre , *tf_win , *pr_tfwin , *pi_tfwin, *s_real_aux,*s_imag_aux;

	cru_fre = dw*nw;

	gt_real    = alloc1float( n );
	gt_imag    = alloc1float( n );
    s_real_aux = alloc1float( n );
    s_imag_aux = alloc1float( n );
	tf_win     = alloc1float( n );
	pr_tfwin   = alloc1float( n );
	pi_tfwin   = alloc1float( n );

	zero1float( gt_real  , n );
	zero1float( gt_imag  , n );
	zero1float( s_real_aux , n );
	zero1float( s_imag_aux , n );
	zero1float( tf_win   , n );
	zero1float( pr_tfwin , n );
	zero1float( pi_tfwin , n );

	for( i=0 ; i<k ; i++ )
	{
		gt_real[i] = gt[i];
	}
	fft( gt_real , gt_imag , n , 1 );
	move_spec( gt_real , gt_imag , n , nw );

	for( i=0 ; i<n ; i++ )
	{
		if( i<=n/2 )
			tf_win[i]   = (fabs(cru_fre)/sqrt(2.0*3.14159265359))*exp(-0.5*cru_fre*cru_fre*tao*tao*i*i);
		else
			tf_win[i]   = tf_win[n-i];
		pr_tfwin[i] = tf_win[i];
	}
	fft( pr_tfwin , pi_tfwin , n , 1 );

	for( i=0 ; i<n ; i++ )
	{
		s_real_aux[i]=crealf( (pr_tfwin[i]+pi_tfwin[i]*I) * (gt_real[i]+gt_imag[i]*I) );
		s_imag_aux[i]=cimagf( (pr_tfwin[i]+pi_tfwin[i]*I) * (gt_real[i]+gt_imag[i]*I) );
	}
	fft( s_real_aux , s_imag_aux , n , -1 );
    for( i=0 ; i<k ; i++ )
    {
        s_real[i] = s_real_aux[i];
        s_imag[i] = s_imag_aux[i];
    }
	
	free1float( gt_real );
	free1float( gt_imag );
	free1float( tf_win );
	free1float( pr_tfwin );
	free1float( pi_tfwin );
    free1float( s_real_aux);
    free1float( s_imag_aux);
}

void fft( float x[], float y[], int n,int sign)
/*< Fourier Transform >*/	
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

void gist_filtered(	float *gt	/* The filtered Trace */,
			float **s_real	/* The Input real part of ST Spectrum */,
			float **s_imag	/* The Input Imag part of ST Spectrum */,
			int n		/* Length of Spectrum */,
			float dt	/* Time Interval of the data */,
			float fre	/* The Window Parameter, Central frequency */,
			float variance	/* The Window Parameter, Variance of the Window */)
/*< Compute frequency-depended seismic data using inverse S-Transform >*/
{
        int i,j;
        float *tf_win;
        float *gt_real,*gt_imag;
        float *win;
	int frequency=(int)(fre*(n*dt));

        tf_win = alloc1float( n/2+1 );
        win    = alloc1float( n );
        gt_real= alloc1float( n );
        gt_imag= alloc1float( n );

        zero1float( tf_win , n/2+1 );
        zero1float( gt_real , n );
        zero1float( gt_imag , n );
        zero1float( win , n );

        if( frequency<0 || frequency>n/2 )
        {
                printf( "Wrong frequency Para!\n" );
                exit(0);
        }

        for( i=0 ; i<n/2+1 ; i++ )
        {
                tf_win[i] = (1./(variance*sqrt(2*3.14159265359)))*exp(-1.*(i-frequency)*(i-frequency)/(2.*variance*variance));
        }
        for( i=0 ; i<n ; i++ )
        {
                if( i<n/2-1 )
                        win[i] = tf_win[n-i-n/2-1];
                else
                        win[i] = tf_win[i-n/2+1];
        }
        for( i=0 ; i<n ; i++ )
        for( j=0 ; j<n ; j++ )
        {
		gt_real[i]+=s_real[i][j]*win[i];
		gt_imag[i]+=s_imag[i][j]*win[i];
//		gt_real[i]+=s_real[i][j];
//	        gt_imag[i]+=s_imag[i][j];
        }

	lh_shift_array( gt_real , n , n/2 );
	lh_shift_array( gt_imag , n , n/2 );

        fft( gt_real , gt_imag , n , -1 );
        for( i=0 ; i<n ; i++ )
                gt[i] = gt_real[i];

        free1float( tf_win );
        free1float( gt_real );
        free1float( gt_imag );
        free1float( win );

}

