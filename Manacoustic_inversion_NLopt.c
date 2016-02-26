/* Inverse Vp, Rou & Qp Using Simulated Annealing Method************* */

/**********************************************************************
INPUT DATA: seismic data: data[nsita][nt];
	    wavelet[];

OUTPUT DATA: inversed parameters: Vp[nt+1], rou[nt+1] , Qp[nt+1], 

Program developed by: Heng Luo in SINOPEC, 2014/07/25

Gm=d

***********************************************************************/

#include <rsf.h>
#include "lh_mp.h"
#include "su_alloc.h"
#include "lh_readwrite.h"
#include "lh_wavelet.h"
#include "lh_bayesian.h"
#include "lh_LSCG.h"
#include <nlopt.h>

struct wave_project
{
        float *wavelet;
        float k;
        float intercept;
        float energy;
        float *correlation;
};

typedef struct
{
	float **G;
	float *d;
	float *m_prior;
	float *lamda;
	int m;
	int n;
}obj_data;

double myfunc( unsigned n , const double *m , double *grad , void *my_func_data )
{
	obj_data *data = (obj_data*)my_func_data;

	int i;
	float obj_value=0.;
	
	float **G	= data->G;
	float *d	= data->d;
	float *m_prior	= data->m_prior;
	float *lamda	= data->lamda;
	int row		= data->m;
	int column	= data->n;
	
	float *Gm , *Gm_d , *m_m_prior, *m_copy;
	float **Gt, *Gt_Gm_d;

	Gm		= alloc1float( row );
	Gm_d		= alloc1float( row );
	m_m_prior	= alloc1float( column );
	m_copy		= alloc1float( column );

	zero1float( Gm , row );
	zero1float( Gm_d , row );
	zero1float( m_m_prior , column );
	zero1float( m_copy , column );

	if( column!=n ){
		printf( "Parameter of Dimension wrong!\n" );
		return 0;
	}

	for( i=0 ; i<column ; i++ )
		m_copy[i] = (float)m[i];

	lh_matrix_mu_vector( G , m_copy , Gm , row , column );
	lh_vector_sub_vector( Gm , d , Gm_d , row );
	lh_vector_sub_vector( m_copy , m_prior , m_m_prior , column );

	for( i=0 ; i<row ; i++ )
		obj_value+=Gm_d[i]*Gm_d[i];
	for( i=0 ; i<column ; i++ )
		obj_value+=m_m_prior[i]*m_m_prior[i]*lamda[i];

	if(grad)
	{

		Gt		= alloc2float( row , column );
		Gt_Gm_d		= alloc1float( column );

		zero2float( Gt , row , column );
		zero1float( Gt_Gm_d , column );

		lh_matrix_transpose( G , Gt , row , column );
		lh_matrix_mu_vector( Gt , Gm_d , Gt_Gm_d , column , row );

		for( i=0 ; i<column ; i++ )
			grad[i] = (double)2.*( Gt_Gm_d[i]+lamda[i]*m_m_prior[i] );
	}

	free1float( Gm );
	free1float( Gm_d );
	free1float( m_m_prior );
	free1float( m_copy );

	if(grad){
		free2float( Gt );
		free1float( Gt_Gm_d );
	}

	return (double)obj_value;
}

int main( int argc , char* argv[] )
{
	sf_init( argc , argv );

	int nfft=512, nt=98 , i , it , iw , j, nfre=nfft/2+1, nsita=40, isita , sita, dsita=1.*SF_PI/180.;
	float amp=10. , fm=50. , t0=0.08, w0=30., dt=0.002;
	float dfre=1.0/(dt*nfft), fre_value;
	float **trace, **logdata, *vp, *rou, *Qp, *initial_vp, *initial_rou, *initial_qp, *vp_inversed, *rou_inversed, *qp_inversed;
	float **Atoms, **Atoms_sub, **L, **G, *lamda, *d, weight=0.01;
	float *wavelet_E, *wavelet_Q, *MP_RVF_rvf , *MP_RVF_wave_rm , *MP_RVF_wave_im;
	float *m_prior;

	double *lb, *ub;
	double *solution;
	double minf;
	nlopt_opt opt;
	obj_data aux_data;

	opt = nlopt_create( NLOPT_LN_COBYLA , 3*(nt+1) );
	/* SPecify the Algorithm used and the DIMENSION! */
	
	lb		= alloc1double( 3*(nt+1) );
	ub		= alloc1double( 3*(nt+1) );
	solution	= alloc1double( 3*(nt+1) );

	m_prior		= alloc1float( 3*(nt+1) );
	wavelet_E	= alloc1float( nfft );
	wavelet_Q	= alloc1float( nfft );
	MP_RVF_rvf	= alloc1float( nfft );
	MP_RVF_wave_rm	= alloc1float( nfft );
	MP_RVF_wave_im	= alloc1float( nfft );
	Atoms_sub	= alloc2float( 2*nt , nt );
	Atoms		= alloc2float( 2*nt*nsita , nt*nsita );
	L		= alloc2float( 3*(nt+1) , 2*nt*nsita );
	G		= alloc2float( 3*(nt+1) , nt*nsita );
	d		= alloc1float( nsita*nt );
	lamda		= alloc1float( 3*(nt+1) );
	initial_vp	= alloc1float( nt+1 );
	initial_rou	= alloc1float( nt+1 );
	initial_qp	= alloc1float( nt+1 );
	vp_inversed	= alloc1float( nt+1 );
	rou_inversed	= alloc1float( nt+1 );
	qp_inversed	= alloc1float( nt+1 );

	zero1double( lb , 3*(nt+1) );
	zero1double( ub , 3*(nt+1) );
	zero1double( solution , 3*(nt+1) );

	zero1float( m_prior , 3*(nt+1) );
	zero1float( wavelet_E		, nfft );
	zero1float( wavelet_Q		, nfft );
	zero1float( MP_RVF_rvf		, nfft );
	zero1float( MP_RVF_wave_rm	, nfft );
	zero1float( MP_RVF_wave_im	, nfft );
	zero2float( Atoms_sub , 2*nt , nt );
	zero2float( Atoms , 2*nt*nsita , nt*nsita );
	zero2float( L , 3*(nt+1) , 2*nt*nsita );
	zero2float( G , 3*(nt+1) , nt*nsita );
	zero1float( d , nsita*nt );
	zero1float( lamda , 3*(nt+1) );
	zero1float( initial_vp , nt+1 );
	zero1float( initial_rou , nt+1 );
	zero1float( initial_qp , nt+1 );
	zero1float( vp_inversed , nt+1 );
	zero1float( rou_inversed , nt+1 );
	zero1float( qp_inversed , nt+1 );

	logdata	= alloc2float( nt+1 , 6 );
	trace	= alloc2float( nt , 40 );
	vp	= alloc1float( nt+1 );
	rou	= alloc1float( nt+1 );
	Qp	= alloc1float( nt+1 );

	zero2float( logdata , nt+1 , 6 );
	zero2float( trace , nt , 40 );
	zero1float( vp , nt+1 );
	zero1float( rou , nt+1 );
	zero1float( Qp , nt+1 );

	lh_read_2d_float_bin( logdata , 6 , nt+1 , "/home/Pitaloveu/DeskTop/log_viscoelastic_resampling.bin");
	lh_read_2d_float_bin( trace , 40 , nt , "/var/tmp/anacoustic_trace.rsf@" );

	for( i=0 ; i<nt+1 ; i++ )
	{
		vp[i] = log(logdata[1][i]);
		rou[i]= log(logdata[3][i]);
		Qp[i] = logdata[4][i];
	}

	lh_write_1d_float_bin( vp , (nt+1) , "vp.bin" );
	lh_write_1d_float_bin( rou, (nt+1) , "rou.bin" );
	lh_write_1d_float_bin( Qp , (nt+1) , "Qp.bin" );

	for( isita=0 ; isita<nsita ; isita++ )
	for( it=0    ; it<nt       ; it++    )
		d[isita*nt+it] = trace[isita][it];

	lh_read_1d_float_bin( initial_qp , nt+1 , "qp_smooth.bin" );
	lh_read_1d_float_bin( initial_vp , nt+1 , "vp_smooth.bin" );
	lh_read_1d_float_bin( initial_rou , nt+1 , "rou_smooth.bin" );

	for( it=0 ; it<nt+1 ; it++ )
	{
//		m_prior[it]		= initial_vp[it];
//		m_prior[nt+1+it]	= initial_rou[it];
		m_prior[2*nt+2+it]	= initial_qp[it];

		m_prior[it]		= vp[it];
		m_prior[nt+1+it]	= rou[it];
//		m_prior[2*nt+2+it]	= Qp[it];
	}

	lh_ricker( wavelet_E , amp , nfft , dt , fm , t0 );

	for( iw=0,fre_value=0. ; iw<nfre ; iw++,fre_value+= dfre )
	{
		if( iw==0 )
			MP_RVF_rvf[0]             = 0;
		else if( iw>0 && iw<nfre-1 )
		{
			MP_RVF_rvf[iw]            = log(fre_value/w0);
			MP_RVF_rvf[nfft-iw]       = MP_RVF_rvf[iw];
		}
		else
			MP_RVF_rvf[iw]            = log(fre_value/w0);
	}
	for( it=0 ; it<nfft ; it++ )
	{
		MP_RVF_wave_rm[it] = wavelet_E[it];
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
		wavelet_Q[it] = MP_RVF_wave_rm[it];
	}

	for( isita=0 ,sita=0 ; isita<nsita ; isita++ ,sita+=dsita )
	for( it=0            ; it<nt       ; it++                 )
	{
		L[isita*2*nt+it][it]			= -0.5/( cos(sita)*cos(sita) );
		L[isita*2*nt+it][it+1]			=  0.5/( cos(sita)*cos(sita) );

		L[isita*2*nt+it][nt+1+it]		= -0.5;
		L[isita*2*nt+it][nt+1+it+1]		=  0.5;

		L[isita*2*nt+nt+it][2*nt+2+it]		= -0.5/( cos(sita)*cos(sita)*SF_PI );
		L[isita*2*nt+nt+it][2*nt+2+it+1]	=  0.5/( cos(sita)*cos(sita)*SF_PI );
	}

	lh_shift_array( wavelet_E , nfft , (int)(t0/dt+0.5) );
	lh_shift_array( wavelet_Q , nfft , (int)(t0/dt+0.5) );

	for( i=0 ; i<nt ; i++ )
	{
		for( j=0 ; j<nt ; j++ )
		{
			Atoms_sub[j][i]		= wavelet_E[j];
			Atoms_sub[j][nt+i]	= wavelet_Q[j];
		}
		lh_shift_array( wavelet_E , nfft , -1 );
		lh_shift_array( wavelet_Q , nfft , -1 );
	}

	for( isita=0 ; isita<nsita ; isita++ )
	for( i=0     ; i<nt        ; i++     )
	for( j=0     ; j<2*nt      ; j++     )
		Atoms[isita*nt+i][isita*2*nt+j] = Atoms_sub[i][j];

	lh_matrix_mu_matrix( Atoms , L , G , nt*nsita , 2*nt*nsita , 3*(nt+1) );

	aux_data.G		= G;
	aux_data.d		= d;
	aux_data.m_prior	= m_prior;
	aux_data.lamda		= lamda;
	aux_data.m		= nt*nsita;
	aux_data.n 		= 3*(nt+1);

	for( i=0 ; i<nt+1 ; i++ )
	{
		lb[i]		= 8.0;
		lb[nt+1+i]	= 7.7;
		lb[2*nt+2+i]	= 0.;

		ub[i]		= 8.6;
		ub[nt+1+i]	= 7.9;
		ub[2*nt+2+i]	= 0.14;
	}

	for( i=0 ; i<3*(nt+1) ; i++ )
	{
//		lamda[i] = weight/pow( (lb[i]+ub[i])/2., 2. );
		solution[i] = ( lb[i]+ub[i] )/2.;
	}
for( i=0 ; i<nt+1 ; i++ )
{
lamda[i] = weight;
lamda[i+nt+1] = weight;
lamda[i+2*nt+2] = weight;
}

	nlopt_set_lower_bounds( opt , lb );
	nlopt_set_upper_bounds( opt , ub );

	nlopt_set_min_objective( opt , myfunc , &aux_data );

	nlopt_set_xtol_rel( opt , 1e-8 );

	if( nlopt_optimize( opt , solution , &minf )<0 )
		printf( "nlopt failed\n" );

	for( i=0 ; i<nt+1 ; i++ )
	{
		vp_inversed[i]	= (float)solution[i];
		rou_inversed[i]	= (float)solution[nt+1+i];
		qp_inversed[i]	= (float)solution[2*nt+2+i];
	}

	lh_write_1d_float_bin( vp_inversed	, nt+1 , "vp_inversed.bin" );
	lh_write_1d_float_bin( rou_inversed	, nt+1 , "rou_inversed.bin" );
	lh_write_1d_float_bin( qp_inversed	, nt+1 , "qp_inversed.bin" );

	nlopt_destroy( opt );

	exit(0);
}
