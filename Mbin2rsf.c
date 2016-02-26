/* bin2rsf transfers the binary data into .rsf style,
the binary data is native_float type and the DIMENSION parameters
are given in a line command */

/* COPYRIGHT RESERVED BY SINOPEC @COP */

#include <rsf.h>
#include "lh_readwrite.h"


int main( int argc , char* argv[] )
{
	sf_init( argc , argv );

	int dim=1;
	int n1,n2=1,n3=1,n4=1;
	float d1=0.,d2=0.,d3=0.,d4=0.;
	float o1=0.,o2=0.,o3=0.,o4=0.;
	char *label1=NULL, *label2=NULL, *label3=NULL, *label4=NULL, *unit1=NULL, *unit2=NULL, *unit3=NULL, *unit4=NULL;
	float *data;
	char *input=NULL;
	sf_file Fr=NULL;

	Fr = sf_output( "out" );

	if(!sf_getint( "n1" , &n1 ))	sf_error( "The first dimension of the binary data should be fixed!\n" );
	/* n1 is the Dimension of the fastest axis */
	if(sf_getint( "n2" , &n2 ))	dim++;
	/* n2 is the Dimension of the second axis */
	if(sf_getint( "n3" , &n3 ))	dim++;
	/* n3 is the Dimension of the third axis */
	if(sf_getint( "n4" , &n4 ))	dim++;
	/* n4 is the Dimension of the forth axis */

	if(!sf_getfloat( "d1" , &d1 ))	;
	/* d1 is the interval of the fastest dimension of the binary data */
	if(!sf_getfloat( "d2" , &d2 ))	;
	/* d2 is the interval of the second dimension of the bianry data */
	if(!sf_getfloat( "d3" , &d3 ))	;
	/* d3 is the interval of the third dimension of the bianry data */
	if(!sf_getfloat( "d4" , &d4 ))	;
	/* d4 is the interval of the forth dimension of the bianry data */

	if(!sf_getfloat( "o1" , &o1 ))	;
	/* o1 is the starting value of the fastest dimension of the bianry data */
	if(!sf_getfloat( "o2" , &o2 ))	;
	/* o2 is the starting value of the second dimension of the bianry data */
	if(!sf_getfloat( "o3" , &o3 ))	;
	/* o3 is the starting value of the third dimension of the bianry data */
	if(!sf_getfloat( "o4" , &o4 ))	;
	/* o4 is the starting value of the forth dimension of the bianry data */

	label1 = sf_getstring( "label1" );
	label2 = sf_getstring( "label2" );
	label3 = sf_getstring( "label3" );
	label4 = sf_getstring( "label4" );

	unit1  = sf_getstring( "unit1" );
	unit2  = sf_getstring( "unit2" );
	unit3  = sf_getstring( "unit3" );
	unit4  = sf_getstring( "unit4" );

	input=sf_getstring( "input" );

	data = sf_floatalloc( n1*n2*n3*n4 );
	lh_read_1d_float_bin( data , n1*n2*n3*n4 , input );	

	if( dim==1 )
	{
		sf_putint	( Fr , "n1" , n1 );
		sf_putfloat	( Fr , "d1" , d1 );
		sf_putfloat	( Fr , "o1" , o1 );
	}
	if( dim==2 )
	{
		sf_putint	( Fr , "n1" , n1 );
		sf_putfloat	( Fr , "d1" , d1 );
		sf_putfloat	( Fr , "o1" , o1 );
		
		sf_putint	( Fr , "n2" , n2 );
		sf_putfloat	( Fr , "d2" , d2 );
		sf_putfloat	( Fr , "o2" , o2 );
	}
	if( dim==3 )
	{
		sf_putint	( Fr , "n1" , n1 );
		sf_putfloat	( Fr , "d1" , d1 );
		sf_putfloat	( Fr , "o1" , o1 );
		
		sf_putint	( Fr , "n2" , n2 );
		sf_putfloat	( Fr , "d2" , d2 );
		sf_putfloat	( Fr , "o2" , o2 );

		sf_putint	( Fr , "n3" , n3 );
		sf_putfloat	( Fr , "d3" , d3 );
		sf_putfloat	( Fr , "o3" , o3 );
	}
	if( dim==4 )
	{
		sf_putint	( Fr , "n1" , n1 );
		sf_putfloat	( Fr , "d1" , d1 );
		sf_putfloat	( Fr , "o1" , o1 );
		
		sf_putint	( Fr , "n2" , n2 );
		sf_putfloat	( Fr , "d2" , d2 );
		sf_putfloat	( Fr , "o2" , o2 );

		sf_putint	( Fr , "n3" , n3 );
		sf_putfloat	( Fr , "d3" , d3 );
		sf_putfloat	( Fr , "o3" , o3 );

		sf_putint	( Fr , "n4" , n1 );
		sf_putfloat	( Fr , "d4" , d1 );
		sf_putfloat	( Fr , "o4" , o1 );
	}
	sf_floatwrite( data , n1*n2*n3*n4 , Fr );

}
