/* THis is a source program to read&write arrays */
/*
  Copyright (C) 2013 Geophysical Insitute of Sinopec,Nanjing,China
   
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
   
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
   
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

#include <stdio.h>
#include <stdlib.h>

void lh_read_1d_float_bin( float *x , int length , char *fname )
/*< read a 1d array float type >*/
{
	int i;
	FILE *fp1;
	if( !(fp1 = fopen( fname , "rb" )) )
	{
		printf( "Error opening %s for reading. Program terminated.\n" , fname );
		exit(1);
	}
	for( i=0;i<length;i++ )
		fread( &x[i] , sizeof( float ) , 1 , fp1 );
	fclose( fp1 );
}

void lh_read_1d_double_bin( double *x , int length , char *fname )
/*< read a 1d array double type >*/
{
	int i;
	FILE *fp1;
	if( !(fp1 = fopen( fname , "rb" )) )
	{
		printf( "Error opening %s for reading. Program terminated.\n" , fname );
		exit(1);
	}
	for( i=0;i<length;i++ )
		fread( &x[i] , sizeof( double ) , 1 , fp1 );
	fclose( fp1 );
}

void lh_write_1d_float_bin( float *x , int length , char *fname )
/*< write a 1d array float type >*/
{
	int i;
	FILE *fp1;
	
	if( !(fp1 = fopen( fname , "wb" )) )
	{
		printf( "Error opening %s for writing. Program terminated.\n" , fname );
		exit(1);
	}
	for( i=0;i<length;i++ )
		fwrite( &x[i], sizeof( float ) , 1 , fp1 );
	fclose(fp1);
}

void lh_write_1d_double_bin( double *x , int length , char *fname )
/*< write a 1d array double type >*/
{
        int i;
        FILE *fp1;

	if( !(fp1 = fopen( fname , "wb" )) )
	{
		printf( "Error opening %s for writing. Program terminated.\n" , fname );
		exit(1);
	}
        for( i=0;i<length;i++ )
                fwrite( &x[i], sizeof( double ) , 1 , fp1 );
        fclose(fp1);
}

void lh_read_2d_float_bin( float **x , int n1 , int n2 , char *fname )
/*< read a 2d array float type >*/
{
	int i,j;
	FILE *fp1;
	if( !(fp1 = fopen( fname , "rb" )) )
	{
		printf( "Error opening %s for reading. Program terminated.\n" , fname );
		exit(1);
	}
	for( i=0;i<n1;i++ )
	for( j=0;j<n2;j++ )	
		fread( &x[i][j] , sizeof(float) , 1 , fp1 );
	fclose(fp1);
}

void lh_read_2d_float_bin_row( float **x , int n1 , int n2 , char *fname )
/*< read a 2d array float type in row direction for su_ximage >*/
{
	int i,j;
	FILE *fp1;
	if( !(fp1 = fopen( fname , "rb" )) )
	{
		printf( "Error opening %s for reading. Program terminated.\n" , fname );
		exit(1);
	}
	for( j=0;j<n2;j++ )
	for( i=0;i<n1;i++ )
		fread( &x[i][j] , sizeof( float ) , 1 , fp1 );
	fclose(fp1);
}

void lh_write_2d_float_bin( float **x , int n1 , int n2 , char *fname )
/*< write a 2d array float type >*/
{
	int i,j;
	FILE *fp1;
	if( !(fp1 = fopen( fname , "wb" )) )
	{
		printf( "Error opening %s for writing. Program terminated.\n" , fname );
		exit(1);
	}
	for( i=0;i<n1;i++ )
	for( j=0;j<n2;j++ )
		fwrite( &x[i][j] , sizeof(float) , 1 , fp1 );
	fclose(fp1);
}

void lh_write_2d_float_bin_row( float **x , int n1 , int n2 , char *fname )
/*< write a 2d array float type in row direction for su_ximage >*/
{
	int i,j;
	FILE *fp1;
	if( !(fp1 = fopen( fname , "wb" )) )
	{
		printf( "Error opening %s for writing. Program terminated.\n" , fname );
		exit(1);
	}
	for( j=0;j<n2;j++ )
	for( i=0;i<n1;i++ )
		fwrite( &x[i][j] , sizeof( float ) , 1 , fp1 );
	fclose( fp1 );
}

void lh_write_3d_float_bin( float ***a , int n1 , int n2 , int n3 , char *fname )
/*< write a 3d data bin to a file in float type >*/
{
        int i,j,k;

        FILE *fp1;
        fp1 = fopen( fname , "wb" );
        for( i=0 ; i<n1 ; i++ )
        for( j=0 ; j<n2 ; j++ )
        for( k=0 ; k<n3 ; k++ )
                fwrite( &a[i][j][k] , sizeof(float) , 1 , fp1 );

        fclose( fp1 );
}

void lh_read_3d_float_bin( float ***a , int n1 , int n2 , int n3 , char *fname )
/*< read a 3d data bin to a array in float type >*/
{
        int i,j,k;

        FILE *fp1;
        fp1=fopen( fname , "rb" );
        for( i=0 ; i<n1 ; i++ )
        for( j=0 ; j<n2 ; j++ )
        for( k=0 ; k<n3 ; k++ )
                fread( &a[i][j][k] , sizeof(float) , 1 , fp1 );

        fclose( fp1 );
}
