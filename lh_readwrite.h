/* This file is automatically generated. DO NOT EDIT! */

#ifndef _lh_readwrite_h
#define _lh_readwrite_h


void lh_read_1d_float_bin( float *x , int length , char *fname );
/*< read a 1d array float type >*/


void lh_read_1d_double_bin( double *x , int length , char *fname );
/*< read a 1d array double type >*/


void lh_write_1d_float_bin( float *x , int length , char *fname );
/*< write a 1d array float type >*/


void lh_write_1d_double_bin( double *x , int length , char *fname );
/*< write a 1d array double type >*/


void lh_read_2d_float_bin( float **x , int n1 , int n2 , char *fname );
/*< read a 2d array float type >*/


void lh_read_2d_float_bin_row( float **x , int n1 , int n2 , char *fname );
/*< read a 2d array float type in row direction for su_ximage >*/


void lh_write_2d_float_bin( float **x , int n1 , int n2 , char *fname );
/*< write a 2d array float type >*/


void lh_write_2d_float_bin_row( float **x , int n1 , int n2 , char *fname );
/*< write a 2d array float type in row direction for su_ximage >*/


void lh_write_3d_float_bin( float ***a , int n1 , int n2 , int n3 , char *fname );
/*< write a 3d data bin to a file in float type >*/


void lh_read_3d_float_bin( float ***a , int n1 , int n2 , int n3 , char *fname );
/*< read a 3d data bin to a array in float type >*/

#endif