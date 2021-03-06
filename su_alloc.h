/* This file is automatically generated. DO NOT EDIT! */

#ifndef _su_alloc_h
#define _su_alloc_h


int *alloc1int(size_t n1);
/*< allocate a 1-d array of ints >*/


void free1int(int *p);
/*< free a 1-d array of ints >*/


int **alloc2int(size_t n1, size_t n2);
/*< allocate a 2-d array of ints >*/


void free2int(int **p);
/*< free a 2-d array of ints >*/


int ***alloc3int(size_t n1, size_t n2, size_t n3);
/*< allocate a 3-d array of ints >*/


void free3int(int ***p);
/*< free a 3-d array of ints >*/


float *alloc1float(size_t n1);
/*< allocate a 1-d array of floats >*/


void free1float(float *p);
/*< free a 1-d array of floats >*/


float **alloc2float(size_t n1, size_t n2);
/*< allocate a 2-d array of floats >*/


void free2float(float **p);
/*< free a 2-d array of floats >*/


float ***alloc3float(size_t n1, size_t n2, size_t n3);
/*< allocate a 3-d array of floats >*/


void free3float(float ***p);
/*< free a 3-d array of floats >*/


double *alloc1double(size_t n1);
/*< allocate a 1-d array of doubles >*/


void free1double(double *p);
/*< free a 1-d array of doubles >*/


double **alloc2double(size_t n1, size_t n2);
/*< allocate a 2-d array of doubles >*/


void free2double(double **p);
/*< free a 2-d array of doubles >*/


double ***alloc3double(size_t n1, size_t n2, size_t n3);
/*< allocate a 3-d array of doubles >*/


void free3double(double ***p);
/*< free a 3-d array of doubles >*/


void zero1int(int *p, size_t n1);
/*< set the arrays to zero >*/


void zero2int(int **p, size_t n1, size_t n2);
/*< set the arrays to zero >*/


void zero3int(int ***p, size_t n1, size_t n2, size_t n3);
/*< set the arrays to zero >*/


void zero1float(float *p, size_t n1);
/*< set the arrays to zero >*/


void zero2float(float **p, size_t n1, size_t n2);
/*< set the arrays to zero >*/


void zero3float(float ***p, size_t n1, size_t n2, size_t n3);
/*< set the arrays to zero >*/


void zero1double(double *p, size_t n1);
/*< set the arrays to zero >*/


void zero2double(double **p, size_t n1, size_t n2);
/*< set the arrays to zero >*/


void zero3double(double ***p, size_t n1, size_t n2, size_t n3);
/*< set the arrays to zero >*/


void zero4float(float ****p, size_t n1, size_t n2, size_t n3, size_t n4);
/*< Zeror the 4d float array >*/


void zero5float(float *****p, size_t n1, size_t n2, size_t n3, size_t n4 , size_t n5 );
/*< Zeror the 5d float array >*/

#endif
