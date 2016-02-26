#include <rsf.h>
#include "tl_bayesian.h"

void tl_normalize(float *a	/*Arrary A*/,
		  float *b 	/*Array B*/,
		  int n		/*length of A*/)
/*< normalize an array using max value >*/
{
	int i; 
	float  max;
	max=fabs(a[0]);
	for (i=0; i<n; i++)
		if (fabs(a[i])>max) max=fabs(a[i]);
	for (i=0; i<n; i++)
		b[i]=a[i]/max;	
}

void tl_normalize2d(float **a   /*Array A*/,
		    float **b   /*Array B*/,
		    int n1 	/*Fast dimention*/,
		    int n2 	/*Slow dimention*/)
/*< normalize a 2D array using max value >*/
{
	int i, j;
	float max;
	max=fabs(a[0][0]);
	for (i=0; i<n2; i++)
	for (j=0; j<n1; j++)
	if (fabs(a[i][j])>max) max=fabs(a[i][j]);

	for (i=0; i<n2; i++)
        for (j=0; j<n1; j++)
	b[i][j]=a[i][j]/max;	
}

