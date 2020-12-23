
#include <math.h>

float minfloat(int n, float *x)
{
    float X = 1.e30;

    int i;

    for(i=0; i < n; i++)
	if(X > fabs(x[i])) X = fabs(x[i]);
    
    return X;
}

float maxfloat(int n, float *x)
{
    float X = -1.e30;

    int i;

    for(i=0; i < n; i++)
	if(X < fabs(x[i])) X = fabs(x[i]);
    
    return X;
}

