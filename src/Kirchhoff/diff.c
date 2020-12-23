#include <math.h>

void diff(int n, float *p, float dt)
{
    float y1, y2; int i;

    y1=p[0];
    p[1]=0.0;
    for(i=1; i<n; i++){
        y2=p[i];
        p[i]=(p[i]-y1)/dt;
        y1=y2;
    }    
}    
    
