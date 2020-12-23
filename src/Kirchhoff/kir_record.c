#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct vel_stress
{
    float u;
    float w;
    float t11;
    float t12;
    float t22;
};

/* kirrecord outputs two components of the wavefield from 
   the finite-difference. Note that the curl curl of Us is 
   actually the dirivative
*/   

kirrecord(p,ldmu,iz,h,dt,nx,nl,nh,kdx,kirout,source)
struct vel_stress *p;
float h, dt;
int iz,nx,nl, nh, kdx,source;
FILE *kirout;
{
    register int i;
    register struct vel_stress *pp;
    float *hold;
    int knx, offset;


    knx = (nx-nl-nh-3)/kdx;
    hold= (float *) malloc(sizeof(float)*knx);

/*
    ldmu=(ta(0)+tb(0))/2.0*h/dt;
    ldmu=(2.6*(5*5-3*3))*2.0;
*/
    if(source==0 || source ==2){    		
	for(i=0; i< knx; i++)		/* (gra div u)x */
	{
	    offset   = nl+1 + i*kdx +iz*nx;
	    pp       = p      + offset;
	    
	    hold[i]  = (pp[  1].t11 - pp[  -1].t11)* (9.0/8.0)
	              -(pp[  2].t11 - pp[  -2].t11)* (1.0/24.0); 
	    hold[i] += (pp[  1].t22 - pp[  -1].t22)* (9.0/8.0)
	              -(pp[  2].t22 - pp[  -2].t22)* (1.0/24.0); 
	    hold[i] /= (ldmu*h);          
	    if(fabs(hold[i])>1.)fprintf(stdout,"Error1 i=%d\n",i);
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* z-derivative */
	{
	    offset   = nl+1 + i*kdx +iz*nx;
	    pp       = p      + offset;
	    
	    hold[i]  = (pp[ nx+1].t11 - pp[ nx-1].t11)
	              -(pp[-nx+1].t11 - pp[-nx-1].t11);
	    hold[i] += (pp[ nx+1].t22 - pp[ nx-1].t22)
	              -(pp[-nx+1].t22 - pp[-nx-1].t22);
	    hold[i] /=  (4*ldmu*h*h);
	    if(fabs(hold[i])>1.)fprintf(stdout,"Error2 i=%d\n",i);
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* (gra div u)z */
	{
	    offset   = nl+1 + i*kdx +iz*nx;
	    pp       = p      + offset;
	    
/*
	    hold[i]  = (pp[  nx].t11 - pp[  -nx].t11)* (9.0/8.0)
	              -(pp[2*nx].t11 - pp[-2*nx].t11)* (1.0/24.0); 
	    hold[i]  = (pp[  nx].t22 - pp[  -nx].t22)* (9.0/8.0)
	              -(pp[2*nx].t22 - pp[-2*nx].t22)* (1.0/24.0); 
	    hold[i] *=  0.5/(ta(0)+tb(0))*dt/h/h;          
*/
	    hold[i]  =(pp[nx].t11-pp[-nx].t11)+(pp[nx].t22-pp[-nx].t22);
	    hold[i] /= (2.*ldmu*h);
	    if(fabs(hold[i])>1.)fprintf(stdout,"Error3 i=%d\n",i);
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* z-derivative */
	{
	    offset   = nl+1 + i*kdx +iz*nx;
	    pp       = p      + offset;
	    
	    hold[i]  = (pp[nx].t11 - 2.*pp[0].t11 + pp[-nx].t11)
	              +(pp[nx].t22 - 2.*pp[0].t22 + pp[-nx].t22); 
	    hold[i] /=  (ldmu*h*h);          
	    if(fabs(hold[i])>1.)fprintf(stdout,"Error4 i=%d\n",i);
	}
	fwrite(hold,sizeof(float),knx,kirout);
    }
    if(source==1 || source ==2){	
	for(i=0; i< knx; i++)		/* (curl curl u)x */
	{
	    offset   = nl+1 + i*kdx +iz*nx;
	    pp       = p      + offset;
	    
            hold[i]  = ( ( pp[1+nx].w-pp[-1+nx].w
                          -pp[1-nx].w+pp[-1-nx].w)/4.0
                        -( pp[nx].u-2.*pp[0].u+pp[-nx].u
                          +pp[2*nx].u-2.*pp[nx].u+pp[0].u
                          +pp[1+nx].u-2.*pp[1].u+pp[1-nx].u
                          +pp[1+2*nx].u-2.*pp[1+nx].u+pp[1].u)/4.0)/h/h; 
	    
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* z-derivative */
	{
	    offset   = nl+1 + i*kdx +iz*nx;
	    pp       = p      + offset;

	    hold[i] =  ( (  pp[1+nx].w-2.*pp[1].w+pp[1-nx].w 
	                  -(pp[nx-1].w-2.*pp[-1].w+pp[-1-nx].w)) 
	                -(  pp[2*nx].u-2.*pp[nx].u+pp[0].u + 
	                    pp[2*nx+1].u-2.*pp[nx+1].u+pp[1].u 
	                  -(pp[nx].u-2.*pp[0].u+pp[-nx].u +
	                    pp[nx+1].u-2.*pp[1].u+pp[-nx+1].u) ))/2./h/h/h;
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* (curl curl u)z */
	{
	    offset   = nl+1 + i*kdx +iz*nx;
	    pp       = p      + offset;

	    hold[i] = ( (pp[nx+1].u-pp[nx].u-pp[1].u+pp[0].u)/2.
	               -(pp[1].w-2.*pp[0].w+pp[-1].w))/h/h;
	    
	}
	fwrite(hold,sizeof(float),knx,kirout);
	for(i=0; i< knx; i++)		/* z-derivative */
	{
	    offset   = nl+1 + i*kdx +iz*nx;
	    pp       = p      + offset;

            hold[i] = ( ((pp[2*nx+1].u-pp[2*nx].u-pp[nx+1].u+pp[nx].u)
                        -(pp[1].u-pp[0].u-pp[-nx+1].u+pp[-nx].u))/4. 
                       +((pp[1+nx].w-2.*pp[nx].w+pp[nx-1].w)
                        -(pp[1-nx].w-2.*pp[-nx].w+pp[-nx-1].w))/2.)/h/h/h;
	
	}
	fwrite(hold,sizeof(float),knx,kirout);
    }	
    free(hold);
}
