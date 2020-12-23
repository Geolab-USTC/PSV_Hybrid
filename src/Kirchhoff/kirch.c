/*************************************************************************/
/*                                                                       */
/*                           kirch.c                                     */
/*                                                                       */
/*     Program for interfacing the waves using Kirchhoff theory          */
/*             Written by: Lianxing Wen                                  */
/*                         Caltech                                       */
/*                                                                       */
/*             Last Modified:  April 23, 1997                            */
/*                                                                       */
/*************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "getpar.h"

#define RAY_LENGTH      2000
#define RAY_NUMBER      200
#define RAD             111.194924748

#define DEBUG

FILE *fdebug;
char debugfile[80];
int debug = 0;

int main(int argc, char *argv[])
{
    /* following for the output from the source */
    char greenfile[80];             /* filename Green's from FD */
    char kirfile[80];               /* filename Green's from receiver*/
    float tstart = 0.0;             /* tstart in the source side*/

    /* following for the resultant outputs */
    char accfile[80];               /* filename for output Accels */
    float *acc;                     /* Resultant X Or Z-Accel */

    /* following for Rayparfile */
    char Rayfile[80];
    int  nstart;
    float  *ts;

    /* general functions */
    float * ReadTs(char *Rayfile, int nstart, int nx, int dn, float dt_reduce);
    float * CreatTs(float tst, float rv, int nx, float dt_reduce);
    void kirchhoff(int, int, float, float *, char *, char *, float *, float );
    void OutputAcc(char *, int, float *,float); /* output accels      */
    float minfloat(int n, float *x);
    FILE *open_file();

    int kdx, knx;
    float fac;
    float reduce_vel = 0.0;
    float dt_reduce = 0.0;

    int   ReadFile = 0;      /* Read Parfile flag   */
    float Gtstart;           /* Tstart for Green's Lib */
    float Greduce_vel=0;     /* Reduce Vel for Green's Lib (degree) */

    int     nx;              /* x-dimension of mesh */
    int     nh;              /* absorbing regions */
    int     nt;              /* number of time steps */
    int     dn;              /* seperation of points */

    float   h;               /* spatial mesh interval */
    float   dt;              /* time digization interval */
    float   xmin;
    float   xmax_WKM;
    float   dist;
    float   dx;

    setpar(argc, argv);
    mstpar("greenfile_kir",   "s", kirfile);
    mstpar("fdkirfile",       "s", greenfile);
    mstpar("ReadFile",        "d", &ReadFile);
    if(ReadFile != 0)
        mstpar("RayParFile",  "s", Rayfile);
    else {
        mstpar("Gtstart",     "f", &Gtstart);
        mstpar("Greduce_vel", "f", &Greduce_vel);
    }
    mstpar("accfile",         "s", accfile);
    getpar("tstart",          "f", &tstart);

    mstpar("h",               "f", &h);
    mstpar("dt_WKM",          "f", &dt);
    mstpar("nt_kir",          "d", &nt);
    mstpar("nx",              "d", &nx);
    mstpar("nh",              "d", &nh);
    mstpar("kdx",             "d", &kdx);
    mstpar("xmin",            "f", &xmin);
    mstpar("xmax_WKM",        "f", &xmax_WKM);
    mstpar("dist",            "f", &dist);
    mstpar("dx_WKM",          "f", &dx);
    getpar("reduce_vel",      "f", &reduce_vel);
    getpar("debug",           "d", &debug);
    if (debug)
        mstpar("debugfile",   "s", debugfile);
    endpar();

    acc     = (float *) malloc (5*nt*sizeof(float));

    knx     = (nx - nh - 3)/kdx;
    dn      = (int) ((h*kdx+0.5*dx)/dx);
    nstart  = (int) ((xmax_WKM -(dist - xmin))/dx);

    if(reduce_vel > 1.e-13)
        dt_reduce = h*kdx/reduce_vel;

    if(ReadFile != 0){
        ts = ReadTs(Rayfile, nstart, knx, dn, dt_reduce);
    } else {
        Gtstart += (dist - xmin)/RAD * Greduce_vel;
	    Greduce_vel *= (-dx * dn /RAD);
        ts = CreatTs(Gtstart, Greduce_vel, knx, dt_reduce);
    }

    /* apply kirchhoff theory */
    fac = kdx*h;
    kirchhoff(knx,nt,dt,ts,greenfile,kirfile,acc,fac);

    tstart += minfloat(knx,ts);

    /* output the accelrations */
    fac = kdx*h;
    OutputAcc(accfile, 5*nt, acc, fac);

    fprintf(stdout, "%f\n",tstart);
    free(acc);
    return 1;
}

float * ReadTs(char *Rayfile, int nstart, int nx, int dn, float dt_reduce)
{
    int    i, k, j;
    float *ts;
    int n; float dist, p, tst;
    FILE  *fop;
    FILE *open_file();

    ts = (float *) malloc (nx *sizeof(float));

    fop = open_file(Rayfile,"r");

    /* skip nstart points */
    for(i=0; i < nstart; i++)
        fscanf(fop, "%d %f %f %f \n", &n, &dist, &p, &tst);

    for(i=0, j=0; i< nx; j += dn, i++){
        if(i > 0){
	    /* skip dn -1 Greens points */
	    for(k = 0; k< dn-1; k++)
            fscanf(fop, "%d %f %f %f\n", &n, &dist, &p, &tst);
        }
        fscanf(fop, "%d %f %f %f\n", &n, &dist, &p, &tst);

	    ts[i] = tst + i * dt_reduce;
    }
    return ts;
}

float *CreatTs(float tst, float rv, int nx, float dt_reduce)
{
    int    i;
    float *ts;

    ts = (float *) malloc (nx *sizeof(float));

    for(i=0; i< nx; i++)
	    ts[i] = tst + i * (rv+dt_reduce);

    return ts;
}

void OutputAcc(char *accfile, int nt, float *ax, float fac)
{
    FILE *fop, *open_file();
    int it;

    fop = open_file(accfile,"w");
    for(it=0; it<nt; it++)
        fprintf(fop, "%e \n",ax[it]*fac);
    fclose(fop);
}

/*************************************************************************/
/*                                                                       */
/*              kirchhoff.c                                              */
/*       kirchhoff integration of two set of Green's function            */
/*                                                                       */
/*       Written by: Lianxing Wen                                        */
/*       Last Modified: April 23, 1997                                   */
/*                                                                       */
/*************************************************************************/

#define TINY   1.e-25

void kirchhoff(int n, int nt, float dt, float *ts, char *greenfile, char *kirfile,
	       float *acc, float fac)
/* This function handles the integration of convolutions of Green's  */
/* functions based on Kirchhoff theory. The Green's functions are    */
/* stored in the greenfile and kirfile. The resultant solution is    */
/* stored in the point acc                                           */
/* Inputs:
    n         -> number of points in the integration line
    nt        -> number of points for Greens function
    dt        -> time step for Greens function
    ts        -> tstarts for kirchhoff green functions
    greenfile -> One of the files where Green's functions are stored
    kirfile   -> other files where Green's functions are stored

    Outputs:
    acc       -> the resultant of Kirchhoff integration
*/
{
    float *green, *green_d;          /* Green's and dG/dz from receiver */
    FILE *fop_recev;

    float *greenx, *greenx_d;        /* Green's and dG/dz from source */
    FILE *fop_source;

    /* general functions */
    void zap(int, float *, float);   /* to initialize a vector */
    void ConvSum(int nt, float dt, float *g1, float *g2, float *a,
		    int nshift, float fac); /* convolve and summation */
    float minfloat(int n, float *x);

    int i, j;
    int fsize;
    FILE *open_file();

    float ts_min;
    float ymax1, ymax2;
    int   nshift;

    int error = 0;

    ts_min =  minfloat(n,ts);

    green     = (float *) malloc (nt*sizeof(float));
    green_d   = (float *) malloc (nt*sizeof(float));
    greenx    = (float *) malloc (nt*sizeof(float));
    greenx_d  = (float *) malloc (nt*sizeof(float));

    zap(5*nt,acc,0.0);

    fsize = sizeof(float);
    fop_source = open_file(greenfile,"rb");
    fop_recev  = open_file(kirfile,  "rb");
    if (debug) {
        fdebug = open_file(debugfile, "wb");
        fwrite(&n, sizeof(int), 1, fdebug);
        fwrite(&nt, sizeof(int), 1, fdebug);
        fwrite(&dt, sizeof(float), 1, fdebug);
    }

    for(i=0; i<n; i++){
        if(error) break;
        ymax1 = ymax2 = 0;
	    nshift = (int) ((ts[i] - ts_min)/dt);

        /* reading the output from source */
        if(fread(greenx,  fsize,nt,fop_source) != nt){
	        fprintf(stderr, "Error reading file %s\n",greenfile);
	        error = 1;
	    }
        if(fread(greenx_d,fsize,nt,fop_source) != nt){
	        fprintf(stderr, "Error reading file %s\n",greenfile);
	        error = 1;
	    }

        if(!error){
	        for(j=0; j<nt; j++)
	            if(fabs(greenx[j]) > ymax1 || fabs(greenx_d[j]) > ymax1)
		        ymax1 = ((fabs(greenx[j]) > fabs(greenx_d[j])) ?
                        fabs(greenx[j]) : fabs(greenx_d[j]));
        }

        /* reading the output from receiver */
        if(fread(green,  fsize,nt,fop_recev) != nt){
	        fprintf(stderr, "Error reading file %s\n",kirfile);
	        error = 1;
	    }
        if(fread(green_d,fsize,nt,fop_recev) != nt){
	        fprintf(stderr, "Error reading file %s\n",kirfile);
	        error = 1;
	    }

        if(!error){
	        for(j=0; j<nt; j++)
	            if(fabs(green[j]) > ymax2 || fabs(green_d[j]) > ymax2)
		            ymax2 = ((fabs(green[j]) > fabs(green_d[j])) ?
                            fabs(green[j]) : fabs(green_d[j]));
        }

        /* If any of the four Green's function is zero,
         * do not sum theme.
         */
        if(!error && ymax1 > TINY  && ymax2 > TINY){
            if (debug) {  /* start time for each trace */
                fwrite(&ts[i], sizeof(float), 1, fdebug);
            }
            ConvSum(nt,dt,green,greenx_d,acc,nshift,  1.0);
            ConvSum(nt,dt,green_d,greenx,acc,nshift, -1.0);
            if (debug) {
                fwrite(acc, sizeof(float), 5*nt, fdebug);
            }
	    }
    }

    free(green);  free(green_d);
    free(greenx); free(greenx_d);
    fclose(fop_source);
    fclose(fop_recev);
}

void ConvSum(int nt, float dt, float *g1, float *g2, float *a, int nshift, float fac)
{
    static int counter = 1;
    int nso  = 1;
    float *convt_out;
    float *convt_in;
//    int i, nt2, nz;
    void sum_vector(int, float *, float *, float, int);
    float maxfloat(int n, float *x);

    int time_domain = 0; /* _flag - time domain convolution */

    convt_out = (float *) malloc(3*nt*sizeof(float));
    convt_in  = (float *) malloc(2*nt*sizeof(float));

    /* do convolution  */
    if( time_domain){
        step_(&nt,&nso,&dt,g1,g2,convt_out);
	    nt     /= 2;
	    nshift /= 2;
    } else {
	    /* frequency domain operations: following from Hiroo Kanamori */
/*	    nt2 = 2*nt;
	    for(i=0; i<nt; i++){
	        convt_in[i+nt] = g1[i];
	        convt_in[i] = 0.0;
        }

        conv_(convt_in,&nt2,g2,&nt,convt_out,&nz);
	    for(i=0; i<nz - nt; i++){
	        j = nt + i;
	        convt_out[i] = convt_out[j] * dt;
        }
*/
        convt1_(g1,&nt,g2,&nt,convt_out,&nso,&dt);
        if (debug) {
            fwrite(convt_out, sizeof(float), nt, fdebug);
        }
    }

    /* do summation  */
    if(fabs(maxfloat(nt,convt_out)) < 1.e20)
        sum_vector(nt,a,convt_out,fac,nshift);
    else
	    fprintf(stderr, "Warning: Bad Greens' functions at point %d(Skip)\n", counter++);

    counter++;

    free(convt_out);
    free(convt_in);
}

/* a = a+b*val */
void sum_vector(int n, float *a, float *b, float val, int nshift)
{
    float *pa=a + nshift, *pb = b;

    if(nshift > 4*n)  n = 5*n - nshift;

    while(n--){
        *pa += *pb * val;
        pa++; pb++;
    }
}

void zap(int n, float *a, float val)
{
    int i;
    for(i=0; i<n; i++) {
        a[i] = val;
    }
}
