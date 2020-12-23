/******************************************************************/
/*                                                                */
/* Generalized Ray Theory for Core and Mantle Phases              */
/*                                                                */
/* Written by: Lianxing Wen, Seismological Lab. Caltech           */
/* Last modified: Sunday, Feb. 10, 1997                           */
/*                                                                */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rayconst.h"
#include "getpar.h"

#define RAD  111.194924748
#define DEBUG

int get_refls(int nb, int top, int bottom, int ncmb, int nrec,
              int *lfinal, int **nen, int **ncon, int **na,
	          int **nray, int raytype);
int modify_model(float sdep, int *nb, int *jo, float *c, float *s,
                 float *d, float *th);
void findcmbs(int jo, float *ss, int *ncmb, int *icmb);
int read_model(char *model, int *jo, float **cc, float **ss, float **dd,
	       float **th);

int main(int argc, char *argv[])
{
    char greenfile[80];         // GRT green's library file
    FILE *fop_green;

    char rayfile[80];           // rayfile for tstart
    float tstart;
    float tstart0 = 0.0;
    float reduce_vel = 0.0;

    float rdep   = 2901.0;  /* receiver depth */
    float sdep   = 0.01;    /* source depth */
    int   nb;               /* source layer */
    int   nrec;             /* receiver layer */

    float dist;             /* epicentral distance */
    float dt     = 0.05;
    int   nt     = 2000;
    int   raytype = SV;

    float xmax;
    float dx;
    int   ncmb;             /* the bottomest layer of mantle */
    int   icmb;             /* the bottomest layer of outer core */
    int   top;              /* toppest reflection layer */
    int   bottom;           /* bottomest reflection layer */

    int   isource  = 0;     /* homogeneous source */
    int   line = 1;         /* line source */
    int   iwh  = 0;         /* for GRT */
    int   nstress = 0;
    int   ispecial = 0;
    int   mtd = 3;
    int   iflat = 0;       /* =1 flat 0: no (spherical geometry)*/
    int   skip = 1;

    float *green, *greenz;

    float *ts, *ReadTs(char *file, int nstart, int nx, int dn);

    FILE *open_file(char *file_name, char *access_mode);
    int   grt_();

    /* following for the ray descriptions */
    int lfinal;
    int *nen, *ncoun, *na, *nray;
    int ReadFile = 1;

    int jo;
    float *cc,  *ss,  *dd, *th;
    char model[80];             // 1D reference model

    int i, nr;

    /* source parameters */
    float moment = 1.0;
    float so[5];
    so[0] = so[1] = so[2] = so[3] = so[4] = 1.0;

    setpar(argc, argv);
    mstpar("model",          "s", model);
    /*  Read models from the file model */
    read_model(model, &jo, &cc, &ss, &dd, &th);
    #ifdef DEBUG
    fprintf(stderr, "jo=%d\n", jo);
    #endif

    mstpar("greenfile_WKM",  "s", greenfile);

    mstpar("dt_WKM",         "f", &dt);
    mstpar("xmax_WKM",       "f", &xmax);
    mstpar("dx_WKM",         "f", &dx);
    mstpar("nt_WKM",         "d", &nt);
    mstpar("nr_WKM",         "d", &nr);
    #ifdef DEBUG
    fprintf(stderr, "dt=%f nt=%d dx=%f nr=%d xmax=%f\n",
            dt, nt, dx, nr, xmax);
    fprintf(stderr, "xmin=xmax-(nr-1)*dx=%f\n", xmax-(nr-1)*dx);
    #endif

    getpar("flat",           "d", &iflat);
    getpar("mtd",            "d", &mtd);
    getpar("raytype",        "d", &raytype);
    getpar("skip",           "d", &skip);

    mstpar("ReadFile",       "d", &ReadFile);
    if(ReadFile != 0){
        mstpar("RayParFile",     "s", rayfile);
    } else {
        getpar("Gtstart",        "f", &tstart0);
        getpar("Greduce_vel",    "f", &reduce_vel);
    }

    getpar("rdep",           "f", &rdep);
    getpar("sdep",           "f", &sdep);
    if(sdep < rdep ){
        modify_model(sdep, &nb,   &jo, cc, ss, dd, th);
        modify_model(rdep, &nrec, &jo, cc, ss, dd, th);
    } else {
        modify_model(rdep, &nrec, &jo, cc, ss, dd, th);
        modify_model(sdep, &nb,   &jo, cc, ss, dd, th);
    }
    #ifdef DEBUG
    fprintf(stderr, "flat=%d\n", iflat);
    fprintf(stderr, "sdep=%f nb=%d\n", sdep, nb);
    fprintf(stderr, "rdep=%f nrec=%d\n", rdep, nrec);
    fprintf(stderr, "modified jo=%d\n", jo);
/*  for (i=0; i<jo; i++)
        fprintf(stderr, "%.3d %.3f %.3f %.3f %.3f\n",
                i+1, cc[i], ss[i], dd[i], th[i]);
*/
    #endif

    /* find layer of CMB and ICB */
    findcmbs(jo, ss, &ncmb, &icmb);
    #ifdef DEBUG
    fprintf(stderr, "cmb=%d icb=%d\n", ncmb, icmb);
    #endif

    /* following parameters for core phases only */
    top     = ncmb + 2;
    if (top < nrec + 1)
        top     = nrec + 1;
    bottom  = icmb -1;
    getpar("top",            "d", &top);
    getpar("bottom",         "d", &bottom);

    endpar();

    get_refls(nb, top, bottom, ncmb, nrec, &lfinal, &nen, &ncoun, &na,
             &nray, raytype);
    #ifdef DEBUG
    fprintf(stderr, "lfinal=%d\n", lfinal);
    #endif

    fop_green = open_file(greenfile,"wb");
    green  =(float *) malloc ((nt+3)*sizeof(float));
    greenz =(float *) malloc ((nt+3)*sizeof(float));

    if(ReadFile != 0) ts = ReadTs(rayfile, 0, nr, skip);

    #ifdef DEBUG
    fprintf(stderr, "size of %s: %ld bytes\n", greenfile, nr*nt*sizeof(float)*2);
    fprintf(stderr, "Gstart=%f Gred_vel=%f\n", tstart0, reduce_vel);
    #endif

    for(i = 0; i < nr; i++){    // loop over receivers
        dist   = xmax - i * dx;
        tstart = ((ReadFile != 0) ? ts[i] : (dist/RAD * reduce_vel + tstart0));
        #ifdef DEBUG
        if(i==0||i==nr-1)
            fprintf(stderr, "i=%d dist=%.2f tstart=%.2f\n", i, dist, tstart);
        #endif
        /* if dist<0 green=0 && grenn_z=0 ? */
        grt_(&isource,&iflat,&ispecial,&mtd,&nb,&jo,cc,ss,dd,th,
	        &lfinal,nen,na,nray,ncoun,&line,&iwh,
            &dt,&nt,&dist,&tstart,&nstress,&moment,so,green,greenz);
        fwrite(green, sizeof(float),nt,fop_green);
        fwrite(greenz,sizeof(float),nt,fop_green);
    }

    free(green); free(greenz);
    fclose(fop_green);
    return 1;
}

float *ReadTs(char *Rayfile, int nstart, int nx, int dn)
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
        if(i){
	    /* skip dn -1 Greens points */
	    for(k = 0; k< dn-1; k++)
                fscanf(fop, "%d %f %f %f\n", &n, &dist, &p, &tst);
        }
        fscanf(fop, "%d %f %f %f\n", &n, &dist, &p, &tst);

	ts[i] = tst;

    }

    return ts;
}
