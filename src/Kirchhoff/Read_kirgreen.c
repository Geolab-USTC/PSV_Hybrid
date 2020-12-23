/*************************************************************************/
/*                                                                       */
/*                          READ_kirgreen.c                              */
/*                                                                       */
/*     Program for reading Green's functions from receivers              */
/*             Written by: Lianxing Wen                                  */
/*                         Caltech                                       */
/*                                                                       */
/*             Last Modified:  Nov 19, 1996                              */
/*                                                                       */
/*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "getpar.h"

#define DEBUG

int nr_WKM;
int main(int argc, char *argv[])
{
    char WKMfile[80];
    char kirfile[80];

    int     nstart;          /* start element */
    int     dn;              /* increment of the elements */
    int     nx;              /* number of meshs */
    int     nh;              /* number of meshs of absorbing zone */
    int     nt_WKM;          /* number of time steps in WKM*/
    int     nt_kir;          /* number of time steps in kirchhoff*/
    float   dt;              /* time step for the library */
    float   h;               /* space interval in FD */
    int     kdx;             /* spacing between receivers in FD */
    float   dx;              /* space interval in WKM */
    float   dist;            /* epicentral distance for the receiver */
    float   xmin;            /* left boundary in FD */
    float   xmax_WKM;        /* epicentral distance of leftest point in WKM*/
    int     taper  = 0;      /* taper or not */
    float   xcenter;
    float   xleft;
    float   xright;
    float   xtaper;
    int     ncenter = 0;
    int     nleft   = 0;
    int     nright  = 0;
    int     ntaper  = 0;

    int Read_kirgreen(char *WKMfile, char *kirfile, int nt_WKM,
                      int nt_kir, int nx, float dt, int nstart, int dn,
                      int ncenter, int nleft, int nright, int ntaper);

    setpar(argc,argv);
    mstpar("greenfile_WKM",   "s", WKMfile);
    mstpar("greenfile_kir",   "s", kirfile);
    mstpar("nt_kir",          "d", &nt_kir);
    mstpar("nt_WKM",          "d", &nt_WKM);
    mstpar("nx",              "d", &nx);
    mstpar("nh",              "d", &nh);
    mstpar("kdx",             "d", &kdx);
    mstpar("dt_WKM",          "f", &dt);
    mstpar("h",               "f", &h);
    mstpar("dx_WKM",          "f", &dx);
    mstpar("nr_WKM",          "d", &nr_WKM);
    mstpar("xmax_WKM",        "f", &xmax_WKM);
    mstpar("xmin",            "f", &xmin);
    mstpar("dist",            "f", &dist);
    getpar("taper",           "d", &taper);
    if(taper > 0){
        mstpar("xcenter",          "f", &xcenter);
        mstpar("xleft",            "f", &xleft);
        mstpar("xright",           "f", &xright);
        mstpar("xtaper",           "f", &xtaper);
    }
    endpar();

    nx        = (nx -nh -3)/kdx;
    dn        = (int) ((h*kdx + 0.5*dx)/dx);
    // nstart is the most distant record of GRT library
    // dist = xmax_WKM - i*dx
    nstart    = (int) (xmax_WKM - (dist - xmin))/dx;
    fprintf(stderr,"nx=%d dn=%d nstart=%d nt_WKM=%d\n",nx,dn,nstart,nt_WKM);
    if(taper > 0){
        ncenter   = (int) ((xmax_WKM -xcenter)/dx );
        nleft     = (int) ((xmax_WKM -xleft)/dx );
        nright    = (int) ((xmax_WKM -xright)/dx );
        ntaper    = (int) (xtaper/dx);
    }
    #ifdef DEBUG
    fprintf(stderr, "xcenter=%f xleft=%f xright=%f xtaper=%f\n",
                xcenter, xleft, xright, xtaper);
    fprintf(stderr, "ncenter=%d nleft=%d nright=%d ntaper=%d\n",
                ncenter, nleft, nright, ntaper);
    #endif

    Read_kirgreen(WKMfile, kirfile,  nt_WKM, nt_kir,
              nx, dt, nstart, dn,ncenter,nleft,nright,ntaper);
    return 1;
}

int Read_kirgreen(char *WKMfile, char *kirfile,  int nt_WKM,
                  int nt_kir, int nx, float dt, int nstart, int dn, int ncenter,
                  int nleft, int nright, int ntaper)
{
    FILE *open_file();
    FILE *fop_WKM;
    FILE *fop_kir;

    float *greenx,*greenx_d;
    int i, j, k;
    int nt;
    float factor;

    fop_WKM   = open_file(WKMfile,"rb");

    if((fop_kir   = open_file(kirfile,"wb+"))==NULL) {
        fprintf(stderr,"open error\n");
        return -1;
    }

    nt = (nt_WKM < nt_kir) ? nt_kir : nt_WKM;
    greenx    = (float *) malloc (nt*sizeof(float));
    greenx_d  = (float *) malloc (nt*sizeof(float));

    for(i=nt_WKM; i<nt; i++){
        greenx[i] = greenx_d[i] = 0.0;
    }

    if(fseek(fop_WKM, 2*nstart*nt_WKM*sizeof(float), SEEK_SET) != 0){
        fprintf(stderr, "fseek error exiting ... \n");
        exit(-1);
    }

    /* starting to read the green's function*/
    for(i=0, j=0; i<nx; i++, j += dn){
        /* skip dn-1 Greens functions */
        if(i > 0) fseek(fop_WKM, sizeof(float)*nt_WKM*(dn-1)*2, SEEK_CUR);

        if (j + nstart < nr_WKM) {
            /* Reading the Green's functions  */
            if(fread(greenx  ,sizeof(float),nt_WKM,fop_WKM) != nt_WKM ){
                fprintf(stderr, "3 Error in reading the file %s at %d\n", WKMfile,j);
                fclose(fop_kir);
                exit(-1);
            }
            if(fread(greenx_d ,sizeof(float),nt_WKM,fop_WKM) != nt_WKM ){
                fprintf(stderr, "4 Error in reading the file %s at %d\n", WKMfile,j);
                fclose(fop_kir);
                exit(-1);
            }
        } else {    // no Green's functions to read
            int it;
            for (it=0; it<nt; it++) {
                greenx[it]   = 0.0;
                greenx_d[it] = 0.0;
            }
        }

        int index = j + nstart;
        if(ncenter > 0){
            factor = 1.0;
            if( index < nleft-ntaper || index > nright+ntaper){
                factor = 0.0;
            } else if (index >= nleft-ntaper && index < nleft){
                factor = (float)(index - (nleft - ntaper)) / ntaper;
            } else if (index > nright && index <= nright + ntaper){
                factor = (float)(nright + ntaper - index) / ntaper;
            }
            for(k=0; k<nt_kir; k++){
                greenx[k] *= factor;
                greenx_d[k] *= factor;
            }
        }
        #ifdef DEBUG
        if ((index == nleft-ntaper) || (index==nleft)
                || (index==nright) || (index==nright+ntaper)) {
            fprintf(stderr, "j=%d index=%d factor=%f\n", j, j+nstart, factor);
        }
        #endif
        fwrite(greenx  ,sizeof(float),nt_kir,fop_kir);
        fwrite(greenx_d,sizeof(float),nt_kir,fop_kir);
    }
    free(greenx); free(greenx_d);
    fclose(fop_kir);
    fclose(fop_WKM);

    return 1;
}
