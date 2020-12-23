/*
 * Program to convert the aserpsvfd output to sac files
 * Author: Dongdong Tian @ USTC
 * Date: 2014-08-02
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "sacio.h"
#include "getpar.h"

#define RAD 111.194924748

int main(int argc, char *argv[])
{
    char greenfile[80];
    float tstart;
    float xmin;
    float h;
    int nx;
    int nl_skip = 0;
    int kdx;
    int nt;
    FILE *fop;
    int nz;
    float dp;

    setpar(argc, argv);
    mstpar("greenfile", "s", greenfile);
    mstpar("tstart",    "f", &tstart);
    mstpar("xmin",      "f", &xmin);
    mstpar("h",         "f", &h);
    mstpar("nx",        "d", &nx);
    getpar("nl_skip",   "d", &nl_skip);
    mstpar("kdx",       "d", &kdx);
    mstpar("nt",        "d", &nt);
    mstpar("dp",        "f", &dp);
    endpar();

    if ((fop = fopen(greenfile, "rb"))==NULL) {
        fprintf(stderr, "Error in opening file %s\n", greenfile);
        exit(-1);
    }

    nl_skip = ((int)(nl_skip+0.01)/kdx)*kdx;
    nx += nl_skip;

    fread(&nz, sizeof(int), 1, fop);

    SACHEAD hd;
    char sacfile[80];
    char comp[10];
    int n34, istress;
    float dist;
    float *data = NULL;

    data = (float*)malloc(sizeof(float)*nt);

    hd = new_sac_head(dp, nt, tstart);
    for (n34=3; n34<=4; n34++) {
        for (istress=0; istress<=1; istress++) {
            if (n34==4 && istress==0) strcpy(comp, "T12");
            if (n34==4 && istress==1) strcpy(comp, "Uz");
            if (n34==3 && istress==0) strcpy(comp, "Ux");
            if (n34==3 && istress==1) strcpy(comp, "T22");

            int ix;
            for (ix=0; ix<nx; ix+=10) {
                dist = xmin + ix*h + 0.5*istress*h;
                fseek(fop, sizeof(float)*nt*ix, SEEK_SET);
                fread(data, sizeof(float), nt, fop);
                hd.dist  = dist;
                hd.gcarc = dist/RAD;
                sprintf(sacfile, "GRT/%s.%04d.X", comp, ix);
                write_sac(sacfile, hd, data);
            }
        }
    }

    for (n34=3; n34<=4; n34++) {
        dist = ((n34==3) ? xmin : xmin+0.5*h);
        hd.dist = dist;
        hd.gcarc = dist/RAD;
        for (istress=0; istress<=1; istress++) {
            if (n34==4 && istress==0) strcpy(comp, "T11");
            if (n34==4 && istress==1) strcpy(comp, "Uz");
            if (n34==3 && istress==0) strcpy(comp, "Ux");
            if (n34==3 && istress==1) strcpy(comp, "T12");

            int iz;
            for (iz=0; iz<=nz; iz+=10) {  /* iz=[0,nz] */
                fseek(fop, sizeof(float)*nt*iz, SEEK_SET);
                fread(data, sizeof(float), nt, fop);
                sprintf(sacfile, "GRT/%s.%04d.Z", comp, iz);
                write_sac(sacfile, hd, data);
            }
        }
    }
    free(data);
    fclose(fop);
    return 0;
}
