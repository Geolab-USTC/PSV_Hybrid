#include <stdio.h>
#include <stdlib.h>
#include "sacio.h"
#include "getpar.h"

int main(int argc, char *argv[])
{
    int nx, nh, kdx;
    float h, dx;
    char kirfile[80];
    char dir[80];
    FILE *fop;
    float nt_kir, nt_WKM, dt;

    setpar(argc, argv);
    mstpar("nx",    "d",    &nx);
    mstpar("nh",    "d",    &nh);
    mstpar("kdx",   "d",    &kdx);
    mstpar("h",     "f",    &h);
    mstpar("dt",    "f",   &dt);
    mstpar("dx_WKM","f",    &dx);
    mstpar("nt_WKM","f",    &nt_WKM);
    mstpar("nt_kir","f",    &nt_kir);
    mstpar("greenfile_kir", "s", kirfile);
    mstpar("dir",   "s",    dir);
    endpar();

    nx = (nx-nh-3)/kdx;

    if ((fop = fopen(kirfile, "rb"))==NULL) {
        fprintf(stderr, "Error in opening file %s\n", kirfile);
        exit(-1);
    }

    int i;
    float *green, *green_d;
    SACHEAD hd;
    char sacfile[80];
    hd = new_sac_head(dt, nt_kir, 0.0);
    green  = (float *)malloc(sizeof(float)*nt_kir);
    green_d= (float *)malloc(sizeof(float)*nt_kir);
    for (i=0; i<nx; i++) {
        fread(green,  sizeof(float), nt_kir, fop);
        sprintf(sacfile, "%s/%04d.SAC", dir, i);
        write_sac(sacfile, hd, green);

        fread(green_d,sizeof(float), nt_kir, fop);
        sprintf(sacfile, "%s/%04d.SACZ", dir, i);
        write_sac(sacfile, hd, green_d);
    }
    free(green);
    free(green_d);
    fclose(fop);
    return 0;
}
