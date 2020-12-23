/*
 * Program to convert the aseriesCMB output to sac files
 * Author: Dongdong Tian @ USTC
 * Date: 2014-08-02
 */

#include <stdio.h>
#include <stdlib.h>
#include "sacio.h"
#include "getpar.h"

#define RAD 111.194924748

int main(int argc, char *argv[])
{
    char WKMfile[80];
    FILE *fop;
    float dt, xmax, dx;
    int nt, nr;
    float Gtstart;
    float Greduce_vel;
    int i;

    setpar(argc, argv);
    mstpar("greenfile_WKM", "s", WKMfile);
    mstpar("dt_WKM",    "f", &dt);
    mstpar("xmax_WKM",  "f", &xmax);
    mstpar("dx_WKM",    "f", &dx);
    mstpar("nt_WKM",    "d", &nt);
    mstpar("nr_WKM",    "d", &nr);
    mstpar("Gtstart",    "f", &Gtstart);
    mstpar("Greduce_vel","f", &Greduce_vel);
    endpar();

    if ((fop=fopen(WKMfile, "rb"))==NULL) {
        fprintf(stderr, "Error in opening Green's file %s\n", WKMfile);
        exit(-1);
    }

    float *green, *greenz;
    green  = (float *)malloc(nt*sizeof(float));
    greenz = (float *)malloc(nt*sizeof(float));

    for (i=0; i<nr; i+=100) {
        float dist, tstart;
        dist = xmax - i *dx;
        tstart = dist/RAD *Greduce_vel + Gtstart;

        fread(green,  sizeof(float), nt, fop);
        fread(greenz, sizeof(float), nt, fop);

        SACHEAD hd;
        char sacfile[80];
        hd = new_sac_head(dt, nt, tstart);
        hd.dist = dist;
        hd.gcarc = dist/RAD;
        sprintf(sacfile, "WKM/%04d.SAC", i);
        write_sac(sacfile, hd, green);
        sprintf(sacfile, "WKM/%04d.SACZ", i);
        write_sac(sacfile, hd, greenz);
    }
    free(green);
    free(greenz);
    fclose(fop);
    return 0;
}
