/*
 * Program to convert the demult2kir output to sac files
 * Author: Dongdong Tian @ USTC
 * Date: 2014-08-03
 */
#include <stdio.h>
#include <stdlib.h>
#include "getpar.h"
#include "sacio.h"

int main(int argc, char *argv[])
{
    int nx, nh, kdx, nt_kir;
    char kirfile[80];
    FILE *fop;
    float dt;
    char dir[80];

    setpar(argc, argv);
    mstpar("nx",    "d", &nx);
    mstpar("nh",    "d", &nh);
    mstpar("kdx",   "d", &kdx);
    mstpar("nt_kir","d", &nt_kir);
    mstpar("kirfile","s", kirfile); // kirfile_x or kirfile_z
    mstpar("dir",    "s", dir);
    mstpar("dt",    "f", &dt);
    endpar();

    int ntr;
    ntr = (int)((nx-nh-3+0.01)/kdx);

    int i, m;
    if ((fop = fopen(kirfile, "rb"))==NULL){
        fprintf(stderr, "Error in opening file %s\n", kirfile);
        exit(-1);
    }

    float *data;
    char sacfile[80];
    SACHEAD hd;
    hd = new_sac_head(dt, nt_kir, 0.0);
    data = (float *)malloc(sizeof(float)*nt_kir);
    for (i=0; i<ntr; i++) {
        for (m=0; m<2; m++) {
            if (m==0) sprintf(sacfile, "%s/%04d.SAC", dir, i);
            if (m==1) sprintf(sacfile, "%s/%04d.SACZ", dir, i);
            fread(data, sizeof(float), nt_kir, fop);
            write_sac(sacfile, hd, data);
        }
    }
    free(data);
    fclose(fop);
    return 0;
}
