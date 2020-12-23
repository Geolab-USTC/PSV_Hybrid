#include <stdio.h>
#include <stdlib.h>
#include "getpar.h"
#include "sacio.h"

int main(int argc, char *argv[])
{
    char debugfile[80];
    char dir[80];
    FILE *fop;
    int n, nt;
    float dt;

    setpar(argc, argv);
    mstpar("debugfile", "s", debugfile);
    mstpar("dir",       "s", dir);
    endpar();

    fop = fopen(debugfile, "rb");
    fread(&n,  sizeof(int), 1, fop);
    fread(&nt, sizeof(int), 1, fop);
    fread(&dt, sizeof(float), 1, fop);
    int i;
    float *ts;
    float *conv, *acc;
    char sacfile[80];
    SACHEAD hd;
    ts = (float *)malloc(sizeof(float)*n);
    conv = (float *)malloc(sizeof(float)*nt);
    acc  = (float *)malloc(sizeof(float)*nt*5);
    for (i=0; i<n; i++) {
        if (fread(&ts[i], sizeof(float), 1, fop)!=1)
            break;
        hd = new_sac_head(dt, nt, ts[i]);

        fread(conv, sizeof(float), nt, fop);
        sprintf(sacfile, "%s/%04d.SAC1", dir, i);
        write_sac(sacfile, hd, conv);

        fread(conv, sizeof(float), nt, fop);
        sprintf(sacfile, "%s/%04d.SAC2", dir, i);
        write_sac(sacfile, hd, conv);

        hd = new_sac_head(dt, 5*nt, ts[i]);
        fread(acc, sizeof(float), 5*nt, fop);
        sprintf(sacfile, "%s/%04d.acc", dir, i);
        write_sac(sacfile, hd, acc);
    }
    free(ts);
    free(acc);
    free(conv);
    fclose(fop);

    return 0;
}
