#include <stdio.h>
#include <stdlib.h>
#include "getpar.h"
#include "sacio.h"

int main(int argc, char *argv[]) {
    char tracefile[80];
    FILE *ftrace;
    char dir[80];

    setpar(argc, argv);
    mstpar("trace_file", "s", tracefile);
    mstpar("dir",        "s", dir);
    endpar();

    if ((ftrace = fopen(tracefile, "rb"))==NULL) {
        fprintf(stderr, "Error in opening file %s\n", tracefile);
        exit(-1);
    }

    int xtrace, ztrace;
    int nt;
    float dt, tstart;
    char ufile[80], wfile[80];
    SACHEAD hdu, hdw;
    float *u, *w;

    fread(&xtrace,  sizeof(int), 1, ftrace);
    fread(&ztrace,  sizeof(int), 1, ftrace);
    fread(&nt,      sizeof(int), 1, ftrace);
    fread(&dt,      sizeof(float),1, ftrace);
    fread(&tstart,  sizeof(float),1, ftrace);

    fprintf(stderr, "%d %d %d %f %f\n", xtrace, ztrace, nt, dt, tstart);

    size_t head = 3*sizeof(int) + 2*sizeof(float);
    size_t per_record = 2*sizeof(float);
    size_t per_step = per_record*xtrace*ztrace;

    u = (float *)malloc(nt*sizeof(float));
    w = (float *)malloc(nt*sizeof(float));

    int ix, iz, it, nn;
    size_t offset;

    for (ix=0; ix<xtrace; ix++) {
        for (iz=0; iz<ztrace; iz++) {
            nn = iz*xtrace+ix;
            sprintf(ufile, "%s/%03d.%03d.U", dir, ix, iz);
            sprintf(wfile, "%s/%03d.%03d.W", dir, ix, iz);
            hdu = new_sac_head(dt, nt, tstart);
            hdw = new_sac_head(dt, nt, tstart);

            for (it=0; it<nt; it++) {
                offset =  it*per_step + nn*per_record + head;
                if (fseek(ftrace, offset, SEEK_SET) != 0) {
                    fprintf(stderr, "fseek error!\n");
                    exit(-1);
                }
                fread(&u[it], sizeof(float), 1, ftrace);
                fread(&w[it], sizeof(float), 1, ftrace);
            }
            write_sac(ufile, hdu, u);
            write_sac(wfile, hdw, w);
        }
    }
    fclose(ftrace);

    return 0;
}
