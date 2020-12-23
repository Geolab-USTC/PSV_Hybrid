#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define NP_MAX 1000
#define XMIN -1.0
#define XMAX 8000.0
#define EARTH_RADIUS 6371.0
#define ICB 5149.5

int main(int argc, char *argv[])
{
    float xstart;         /* start for topos */
    float zbase;          /* base depth of topos relative to GRT-FD interface */
    int num;              /* number of topos */
    int div;              /* distance between each two topo */
    float height, width;  /* control the shape of topo */
    int np;               /* half number of points for a topo */
    float x[NP_MAX], y[NP_MAX];

    int i, j, k;
    float xmov, xstep;
    int layers;

    float vp0, vs0, d0;
    float topovp, topovs, topod;
    float dvp, dvs, dd;
    float thick, high;

    float dep;

    float trans;        /* thickness of transition zone below ICB before flatten */
    float dbase;        /* depth increment of zbase for transition zone before flatten */

    /* input */
    scanf("%f %f", &xstart, &zbase);
    scanf("%f %f", &height, &width);
    scanf("%d %d", &num, &div);
    scanf("%d", &np);
    scanf("%f %f %f", &vp0, &vs0, &d0);
    scanf("%f %f %f", &topovp, &topovs, &topod);
    scanf("%d", &layers);
    scanf("%f", &trans);

    float vp[100], vs[100], den[100];

    if (trans>0.1) { // case 3
        int i;
        for (i=0; i<layers+1; i++)
            scanf("%f %f %f", &vp[i], &vs[i], &den[i]);
    } else {    // case 1 or case 2
        /* parameters increment */
        dvp = (topovp - vp0) / layers;
        dvs = (topovs - vs0) / layers;
        dd  = (topod  - d0 ) / layers;

        for (i=0; i<layers+1; i++) {
            vp[i] = vp0 + dvp * i;
            vs[i] = vs0 + dvs * i;
            den[i] = d0  + dd  * i;
        }
    }

    xstep = width / (2*np);
    thick = height/layers;  // thickness for each layer
    dbase = trans/layers;   // zbase increment

    fprintf(stdout, "%d\n", layers + 1);

    int cnt;
    for (cnt=0; cnt<layers; cnt++) {
        high = height - cnt * thick;    /* height of topo for each layer */

        k = 0;
        for (i=0; i<num; i++) {
            float x0;
            x0 = xstart + (div + width) * i;  /* start position for each topo */
            for (j=0; j<2*np + 1; j++) {
                xmov = xstep * j;
                x[k] = x0 + xmov;
                y[k] = -high*exp(-pow(((xmov-width/2.0)/(width/5.0)), 2.0));
                y[k] += zbase;
                k++;
            }
        }
        if (k>NP_MAX) {
            fprintf(stderr, "Too many points.\n");
            exit(0);
        }

        // parameters
        fprintf(stdout, "%d %f %f %f\n", (2*np+1)*num+2, vp[cnt], vs[cnt], den[cnt]);

        // X array
        fprintf(stdout, "%f ", XMIN);
        for (i=0; i<(2*np+1)*num; i++)
            fprintf(stdout, "%f ", x[i]);
        fprintf(stdout, "%f ", XMAX);
        fprintf(stdout, "\n");

        // Y array
        fprintf(stdout, "%f ", zbase);
        for (i=0; i<(2*np+1)*num; i++)
            fprintf(stdout, "%f ", y[i]);
        fprintf(stdout, "%f ", zbase);
        fprintf(stdout, "\n");

        dep = ICB + (cnt + 0.5)* dbase;     /* depth before flatten */

        zbase += dbase * EARTH_RADIUS/(EARTH_RADIUS-dep);   /* PLUS not MINUS */
        fprintf(stderr, "%f\n", zbase);
    }

    fprintf(stdout, "%d %f %f %f\n", 2, vp[layers], vs[layers], den[layers]);
    fprintf(stdout, "%f %f \n", XMIN, XMAX);
    fprintf(stdout, "%f %f \n", zbase, zbase);

    return 0;
}
