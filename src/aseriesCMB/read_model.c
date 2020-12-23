#include <stdio.h>
#include <stdlib.h>

int read_model(char  *model,     /* model file */
               int   *jo,        /* number of layers */
               float **cc,       /* p-wave velocities */
               float **ss,       /* s-wave velocities */
               float **dd,       /* densities */
               float **th        /* thicknesses */
               )
{
    FILE *fop_model;
    float *c, *s, *d, *tth;
    int j;
    FILE *open_file(char *file_name, char *access_mode);

    fop_model = open_file(model,"r");

    fscanf(fop_model,"%d\n",jo);
    *cc  = (float *) malloc ((*jo+3)*sizeof(float));
    *ss  = (float *) malloc ((*jo+3)*sizeof(float));
    *dd  = (float *) malloc ((*jo+3)*sizeof(float));
    *th  = (float *) malloc ((*jo+3)*sizeof(float));
    c =*cc; s=*ss; d=*dd; tth=*th;
    for(j=0; j<*jo; j++)
        fscanf(fop_model,"%f %f %f %f\n",&c[j],&s[j],&d[j],&tth[j]);
    fclose(fop_model);

    return 1;
}
