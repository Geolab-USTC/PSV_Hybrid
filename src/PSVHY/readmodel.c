#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RAY_LENGTH        50000
#define RAY_NUMBER        20000 

FILE *open_file(char *file_name, char *access_mode)
{
    FILE *fop;
 
    if((fop=fopen(file_name,access_mode))== NULL){
        fprintf(stderr,"Error Opening file %s with access %s\n",
                        file_name, access_mode);
        exit(-1);
    }
    return fop;
}                        

int nfinal(int jo,float *ss, int line)
{
    float *s=ss, s1, s2;

    s += line; jo -= line;   /* skip the first line of model */
    s1 = s2 = *s; 
    while(jo--){
        s++; s2 = *s;
        if(fabs(s2-s1) > 0.0001) break;
        s1 = s2;
    }
    return s-ss+1; // layer number start from 1
}

void read_model(int flat, char *model,
        int *jo, int *nb, int *lfinal,
        float **cc, float **ss, float **dd, float **tth,
        int **nen, int **ncoun, int **na, int **nray)
{
    FILE *fop_model;
    float *c, *s, *d, *th;
    int *nen0, *ncoun0, *na0, *nray0;
    int j, k, n;
    double *qrec;
    double flatten(int, float *);

    fop_model= open_file(model,"r");

    /* read the model parameters from "model" */
    fscanf(fop_model,"%d %d\n",jo,nb);
    *cc  = (float *) malloc ((*jo+1)*sizeof(float));
    *ss  = (float *) malloc ((*jo+1)*sizeof(float));
    *dd  = (float *) malloc ((*jo+1)*sizeof(float));
    *tth = (float *) malloc ((*jo+1)*sizeof(float));
    c = *cc; 
    s = *ss;
    d = *dd; 
    th= *tth;
    for(j=0; j<*jo; j++)
        fscanf(fop_model,"%f %f %f %f\n",&c[j],&s[j],&d[j],&th[j]);

    if(flat == 0){
	    /* flattening the earth */
	    qrec  = (double *) malloc ((*jo+1)*sizeof(double));
	    for(j=0; j < *jo; j++)
	        qrec[j] = flatten(j,th);

	    for(j=0; j < *jo; j++){
	        c[j]  *= qrec[j];
	        s[j]  *= qrec[j];
	        d[j]  *= qrec[j];
	        th[j] *= qrec[j];
        }
        free(qrec);
    }

    /* read ray parameters for the rays at n=3,4 */
    fscanf(fop_model,"%d\n",lfinal);
    *nen   = (int *) malloc (RAY_NUMBER * sizeof(int));
    *ncoun = (int *) malloc (RAY_NUMBER * sizeof(int));
    *na    = (int *) malloc (RAY_LENGTH * sizeof(int));
    *nray  = (int *) malloc (RAY_LENGTH * sizeof(int));

    nen0=*nen; ncoun0=*ncoun; na0=*na; nray0=*nray;
    
    k=0;
    for(n=0; n<*lfinal; n++){
        fscanf(fop_model,"%d",&nen0[n]);
        for(j=0; j<nen0[n]; j++)
            fscanf(fop_model," %d",&na0[k++]);
        fscanf(fop_model,"\n");     
        k -= nen0[n];

        fscanf(fop_model,"%d",&ncoun0[n]);
        for(j=0; j<nen0[n]; j++)
            fscanf(fop_model," %d",&nray0[k++]);
        fscanf(fop_model,"\n");     
    }
    fclose(fop_model);
}

#define EARTH_RADIUS     6371.0 
double flatten(int layer, float *tth)
{
    int i;
    double thickness = 0.0;

    for(i=0; i<layer; i++)
        thickness += tth[i];
    thickness += 0.5*tth[layer]; 

    return (EARTH_RADIUS/(EARTH_RADIUS-thickness));
}
