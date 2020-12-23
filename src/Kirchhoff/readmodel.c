#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define RAY_LENGTH        2000
#define RAY_NUMBER        200 

FILE *open_file(file_name,access_mode)
char *file_name, *access_mode;

{
    FILE *fop;
 
    if((fop=fopen(file_name,access_mode))== NULL){
        fprintf(stderr,"Error Opening file %s with access %s\n",
                        file_name, access_mode);
	perror("jjj");
        exit(-1);
    }
    return fop;
}                        

int nfinal(jo,ss)
int jo; float *ss;

{
    float *s=ss;

    while(jo--)
      if(*s++<0.01)break;

    return s-ss;
}

void read_model(model,jo,nb,lfinal,cc,ss,dd,tth,nen,ncoun,na,nray)
char *model;
int *jo, *nb; float **cc, **ss, **dd, **tth;
int *lfinal, **nen, **ncoun, **na, **nray;
{

    FILE * fop_model;
    float *c, *s, *d, *th;
    int *nen0, *ncoun0, *na0, *nray0;
    int j, k, n;

    fop_model= open_file(model,"r");

/* read the model parameters from "model" */
    fscanf(fop_model,"%d %d\n",jo,nb);
    *cc  = (float *) malloc ((*jo+1)*sizeof(float));
    *ss  = (float *) malloc ((*jo+1)*sizeof(float));
    *dd  = (float *) malloc ((*jo+1)*sizeof(float));
    *tth = (float *) malloc ((*jo+1)*sizeof(float));
    c =*cc; s=*ss; d=*dd; th=*tth;
    for(j=0; j<*jo; j++)
        fscanf(fop_model,"%f %f %f %f\n",&c[j],&s[j],&d[j],&th[j]);

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


#ifdef SUN 
#define SEEK_CUR 1
#define SEEK_SET 0
#endif

int read_record(fop,offset,u,num,size,rec_len)
FILE *fop; int size, num; long int offset, rec_len;
float *u;

{
     float *pu=u;
     int  rsize; 
     
     if(fseek(fop,offset,SEEK_SET)!=0){
         fprintf(stderr,"fseek Error exiting......\n"); 
         return (-1);
     }  
     if((rsize=fread(&pu[0],size,1,fop))!=1){
         fprintf(stderr,"Error in reading the file\n"); 
         return (-1);
     }    
     num--;
     
     while(num--){
         pu++; offset += rec_len;
         if(fseek(fop,offset,SEEK_SET)!=0){
             fprintf(stderr,"fseek Error exiting......\n"); 
             return (-1);
         }  
         if((rsize=fread(&pu[0],size,1,fop))!=1){
             fprintf(stderr,"Error in reading the file\n"); 
              return (-1);
         }    
     }
     return (1);
}     

