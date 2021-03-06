/******************************************************************/
/*                                                                */
/* GRT for P-SV system                                            */
/*                                                                */
/* Written by: Lianxing Wen, Seismological Lab. Caltech           */
/* Last modified: Sunday, Feb. 10, 1997                           */
/*                                                                */
/******************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "getpar.h"
#include "nrutil.h"

#define DEBUG

#define RAY_LENGTH      1000000
#define RAY_NUMBER      40000

#define DIRECT          0
#define REFLECTED       1
#define TOTAL           2
#define DOWN           -1
#define UP              1

int     nx;              /* x-dimension of mesh */
int     nz;              /* z-dimension of mesh */
int     nd;              /* z-dimension of mesh */
int     nf         =10;  /* number of fluid meshs */
int     nl_skip = 0;     /* number of meshs skip in the left output */
int     nt;              /* number of time steps */

float   h;               /* spatial mesh interval */
float   dp;              /* time digization interval */
float   xmin;            /* the left edge of the FD region */

int     layerrefl = 0;   /* add the reflection of layers (other than CMB */
int     convert   = 1;   /* ignore convert phaser or not */

/* following for the GRT */
float tstart = 0.0;
char raymodel[80], greenfile[80];

/* following for the Gauss sources */
float sgm, ts, *source;
float theta, dip, lamda, azmuth, azmuth1;

void getrays(int lfinal, int *nen, int *na, int *nray, int *ncoun,
        int zlay, int nrec, int nfil, int jo, int nb,
        int *lfinal0, int *nen0, int *ncoun0, int *na0, int *nray0,
        int jo0, int *nb0, int bounce, int raytype);
int BounceRays(int nbounce, int **rtype, int P_SV, int SV_P);
int MultiBounces(int nrec, int nfil, int zlay, int k, int bounce,int *type,
        int *na, int *nray, int *ncoun, int *nen);
int power(int a, int n);

int main(int argc, char *argv[])
{
    FILE *fop_green, *open_file();
    int n34, istress, nstress, isource;
    int ispecial;
    int nrec, nfil, zlay, nfil0, zlaymax;
    int i, n;
    float  xx, thickness, dep;
    int nso, nconv, l2n, nx1, ncent, fscl;
    float *so0, *so, *dso, *so1;
    float *green;

    /* following for the model */
    int jo, nb, jo0, nb0;
    float *cc,  *ss,  *dd, *tth;
    float *cc0, *ss0, *dd0,*tth0;
    /* following for the ray descriptions */
    int lfinal, lfinal0;
    int *nen, *ncoun, *na, *nray;
    int *nen0,*ncoun0,*na0,*nray0;
    int bounce = 3;
    int raytype;
    int kdx = 1;
    int mtd = 3;

    int der = 0;
    int iflat = 1;     /* =1 flat 0: no (spherical geometry)*/

    int grt_(), gauss_(), convt0_(), convt_(), fault_();
    void read_model(int flat, char *model,
            int *jo, int *nb, int *lfinal,
            float **cc, float **ss, float **dd, float **tth,
            int **nen, int **ncoun, int **na, int **nray);
    int nfinal(int jo,float *ss, int line);
    int modify_model(float dep, int nrec, int nfil, int jo,
            float *cc, float *ss, float *dd, float *tth,
            int *jo0, float *cc0, float *ss0, float *dd0, float *tth0);

    setpar(argc,argv);
    mstpar("raymodel",  "s", raymodel);
    mstpar("greenfile", "s", greenfile);
    getpar("tstart",    "f", &tstart);

    mstpar("h",         "f", &h);
    mstpar("dp",        "f", &dp);
    mstpar("nt",        "d", &nt);
    mstpar("nx",        "d", &nx);
    mstpar("nf",        "d", &nf);
    getpar("nlskip",    "d", &nl_skip);
    getpar("kdx",       "d", &kdx);
    getpar("mtd",       "d", &mtd);
    getpar("der",       "d", &der);
    mstpar("nd",        "d", &nd);
    mstpar("xmin",      "f", &xmin);

    getpar("flat",      "d", &iflat);
    getpar("layerrefl", "d", &layerrefl);
    getpar("convert",   "d", &convert);
    getpar("bounce",    "d", &bounce);

    mstpar("sgm",       "f", &sgm);
    mstpar("ts",        "f", &ts);
    mstpar("theta",     "f", &theta);
    mstpar("dip",       "f", &dip);
    mstpar("lamda",     "f", &lamda);
    mstpar("azmuth",    "f", &azmuth);
    endpar();

    /* calculate the radiation patterns */
    so0 = (float *)malloc(5*sizeof(float));
    fault_(&theta,&dip,&lamda,&azmuth,so0);

    azmuth1 = azmuth + 180;
    so1 = (float *)malloc(5*sizeof(float));
    fault_(&theta,&dip,&lamda,&azmuth1,so1);

    #ifdef DEBUG
    fprintf(stderr, "Radation patter for az=%f: %f %f %f %f %f\n",
            azmuth ,so0[0],so0[1],so0[2],so0[3],so0[4]);
    fprintf(stderr, "Radation patter for az=%f: %f %f %f %f %f\n",
            azmuth1,so1[0],so1[1],so1[2],so1[3],so1[4]);
    #endif

    /* FFT the Gauss source function */
    nso = (int) (2.*ts/dp);     // number of points of source function
    source =(float *) malloc (nso*sizeof(float));
    gauss_(&sgm,&dp,&ts,source,&nso);
    #ifdef DEBUG
    FILE *fgauss;
    fgauss = fopen("gauss.sac", "w");
    for (i=0; i<nso; i++)
        fprintf(fgauss, "%f\n", source[i]);
    fclose(fgauss);
    #endif

    /* what's this? */
    dso =(float *) malloc (80400*sizeof(float));
    convt0_(&nt,&nso,source,dso,&dp,&nconv,&l2n,&nx1,&ncent,&fscl);
    free(source);

    /* modify the left boundary */
    nl_skip = ((int)((nl_skip+0.001)/kdx))* kdx;
    nx   += nl_skip;
    xmin -= nl_skip *h;

    /*
     * Read models from the file model
     *
     * jo: total number of layers
     * nb: source layer
     * lfinal: total number of rays
     * cc, ss, dd, tth: Vp, Vs, density, thickness
     * nen: number of segments for each ray
     * ncoun: type of each ray, UP or DOWN
     * na: ray path for each ray
     * nray: wave type of each segents for rays
     */
    #ifdef DEBUG
    fprintf(stderr, "Read raymodel from %s...\n", raymodel);
    #endif
    read_model(iflat,raymodel,&jo,&nb,&lfinal,&cc,&ss,&dd,&tth,&nen,&ncoun,&na,&nray);

    for(n=0; n<lfinal; n++)
        ncoun[n]=((na[nen[0]-2]<=nb) ? UP : DOWN);  // why 0 ?
        // nen-2 is the layer of last non-1 segment

    #ifdef DEBUG
    fprintf(stderr, "total number of layer = %d\n", jo);
    fprintf(stderr, "source layer = %d\n", nb);
    fprintf(stderr, "number of rays = %d\n", lfinal);
    /*
    int di;
    fprintf(stderr, "Reference Model:\n");
    for(di=0;di<jo;di++) {
        // count from 1 ?
        // n represent the top interface of layer n ?
        fprintf(stderr, "%03d %.3f %.3f %.3f %.3f\n",
                di+1, cc[di], ss[di], dd[di], tth[di]);
    }

    for(di=0; di<lfinal; di++) {
        fprintf(stderr, "Ray %d: %d segments %d->%d\n",
                di, nen[di], na[0], na[nen[di]-2]);
    }
    */
    #endif

    /* The layer number where n=3,4 are, the layer of the interface*/
    nrec = na[nen[0]-2];    // the layer of the interface
    /* Nfil is the layer number where fluid or free surface is */
    nfil = nfinal(jo,ss,nrec);

    #ifdef DEBUG
    fprintf(stderr, "receiver layer=%d\n", nrec);
    fprintf(stderr, "ICB layer=%d\n", nfil);
    #endif

    /* Evaluate the thickness of the FD region; modify the thickness
       of the last layer, or order to fit the actual region into the
       FD grids, nz is the number of total HALF grids in the solid
    */
    // nrec is the bottom interface of the layer;
    // nfil is the top interface of the layer
    thickness=0.0;
    if(nfil>nrec){
        for(n=nrec; n<nfil-1; n++)
	        thickness += tth[n];
    } else {
        for(n=nfil; n<nrec-1; n++)
	        thickness += tth[n];
    }

    /* modify the thickness of last layer above ICB to fit the actual FD region
     * so the ICB may be deeper than the reference model.
     * The depth of ICB has a error less than 0.5*h
     */
    nz = (int) ((thickness+1.e-3)/(0.5*h));
    n =  ((nfil>nrec) ? nfil-2 : nfil);
    tth[n]= 0.5*nz*h- (thickness-tth[n])+0.0001+0.5*h;
    thickness = nz*0.5*h + 0.5 *h;

    n  = (int) (tth[nrec-1]/h) - 3;
    if( nd > n ) nd = n;

    #ifdef DEBUG
    fprintf(stderr,"!!!!!  thickness of FD domain = %f !!!!\n",thickness);
    fprintf(stderr,"number of half grid: nz=%d\n", nz);
    #endif

    /* This set of points are derived from the initial model and rays */
    cc0   = (float *) malloc ((jo+1)*sizeof(float));
    ss0   = (float *) malloc ((jo+1)*sizeof(float));
    dd0   = (float *) malloc ((jo+1)*sizeof(float));
    tth0  = (float *) malloc ((jo+1)*sizeof(float));

    nen0   = (int *) malloc (RAY_NUMBER * sizeof(int));
    ncoun0 = (int *) malloc (RAY_NUMBER * sizeof(int));
    na0    = (int *) malloc (RAY_LENGTH * sizeof(int));
    nray0  = (int *) malloc (RAY_LENGTH * sizeof(int));

    nz = (nz+1+2*nf)/2-1 + nd;
    // nz: half grids from receiver to ICB
    // nf: grid below ICB
    // nd: reflection layer above receiver
    #ifdef DEBUG
    fprintf(stderr, "number of z grids: nf=%d nd=%d nz=%d\n", nf, nd, nz);
    #endif

    /* output Green's function */
    fop_green = open_file(greenfile,"wb");
    green =(float *) malloc ((nt+3)*sizeof(float));
    fwrite(&nz,sizeof(int),1,fop_green);

    #ifdef DEBUG
    long int size = (nx+nz+1)*4*sizeof(float);
    fprintf(stderr,"nt=%d\n", nt);
    fprintf(stderr,"size of %s = %ld\n", greenfile, sizeof(int)+size*nt);
    #endif

    nz -= nd;

    dep = (nz -1.0)*h;
    zlaymax = modify_model(dep,nrec,nfil,jo,cc,ss,dd,tth,&jo0,cc0,ss0,dd0,tth0);
    zlaymax -= 1;

    /* Calculate the initial responses at n=3 and 4 */
    isource=1;
    iflat = 1;   /* take care of GRT calculation:model has been flattened*/
    ispecial = 0;
    for(n34=3; n34<=4; n34++){
	    dep = ((n34 == 3) ? -0.75*h : -0.25 *h);
        zlay = modify_model(dep,nrec,nfil,jo,cc,ss,dd,tth,&jo0,cc0,ss0,dd0,tth0);

        nfil0 = nfinal(jo0,ss0,nrec);
        if(layerrefl == 1) nfil0 = zlaymax;
        #ifdef DEBUG
	    fprintf(stderr,"n34=%d dep=%f nfil0=%d zlay=%d\n", n34, dep, nfil0, zlay);
        #endif

        for(istress=0; istress<=1; istress++){
        /*   nstress =0        w   velocity
             nstress =12       T12 stress
             nstress =1        u   velocity
             nstress =22       T22 stress
        */
            if(n34==4 && istress==0) nstress=12;
            if(n34==4 && istress==1) nstress=0;
            if(n34==3 && istress==0) nstress=1;
            if(n34==3 && istress==1) nstress=22;

            #ifdef DEBUG
	        fprintf(stderr, "Rays for DIRECT wave, nstress=%d\n",nstress);
            #endif

            getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil0,jo,nb,
         	    &lfinal0,nen0,ncoun0,na0,nray0,jo0,&nb0,bounce,DIRECT);

            int ix;
            for(ix=0; ix<nx; ix++){
                xx=xmin+ix*h+0.5*istress*h;
                so = (xx > 0) ? so0 : so1;
                grt_(&isource,&iflat,&ispecial,&mtd,&nb0,&jo0,cc0,ss0,dd0,tth0,
		             &lfinal0,nen0,na0,nray0,ncoun0,
                     &dp,&nt,&xx,&tstart,&nstress,so,&der,green);
                convt_(green,&nt,green,dso,&nconv,&l2n,&nx1,&ncent,&fscl);

                fwrite(green,sizeof(float),nt,fop_green);
            }
        }
    }

    for(n34=3; n34<=4; n34++){
        xx = ( (n34==3) ? xmin : xmin+0.5*h);
        so = (xx > 0) ? so0 : so1;
        #ifdef DEBUG
	    fprintf(stderr,"n34=%d xx= %f \n", n34, xx);
        #endif

	    for(istress=0; istress<=1; istress++){
            /*  nstress =0        w   velocity
                nstress =12       T12 stress
                nstress =1        u   velocity
                nstress =11       T11 stress
            */
            if(istress==1 && n34==4) nstress =0;
            if(istress==0 && n34==3) nstress =1;
            if(istress==0 && n34==4) nstress =11;
            if(istress==1 && n34==3) nstress =12;

            #ifdef DEBUG
	        fprintf(stderr, "Rays for TOTAL/Reflected wave, nstress=%d\n",nstress);
            #endif

	        ispecial =  0;
            int iz;
            for(iz=-nd; iz<=nz; iz++){
                dep =(iz-1.0)*h+istress*0.5*h+0.25*h;

		        if(iz < 0 || (iz == 0 && istress == 0))
		            raytype = REFLECTED;
                else
		            raytype = TOTAL;

                if(xx < 0 && raytype == REFLECTED)
                    raytype = DIRECT;

                if(xx > 0 || raytype == DIRECT){
                    zlay = modify_model(dep,nrec,nfil,jo,cc,ss,dd,tth,&jo0,cc0,ss0,dd0,tth0);
                    nfil0 = nfinal(jo0,ss0,nrec);
                    if(layerrefl == 1)  nfil0 = zlaymax;

                    getrays(lfinal,nen,na,nray,ncoun,zlay,nrec,nfil0,jo,nb,
                        &lfinal0,nen0,ncoun0,na0,nray0,jo0,&nb0,bounce,raytype);

                    grt_(&isource,&iflat,&ispecial,&mtd,&nb0,&jo0,cc0,ss0,dd0,tth0,
		                 &lfinal0,nen0,na0,nray0,ncoun0,
                         &dp,&nt,&xx,&tstart,&nstress,so,&der,green);
                    convt_(green,&nt,green,dso,&nconv,&l2n,&nx1,&ncent,&fscl);

                    if(raytype == DIRECT){
                        for(i = 0; i<nt; i++)
                            green[i] *= -1.0;
                    }
                } else {
                    for(i = 0; i<nt; i++)
                        green[i] = 0.0;
                }
                fwrite(green,sizeof(float),nt,fop_green);
            }
        }
    }

    fclose(fop_green);
    return 0;
}

/* getrays determins the ray parameters for the responses
   at zlay layer. The input is the ray group for the n=3,4,
       which are: lfinal, nen, na, nray
                  zlay: the layer where the receive is
                  nrec: the layer where n=3,4 lines are
                  nfil: the layer where solid-liquid interface is
                  type =0: direct wavefield
                       =1: reflected wavefield
                       =2: total wavefield
*/

void getrays(int lfinal, int *nen, int *na, int *nray, int *ncoun,
            int zlay, int nrec, int nfil, int jo, int nb,
            int *lfinal0, int *nen0, int *ncoun0, int *na0, int *nray0,
            int jo0, int *nb0, int bounce, int raytype)
/*int lfinal, *nen, *na, *nray, *ncoun;
int zlay, nrec, nfil, bounce, raytype, jo, nb, jo0, *nb0;
int *lfinal0, *nen0, *na0, *nray0, *ncoun0;
*/
{
    int ray, k, n, j, rays1, rays2, nenn, rays, incre=0;
    int P_SV, SV_P;

    int **rtype, *type;
    int nbounce, i, m, num;
    int nref, nref1, nref2;

    *nb0=nb;
    if(nfil<nrec+1 && jo != jo0){
        *nb0=nb+1; nrec++; incre++;
    }

    n=k=0;
    for (ray=0; ray<lfinal; ray++){
        P_SV = nray[nen[ray]-2];
	    SV_P = ( (P_SV==3) ? 5 : 3);

	    switch(raytype){

	    case DIRECT:
	        rays1=0; rays2=0; break;
	    case REFLECTED:
	        rays1=1; rays2=1; break;
	    case TOTAL:
	        rays1=0; rays2=1; break;
        }

        /*	if(zlay==nfil) rays2=0;  */

    	rtype = imatrix(0,bounce,1,power(2,bounce+2));
	    type  = ivector(0,bounce);

        for(rays=rays1; rays<=rays2; rays++){
            if(rays == 0){
                nenn=0;
                /* keep the ray paths to n=3, 4 */
                for(j=0; j<nen[ray]-1; j++){
                    na0[k]   =na[j]+incre;
                    nray0[k] =nray[j];
                    k++;
                    nenn++;
                }

                /* calculate the direct wave from the source */
		        nbounce  = 1;
		        type[0] = P_SV;
		        nen0[n] = nenn;

		        k  = MultiBounces(nrec,nfil,zlay,k,nbounce,type,
				   na0,nray0,&ncoun0[n],&nen0[n]);
                n++;
            } else {
                /* reflected waves from the solid-liquid boundary */
                for(nbounce=2; nbounce<=bounce; nbounce++){
		            if((zlay == nfil) && ((nbounce % 2) == 0))
			            continue;     /* skip the uparriving ray */
		            if((zlay == nrec) && ((nbounce % 2) != 0))
			            continue;     /* skip the downarriving ray */

                    num = BounceRays(nbounce,rtype,P_SV,SV_P);

		            for(j=1; j<=num; j++){
                        /* modified for only valid when nfil > nrec */
                        /* not tested for the other cases */
                        nref1 = nfil;
                        nref2 = (layerrefl == 1) ? zlay +1 : nfil + 1;
                        for(nref = nref1; nref > nref2; nref--){
                            nenn=0;

                            /* keep the ray paths to n=3, 4 */
                            for(m=0; m<nen[ray]-1; m++){
                                na0[k]   =na[m]+incre;
                                nray0[k] =nray[m];
                                k++;
                                nenn++;
                            }
			                for(i=0; i<nbounce; i++){
			                    type[i]  = rtype[i][j];
//			                    printf("%d ",type[i]);
                            }
//			                printf("\n");
		    		        nen0[n] = nenn;
		                    k  = MultiBounces(nrec,nref,zlay,k,nbounce,type,
				                    na0,nray0,&ncoun0[n],&nen0[n]);
                            n++;
                        }
                    }
                }
            }
        }
	    free_imatrix(rtype,0,bounce,1,power(2,bounce+2));
	    free_ivector(type,0,bounce);
    }
    *lfinal0=n;
}

/*
    int BounceRays(int nbounce, int ** rtype, P_SV SV_P

    rtype[i][j]  i -> ray segments j -> number of rays

*/
int BounceRays(int nbounce, int **rtype, int P_SV, int SV_P)
{
    int i, j, m, num;

    /* first fundamental ray */
    for(i=0; i<nbounce; i++)
         rtype[i][1] = P_SV;

    num = 1;

    if(convert == 1){
        for(i=1; i< nbounce; i++){
            for(j=1; j<=num; j++){
	            for(m=0; m<nbounce; m++) rtype[m][j+num] = rtype[m][j];
                rtype[nbounce - i][j+num] = SV_P;
            }
	        num *= 2;
        }
    }
    return num;
}

int MultiBounces(int nrec, int nfil, int zlay, int k, int bounce,int *type,
        int *na, int *nray, int *ncoun, int *nen)
{
    int j, i; int start, end1, end2, temp; int going;
    int nen0, k0;

    start = ((nfil > nrec) ? nrec + 1 : nrec);
    end1  = ((nfil > nrec) ? nrec + 1 : nrec);
    end2  = ((nfil > nrec) ? nfil - 1 : nfil +1);
    going = ((nfil > nrec) ? DOWN     : UP);
    if((zlay == nrec) && (nfil > nrec))
	end1 = nrec + 2;

    nen0 = *nen;  k0 = k;
    i = 0;
    while ( i < bounce){
	    if(i == bounce-1)	        /* last segement of the ray */
	        end2  = ((going == UP) ? zlay+1 : zlay);

        if(going == DOWN){
            for(j=start; j<=end2; j++){
                na[k0]   = j;
                nray[k0] = type[i];
                if((nfil > nrec) && (j >= nfil))
                    nray[k0] = 5;    /* fixed to P wave */
                if((nfil < nrec) && (j <= nfil))
                    nray[k0] = 5;    /* fixed to P wave */
                k0++; nen0++;
            }
            going = UP;
        } else {
	        for(j=start; j>=end2; j--){
	            na[k0]   = j;
	            nray[k0] = type[i];
		        if((nfil > nrec) && (j >= nfil))
		            nray[k0] = 5;    /* fixed to P wave */
		        if((nfil < nrec) && (j <= nfil))
		            nray[k0] = 5;    /* fixed to P wave */
	            k0++; nen0++;
            }
	        going = DOWN;
        }
	    if(i==0){
	        start = end2; end2 = end1;
        } else {
	        temp = start; start = end2; end2 = temp;
	    }
	    i++;
    }
    na[k0]   = 1;
    nray[k0] = 5;
    k0++; nen0++;
    *ncoun  = -going;
    *nen    =  nen0;
    return k0;
}

/* modify the model parameters for each point at m=0*/
/* returns the layer number of the receiver */

int modify_model(float dep, int nrec, int nfil, int jo,
        float *cc, float *ss, float *dd, float *tth,
        int *jo0, float *cc0, float *ss0, float *dd0, float *tth0)
{
    int layer, zlay;
    float thickness=0.0;

    if(nrec<nfil){
        for(layer=nrec; layer<nfil+500; layer++){
	        thickness += tth[layer];
	        if(thickness >= dep) break;
        }
        zlay = ((dep<0.00001) ? nrec : layer+1);
    } else {
        for(layer=nrec-2; layer>nfil-1; layer--){
	        thickness += tth[layer];
	        if(thickness >= dep) break;
        }
        zlay = ( (dep<0.00001) ? nrec-1 : layer+1);
    }

    if(zlay == nrec) thickness = 0.0;

    if(dep > thickness)zlay = nfil;

    if(fabs(thickness-dep) <0.0001 || zlay==nfil){
        for(layer=0; layer<jo; layer++){
            cc0[layer]  =  cc[layer];
            ss0[layer]  =  ss[layer];
            dd0[layer]  =  dd[layer];
            tth0[layer] = tth[layer];
        }
	    if(zlay==nfil)
            tth0[nfil-1] = dep-thickness;
        *jo0 =jo;
    } else {
        for(layer=0; layer<zlay; layer++){
            cc0[layer]  =  cc[layer];
            ss0[layer]  =  ss[layer];
            dd0[layer]  =  dd[layer];
            tth0[layer] = tth[layer];
        }

        cc0[zlay]  =  cc[zlay-1]+0.00001;
        ss0[zlay]  =  ss[zlay-1]+0.00001;
        dd0[zlay]  =  dd[zlay-1]+0.00001;
        tth0[zlay]   = (nrec<nfil) ? thickness-dep : tth[zlay-1]-(thickness-dep);
        tth0[zlay-1] = (nrec>nfil) ? thickness-dep : tth[zlay-1]-(thickness-dep);

        for(layer=zlay+1; layer<jo+1; layer++){
            cc0[layer]  =  cc[layer-1];
            ss0[layer]  =  ss[layer-1];
            dd0[layer]  =  dd[layer-1];
            tth0[layer] = tth[layer-1];
        }

        *jo0 =jo+1;
    }

    return zlay;
}

int power(int a, int n)
{
    int j, num;

    num = a;
    for(j=2; j<=n; j++)
	    num *= a;

    return num;
}
