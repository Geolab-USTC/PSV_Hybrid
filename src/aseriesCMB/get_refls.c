/*************************************************************************/
/*                                                                       */
/*                              get_refls                                */
/*         get the ray descriptions for a certain specification          */
/*                                                                       */
/*************************************************************************/

/* This subroutine returns the ray descriptions of a certain type of waves
   It generates the reflections from source layer to the receiver layer.
   The Reflections include the first reflection on the layers from top
   to the bottom.

   Modify History:
       7 Mar 1997    Lianxing Wen     Initial Revision

*/

#include <stdio.h>
#include <stdlib.h>
#include "rayconst.h"

int get_refls(int nb,         /* source layer */
              int top,        /* toppest layer of reflection */
              int bottom,     /* bottomest layer of reflection */
              int ncmb,       /* the bottomest layer of the mantle */
              int nrec,       /* receiver layer */
              int *lfinal,    /* number of rays */
              int **nen,      /* number of segments of rays*/
              int **ncoun,    /* type of going */
              int **na,       /* passing layers of the ray */
              int **nray,     /* type of passing ray */
              int raytype     /* raytype */
)
{
    int n, k, layer, j, nenn;
    int rtypedown, rtypeup;

    int *nen0, *ncoun0, *na0, *nray0;

    *nen    = (int *) malloc (RAY_NUMBER * sizeof(int));
    *ncoun  = (int *) malloc (RAY_NUMBER * sizeof(int));
    *na     = (int *) malloc (RAY_LENGTH * sizeof(int));
    *nray   = (int *) malloc (RAY_LENGTH * sizeof(int));

    nen0 = *nen; ncoun0 = *ncoun; na0 = *na; nray0 = *nray;

    if(nen0 == NULL || ncoun0 == NULL || na0 == NULL || nray0 == NULL){
	    fprintf(stderr, "Can not allocate memory\n");
	    exit(-1);
    }

    rtypedown = rtypeup = raytype;

    if(raytype == SH && bottom > ncmb){
        fprintf(stderr, "Warning SH wave can not propagate in the fluid\n");
        fprintf(stderr, "On contributions in the mantle are calculated\n");
        bottom = ncmb;
    }
    if(raytype == PS){
	    rtypedown = P;
	    rtypeup   = SV;
    }
    if(raytype == SP){
	    rtypedown = SV;
	    rtypeup   = P;
    }

    n = k = 0;
    for(layer = top; layer <= bottom; layer++){
        nenn = 0;
	    /* for downgoing */
        for(j = nb+1; j<=layer; j++){
            na0[k]    = j;
            nray0[k]  = ((j > ncmb) ? P : rtypedown);
            k++;
            nenn ++;
        }
	    if(layer >= nrec){
            for(j = layer; j >= nrec; j--){
                na0[k]    = j;
                nray0[k]  = ((j > ncmb) ? P : rtypeup);
                k++;
                nenn ++;
            }
        } else {
            for(j = layer+1; j <= nrec; j++){
                na0[k]    = j;
                nray0[k]  = ((j > ncmb) ? P : rtypeup);
                k++;
                nenn ++;
            }
        }
        na0[k]   = 1;
        nray0[k] = 3;
        k++; nenn++;

        ncoun0[n]  = (layer >= nrec) ? UP : DOWN;      /* ray if uparriving */
        nen0[n]    = nenn;
        n ++;

	    if(n >= RAY_NUMBER){
	        fprintf(stderr, "Warning RAY_NUMBER is too small\n");
	        break;
        }
	    if(k + 2* layer -nb - nrec + 2 >= RAY_LENGTH){
	        fprintf(stderr, "Warning RAY_LENGTH is too small\n");
	        break;
        }
    }
    *lfinal = n;

    return 1;
}
