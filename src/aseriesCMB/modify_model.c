/*************************************************************************/
/*                                                                       */
/*                         modify_model.c                                */
/*         Modify the model to fit the GRT calculation                   */
/*                                                                       */
/*************************************************************************/

/* This subroutine modifies the model to fit the Generalized Ray code,
   It modifies the model and finds the source and receiver layer, the
   the bottomest layer in the mantle.

   Modify History:
       7 Mar 1997    Lianxing Wen     Initial Revision

*/
#include<stdio.h>
/* modify the model parameters for each source -receiver pair*/
int modify_model(float sdep,       /* source depth */
                 int   *nb,        /* source layer */
                 int   *jo,        /* number of layers of final layer */
                 float *cc,        /* p-wave velocities */
                 float *ss,        /* s-wave velocities */
                 float *dd,        /* densities */
                 float *th         /* thicknesses */
                )
{
    int layer;

    float thickness = 0.0;
    float tiny = 1.e-6;

    for(layer=0; layer<*jo; layer++){
        thickness += th[layer];
        if(thickness > sdep) break;
    }
    if(sdep < tiny) layer = 0;

    if(layer == 0){
        *nb = 1;
    } else if(sdep == thickness - th[layer]){
       /* no modification needed */
        *nb  = layer;
    } else {
        /* source is between layer -> layer + 1  */
        *nb  = layer + 1;
        for(layer = *jo; layer > *nb; layer--){
            cc[layer]  = cc[layer -1];
            ss[layer]  = ss[layer -1];
            dd[layer]  = dd[layer -1];
            th[layer]  = th[layer -1];
        }
        layer = *nb;
        cc[layer]   =  cc[layer-1] + tiny;
        ss[layer]   =  ss[layer-1] + tiny;
        dd[layer]   =  dd[layer-1] + tiny;
        th[layer]   =  thickness - sdep;

        th[layer-1] = th[layer -1] - th[layer];
        (*jo)++;
    }
    return 1;
}
