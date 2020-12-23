#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "getpar.h"

int line2point_();
char iname[20];
int n;
float amom,dt1,dt2,dt3,rhos,dt;

int main(int ac, char **av)
{
    int	ite  = 0;
    int	key  = 1;
    int	ndis2 = 1;
    int idev  = 0;
    float po    = 1.0;
    float dist  = 1.0;

	setpar(ac,av);

    mstpar("accfile",       "s", iname);
    mstpar("Moment",        "f", &amom);

    mstpar("Conv_Source",   "d", &ite);
    mstpar("Dt1",           "f", &dt1);
    mstpar("Dt2",           "f", &dt2);
    mstpar("Dt3",           "f", &dt3);
    mstpar("beta",          "f", &rhos);
    getpar("point_source",  "d", &key);
    mstpar("nt_kir",        "d", &n);
    mstpar("dt_WKM",        "f", &dt);
    getpar("po",            "f", &po);
    getpar("dist",          "f", &dist);

    getpar("nrecv",         "d", &ndis2);
    getpar("finaldev",      "d", &idev);
    endpar();

    line2point_(iname,&amom,&ite,&dt1,&dt2,&dt3,&rhos,
                 &key,&n,&dt,&po,&dist,&ndis2,&idev);
    return 0;
}
