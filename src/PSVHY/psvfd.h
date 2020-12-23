int NTSIZE = 20000;           /* starting number of table entries */

/* table entries for material parameters, each grid point */
#define    ta(x)      na[pt[x]+0]
#define    tb(x)      na[pt[x]+1]
#define    tc(x)      na[pt[x]+2]
#define    td(x)      na[pt[x]+3]
#define    te(x)      na[pt[x]+4]

struct vel_stress
{
   float u;
   float w;
   float t11;
   float t12;
   float t22;
};
