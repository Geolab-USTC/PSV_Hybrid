void findcmbs(int jo, float *ss, int *ncmb, int *icmb)
{
    float *s=ss;

    s++; jo--;   /* skip the first line of model */
    while(jo--)
        if(*s++ < 0.01)break;
    *ncmb = s - ss - 1;

    while(jo--)
        if(*s++ > 0.01)break;
    *icmb = s - ss - 1;
}
