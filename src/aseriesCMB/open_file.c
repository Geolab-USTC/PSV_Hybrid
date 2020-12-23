#include <stdio.h>
#include <stdlib.h>

FILE *open_file(char *file_name, char *access_mode)
{
    FILE *fop;

    if((fop=fopen(file_name, access_mode)) == NULL) {
        fprintf(stderr, "Error Opening file %s with access %s\n",
                        file_name, access_mode);
        exit(-1);
    }
    return fop;
}
