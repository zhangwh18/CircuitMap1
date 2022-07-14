#include <stdio.h>

int main(int argc, char *argv[])
{
    FILE *fp;
    fp = fopen("samrat.txt", "w");
    perror("fopen");
    fprintf(fp, "this is a demo text");
    fclose(fp);
}