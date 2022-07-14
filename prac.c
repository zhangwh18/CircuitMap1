// Program to illustrate the getopt()
// function in C

#include <stdio.h>
#include <unistd.h>

int main(int argc, char *argv[])
{
    int opt;

    // put ':' in the starting of the
    // string so that program can
    // distinguish between '?' and ':'
    // printf("2\n");
    opt = getopt(argc, argv, ": if: lrx");
    // while ((opt = getopt(argc, argv, ": if: lrx")) != -1)
    // {
    // printf("1\n");
    // switch (opt)
    // {
    // case 'i':
    // case 'l':
    // case 'r':
    //     printf("option: %d\n", (int)(opt));
    //     break;
    // case 'f':
    //     printf("filename: %s\n", optarg);
    //     break;
    // case ':':
    //     printf("option needs a value\n");
    //     break;
    // case '?':
    //     printf("unknown option: %c\n", (int)(optopt));
    //     break;
    // }
    // }

    // optind is for the extra arguments
    // which are not parsed
    // for (; optind < argc; optind++)
    // {
    //     printf("extra arguments: %s\n", argv[optind]);
    // }
    int z;
    // z = atoi(argv[0]);
    sscanf(argv[1], "%d", &z);
    printf("%d\n", z);
    return 0;
}
