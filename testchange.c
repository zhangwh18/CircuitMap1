#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

/*参数： filepath - 要修改的文件   str - 修改为的字符串 */
void ModifyFile(char *filepath, char *str)
{

    char linebuffer[512] = {0};
    char buffer1[512] = {0};
    char buffer2[512] = {'a'};
    // char buf[512] = {0};
    char cmd[256] = "rm temp.txt";
    printf("1\n");
    FILE *fp = fopen(filepath, "r");
    FILE *fpttmp = fopen("temp.txt", "wt+");
    printf("1\n");

    while (fgets(linebuffer, 512, fp))
    {
        printf("%c\n", buffer2[0]);

        sscanf(linebuffer, "%[^_]_%[^_]", buffer1, buffer2);
        if (!strcmp("11", buffer1))
        {
            printf("1\n");
            memset(linebuffer, '\0', strlen(linebuffer) + 1);
            memcpy(linebuffer, str, strlen(str) + 1);
            // linebuffer[7]  = '\0';
            fprintf(fpttmp, "%s", linebuffer);
            fprintf(fpttmp, "\n");
        }
        else
        {
            fprintf(fpttmp, "%s", linebuffer);
        }
    }
    fclose(fp);
    fclose(fpttmp);

    //清空文件中的内容
    fp = fopen(filepath, "w");
    fclose(fp);

    fp = fopen(filepath, "wt+");
    fpttmp = fopen("temp.txt", "r");

    while (fgets(linebuffer, 512, fpttmp))
    {

        fprintf(fp, "%s", linebuffer);
    }
    fclose(fp);
    fclose(fpttmp);
    system(cmd);
}

int main()
{
    char name[128] = "wugaoquan";
    ModifyFile("mm.txt", name);
    printf("finished\n");

    return 0;
}