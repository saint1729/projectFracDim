#include <stdio.h>
#include <string.h>

int main(int argc, char* argv[]) {

    /*char date[50];
    int i = 0;
    char slash[1] = "/";
    date[0] = argv[1];*/
    
    char str1[]= "To be or not to be";
    char str2[40];
    char str3[40];

    char str4[50];
    
    /* copy to sized buffer (overflow safe): */
    strncpy ( str2, argv[1], sizeof(str2) );

    /* partial copy (only 5 chars): */
    strncpy ( str3, str2, 5 );
    str3[5] = '\0';   /* null character manually added */
    strcat (str3, "_");

    puts (str1);
    puts (str2);
    puts (str3);
    puts ("_");
    //puts (str4);
    
    
    //printf("%s", date);
    /*while (argv[1][i] != slash[1]) {
        date[i] = argv[1][i];
        printf ("%c\n", date[i]);
        i ++;
    }/**/
    

return 0;
}