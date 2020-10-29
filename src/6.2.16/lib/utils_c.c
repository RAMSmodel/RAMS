#include "utils_sub_names.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <glob.h>
#include <errno.h>

/*****************************************************************************/
 void c_listfile (char*pattern,char*fstring,int patlen,int fslen) {
   glob_t pglob;
   int i, nc, nb, err;
   char colon[]=":";  

/* Get a directory listing of the files with "pattern".
   Concat all names into a colon-delimted string.      */

   glob(pattern, GLOB_ERR, NULL, &pglob);
   
   nb=0;
   for (i = 0; i < pglob.gl_pathc; i++) {
      nc = strlen(pglob.gl_pathv[i]);
      
      if (nb+nc+1 > fslen) {
         printf("\n c_listfile: Fatal error: String length not long enough\n\n");
         exit(1);
      }
      
      /*printf("found: %i %d %d %s\n",fslen,nb,nc,pglob.gl_pathv[i]);*/

      memcpy (&fstring[nb],pglob.gl_pathv[i],nc);
      nb = nb + nc;
      memcpy (&fstring[nb],colon,1);
      nb = nb + 1;
      
   }
   globfree(&pglob);

}

/*****************************************************************************/
 void irsleep (int*seconds)
{
   extern int sleep(int);

/*MAY NEED MACHINE DEPENDENT ARGUMENTS OR USE STATEMENTS FOR THIS ROUTINE*/

#if defined (PC_LINUX1)
   sleep (*seconds);
#else
   print*,"You specified machine/OS other than PC_LINUX1"
   print*,"You need to modify utils_c.c to add your machine/OS"
   stop
#endif

   return;
}
