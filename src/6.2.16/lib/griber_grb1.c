#include "sub_gribnames.h"

#include <stdio.h>
#include <math.h>
#include <string.h>

/**********************************************************************/
#define BUFSIZE 2000000
#define MAXLINES 40000
#define CMDLEN 512
#define miss -999
#define amiss -999.

static char gribfile[CMDLEN];

static int i,nlines;
static char cmd[CMDLEN], buf[BUFSIZE], *lines[MAXLINES];
static FILE *pipe;

static char projection[15], funit[4];
static int  bigdate, rec, pos;
static int center,subcenter,process,table;
static float la1=amiss, lo1=amiss, orient=amiss, lov=amiss, latin1=amiss
            ,latin2=amiss, dx=amiss, dy=amiss;
static float la2=amiss, lo2=amiss;
static int rnum=miss,rpos=miss,kpds5=miss,kpds6=miss,kpds7=miss
          ,grid=miss,nx=miss,ny=miss,mode=miss,scan=miss,nrec=-1
          ,irecs[MAXLINES],levs[MAXLINES],ifhrs[MAXLINES]
          ,idates[MAXLINES];
static char var[10],lev[20],typ[20],req[20],fields[MAXLINES][20];

static char *progname="dummy";

/**********************************************************************/
/**** Saleeby (3-16-06) Add Variables here that you need to extract ***/
/**********************************************************************/
static char *vars[]= {":GP:",":UGRD:",":VGRD:",":TMP:",":HGT:",":RH:"
                     ,":Z:",":U:",":V:",":T:",":R:",":Q:"
                     ,":WEASD:",":SNOD:",":SOILW:",":TSOIL:",":VSOILM:"
                     ,":SPFH:",":PRMSL:",":PRES:",":WTMP:",":SWVL1:"
                     ,":SWVL2:",":STL1:",":STL2:",":SD:",NULL};

 void grib_queryproj ();
 void grib_queryfield ();

/**********************************************************************/
 void grib_query (char *filein, int fileinlen) {

/* Initialize the integer arrays returned to fortran */
for (i=1; i<=MAXLINES; i++) {
   irecs[i] =miss;
   levs[i]  =miss;
   ifhrs[i] =miss;
   idates[i]=miss;
}

/* Add a null character to force string termination */
filein[fileinlen]='\0';

strncpy(gribfile,filein,fileinlen);

grib_queryproj ();

grib_queryfield ();

}

/**********************************************************************/
 void grib_queryc (char *req, char *ret, int reqlen, int retlen) {

if (!strncmp(req,"projection",reqlen)) {
   memset(ret,' ',retlen);
   strncpy(ret,projection,strlen(projection));
} else if (!strncmp(req,"funit",reqlen)) {
   memset(ret,' ',retlen);
   strncpy(ret,funit,strlen(funit));
} else if (!strncmp(req,"fields",reqlen)) {
   for (i=0; i<nrec; i++) {
      memset(ret+20*i,' ',retlen);
      strncpy(ret+20*i,fields[i],strlen(fields[i]));
   }
} else {
   fprintf(stderr,"variable %s not available\n",req);
}

}

/**********************************************************************/
 void grib_queryf (char *req, float *ret, int reqlen) {

if (!strncmp(req,"lat1",reqlen)) {
   *ret=la1;
} else if (!strncmp(req,"lon1",reqlen)) {
   *ret=lo1;
} else if (!strncmp(req,"lat2",reqlen)) {
   *ret=la2;
} else if (!strncmp(req,"lon2",reqlen)) {
   *ret=lo2;
} else if (!strncmp(req,"lov",reqlen)) {
   *ret=lov;
} else if (!strncmp(req,"orient",reqlen)) {
   *ret=orient;
} else if (!strncmp(req,"latin1",reqlen)) {
   *ret=latin1;
} else if (!strncmp(req,"latin2",reqlen)) {
   *ret=latin2;
} else if (!strncmp(req,"dx",reqlen)) {
   *ret=dx;
} else if (!strncmp(req,"dy",reqlen)) {
   *ret=dy;
} else {
   fprintf(stderr,"variable %s not available\n",req);
}

}

/**********************************************************************/
 void grib_queryi (char *req, int *ret, int reqlen) {

if (!strncmp(req,"center",reqlen)) {
   *ret=center;
} else if (!strncmp(req,"grid",reqlen)) {
   *ret=grid;
} else if (!strncmp(req,"nx",reqlen)) {
   *ret=nx;
} else if (!strncmp(req,"ny",reqlen)) {
   *ret=ny;
} else if (!strncmp(req,"nrec",reqlen)) {
   *ret=++nrec;
} else if (!strncmp(req,"longdate",reqlen)) {
   *ret=bigdate;
} else if (!strncmp(req,"mode",reqlen)) {
   *ret=mode;
} else if (!strncmp(req,"iscan",reqlen)) {
   *ret=scan;
} else if (!strncmp(req,"irecs",reqlen)) {
   for (i=0; i<=nrec; i++) {
      *(ret+i)=irecs[i];
   }
} else if (!strncmp(req,"idates",reqlen)) {
   for (i=0; i<=nrec; i++) {
      *(ret+i)=idates[i];
   }
} else if (!strncmp(req,"ifhrs",reqlen)) {
   for (i=0; i<=nrec; i++) {
      *(ret+i)=ifhrs[i];
   }
} else if (!strncmp(req,"levs",reqlen)) {
   for (i=0; i<=nrec; i++) {
      *(ret+i)=levs[i];
   }
} else {
   fprintf(stderr,"variable %s not available\n",req);
}

}

/**********************************************************************/
 void grib_queryproj () {
  
/* Verbose query for projection for 1st record*/
sprintf(cmd,"./rams_wgrib -V %s",gribfile);
if ((pipe = popen(cmd, "r")) != NULL) {

   /*printf("executed command=%s\n",cmd);*/
   /* Read the wgrib Verbose description of the first record */
   /*if (fread(buf, 1, BUFSIZE, pipe) >= BUFSIZE) {*/
   /* Read the wgrib Verbose description of the first 1000 bytes which
      is long enough to contain the first record's info */
   if (fread(buf, 1, 1000, pipe) >= BUFSIZE) {
      fprintf(stderr,"Maximum BUFSIZE exceeded\n");
      fprintf(stderr,"Increase BUFSIZE in grib_query\n");
      /*exit(1); Saleeby, Problem with exits in compile*/
   }

   /* split into "lines" */
   i=0;
   nlines=1;
   lines[i++]=strtok(buf,"\n");
   while ( (lines[i++] = strtok(NULL,"\n")) != NULL) {
      nlines++;
      if (nlines > MAXLINES) {
         fprintf(stderr,"Maximum MAXLINES exceeded\n");
         fprintf(stderr,"Increase MAXLINES in grib_query\n");
         /*exit(1);*/
      }
   }
   printf("%i: %s\n",0,lines[0]);
   printf("%i: %s\n",1,lines[1]);
   printf("%i: %s\n",2,lines[2]);
   printf("%i: %s\n",3,lines[3]);
   printf("%i: %s\n",4,lines[4]);
   printf("%i: %s\n",5,lines[5]);
   printf("%i: %s\n",6,lines[6]);
   printf("%i: %s\n",7,lines[7]);
   printf("%i: %s\n",8,lines[8]);
   printf("%i: %s\n",9,lines[9]);
   printf("%i: %s\n",10,lines[10]);

   /* Parse the projection information */
   sscanf(lines[0],"rec %d:%d:date %d %s kpds5=%d kpds6=%d kpds7=%d %*[^ ] grid=%d"
                   ,&rec,&pos,&bigdate,var,&kpds5,&kpds6,&kpds7,&grid);
   sscanf(lines[2]," timerange %*d P1 %*d P2 %*d TimeU %*d nx %d ny %d"
                   ,&nx,&ny);
   sscanf(lines[3]," center %d subcenter %d process %d Table %d"
                   ,&center,&subcenter,&process,&table);
   sscanf(lines[4]," %[^:]",projection);

   /***
   printf("   Date:%d  Projection:%s\n",bigdate,projection);
   printf("   kpds5:%d\n",kpds5);
   printf("   kpds6:%d\n",kpds6);
   printf("   kpds7:%d\n",kpds7);
   printf("   grid:%d\n",grid);
   printf("   nx:%d\n",nx);
   printf("   ny:%d\n",ny);
   printf("   center:%d\n",center);
   printf("   subcenter:%d\n",subcenter);
   printf("   process:%d\n",process);
   printf("   table:%d\n",table);
   ***/

   if (!strncmp(projection,"polar",5)) {
      sscanf(lines[4],"%*[^:]: Lat1 %f Long1 %f Orient %f",&la1,&lo1,&orient);
      sscanf(lines[5]," %*s %*s (%*d x %*d) Dx %f Dy %f scan %d",&dx,&dy,&scan);
      /* Assume the mode is 8 */
      mode=8;
   }

   if (!strncmp(projection,"Lambert Conf",12)) {
      sscanf(lines[4],"%*[^:]: Lat1 %f Lon1 %f Lov %f",&la1,&lo1,&lov);
      sscanf(lines[5]," Latin1 %f Latin2 %f",&latin1,&latin2);
      sscanf(lines[6]," %*s %*s (%*d x %*d) Dx %f Dy %f scan %d mode %d",&dx,&dy,&scan,&mode);
   }

   if (!strncmp(projection,"latlon",12)) {
      sscanf(lines[4]," %*[^:]: lat %f to %f by %f",&la1,&la2,&dx);
      sscanf(lines[5]," long %f to %f by %f, (%*d x %*d) scan %d mode %d"
         ,&lo1,&lo2,&dy,&scan,&mode);
   }

}
pclose(pipe);

}

/**********************************************************************/
 void grib_queryfield () {

int k;
char *var;
char fhr[20];
int num,rec,plev;

sprintf(cmd,"./rams_wgrib -s %s",gribfile);

if ((pipe = popen(cmd, "r")) != NULL) {

   /* Read the wgrib short description of the first record */
   if (fread(buf, 1, BUFSIZE, pipe) >= BUFSIZE) {
      fprintf(stderr,"Maximum BUFSIZE exceeded\n");
      fprintf(stderr,"Increase BUFSIZE in grib_query\n");
      /*exit(1);*/
   }
   /*printf("done\n");*/

   /* split into "lines" */
   i=0;
   nlines=1;
   lines[i++]=strtok(buf,"\n");
   while ( (lines[i++] = strtok(NULL,"\n")) != NULL) {
      nlines++;
      if (nlines > MAXLINES) {
         fprintf(stderr,"Maximum MAXLINES exceeded\n");
         fprintf(stderr,"Increase MAXLINES in grib_query\n");
         /*exit(1);*/
      }
   }

   /* Loop through each grib inventory record */
   for (i=0; i<nlines; i++) {
/******************************************************************/
/*** Saleeby (3-16-06) Add Units string here that is found in the
     wgrib inventory. "mb" is for pressure levels, "sfc" is for
     surface fields, "MSL" is for mean sea level pressure fields
     and "cm down" is for soil level fields.
*******************************************************************/
      if (strstr( lines[i]," mb:")  ||
          strstr( lines[i],":sfc:") ||
          strstr( lines[i],":MSL:") ||
          strstr( lines[i]," cm down:")) {


         k=0;
         while (var=vars[k++]) {
            if (strstr(lines[i],var)) {
               nrec++;
               sscanf(lines[i],"%d:%d:d=%d:%[^:]:%[^:]:%[^:]"
                     ,&irecs[nrec],&rec,&idates[nrec],typ,lev,fhr);
               /*printf("typ=%s lev=%s fhr=%s\n",typ,lev,fhr);*/

               /* Disregard the fields like "0-0hr acc" */
               if (strstr(fhr,"acc")) continue;

               /* Determine the units and time of the forecast */
               if (strstr(fhr,"anl")) {
                  strncpy(funit,"hr",2);
                  ifhrs[nrec]=0;
               } else if (strstr(fhr,"hour")) {
                  strncpy(funit,"hour",4);
                  sscanf(fhr,"%dhour fcst",&ifhrs[nrec]);
               } else if (strstr(fhr,"hr")) {
                  strncpy(funit,"hr",2);
                  sscanf(fhr,"%dhr fcst",&ifhrs[nrec]);
               } else if  (strstr(fhr,"min")) {
                  strncpy(funit,"min",3);
                  sscanf(fhr,"%dmin fcst",&ifhrs[nrec]);
               } else if  (strstr(fhr,"sec")) {
                  strncpy(funit,"sec",3);
                  sscanf(fhr,"%dsec fcst",&ifhrs[nrec]);
               } else if  (strstr(fhr,"day")) {
                  strncpy(funit,"day",3);
                  sscanf(fhr,"%dday fcst",&ifhrs[nrec]);
               }

/**********************************************************************/
/**** Saleeby (3-16-06) Add Variables types here if other than these **/
/**********************************************************************/
               if (strstr(lev,"sfc")) {
                  levs[nrec]=0;
               } else if (strstr(lev,"MSL")) {
                  levs[nrec]=0;
               } else if (strstr(lines[i],"cm down")) {
                  sscanf(lev,"%*d%d",&levs[nrec]);
               } else if (strstr(lines[i],"mb")) {
                  sscanf(lev,"%d",&levs[nrec]);
               }
               strncpy(fields[nrec],typ,strlen(typ));
            }
         }
      }
   }
}


pclose(pipe);

}

/**********************************************************************/
 void grib_get (float *data, int *readrec) {

int res;

/* Read the wgrib short description of the first record */
sprintf(cmd,"./rams_wgrib -d %d -nh -o - %s",*readrec,gribfile);

if ((pipe = popen(cmd, "r")) != NULL) {
   /* Read the wgrib results */
   res=fread(data, 4, nx*ny, pipe);
   /*printf(" Record %3d: %5d floats %i  %i\n",rec,res,nx,ny);*/

   pclose(pipe);
}

}
