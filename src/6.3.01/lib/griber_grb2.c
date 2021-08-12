#include "sub_gribnames.h"

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

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

static char projection[25], funit[4];
static int  bigdate, rec, pos;
static int center,subcenter,process,table;
static float la1=amiss, lo1=amiss, orient=amiss, lov=amiss, latin1=amiss
            ,latin2=amiss, latind=amiss, dx=amiss, dy=amiss;
static float la2=amiss, lo2=amiss, frecs[MAXLINES];
static int rnum=miss,rpos=miss,kpds5=miss,kpds6=miss,kpds7=miss
          ,grid=miss,nx=miss,ny=miss,mode=miss,scan=miss,nrec=-1
          ,levs[MAXLINES],ifhrs[MAXLINES]
          ,idates[MAXLINES];
static char var[10],lev[20],typ[20],req[20],fields[MAXLINES][20];

static char *progname="dummy";

/**********************************************************************/
/**** Saleeby (3-16-06) Add Variables here that you need to extract ***/
/**********************************************************************/
static char *vars[]= {":GP:",":UGRD:",":VGRD:",":TMP:",":HGT:",":RH:"
                     ,":Z:",":R:",":T:",":U:",":V:",":Q:"
                     ,":WEASD:",":SNOD:",":SOILW:",":TSOIL:",":VSOILM:"
                     ,":SPFH:",":PRMSL:",":PRES:",":WTMP:",NULL};

 void grib_queryproj ();
 void grib_queryfield ();

/**********************************************************************/
 void grib_query (char *filein, int fileinlen) {

/* Initialize the integer arrays returned to fortran */
for (i=1; i<=MAXLINES; i++) {
   frecs[i] =amiss;
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
} else if (!strncmp(req,"frecs",reqlen)) {
   for (i=0; i<=nrec; i++) {
      *(ret+i)=frecs[i];
   }
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
sprintf(cmd,"./wgrib2 -V -end %s",gribfile);
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

   /* Parse the projection information */
   sscanf(lines[0],"%d:%d:vt=%d:%*[^:]:%*[^:]:%s"
                   ,&rec,&pos,&bigdate,var);
   sscanf(lines[3]," %[^:]",projection);

   if (!strncmp(projection,"polar",5)) {
      sscanf(lines[4],"%*[^:]: Lat1 %f Long1 %f Orient %f",&la1,&lo1,&orient);
      sscanf(lines[5]," %*s %*s (%*d x %*d) Dx %f Dy %f scan %d",&dx,&dy,&scan);
      /* Assume the mode is 8 */
      mode=8;
   }

   if (!strncmp(projection,"Lambert",7)) {
      sscanf(lines[4],"%*s %f %*s %f %*s %f",&la1,&lo1,&lov);
      sscanf(lines[5],"%*s %f %*s %f %*s %f",&latind,&latin1,&latin2);
      sscanf(lines[7],"%*s %*s (%d x %d) Dx %f m Dy %f m mode %d",&nx,&ny,&dx,&dy,&mode);
   }

   if (!strncmp(projection,"lat-lon",7)) {
      sscanf(lines[3]," %*[^:]:(%d x %d)",&nx,&ny);
      sscanf(lines[4]," lat %f to %f by %f",&la1,&la2,&dx);
      sscanf(lines[5]," lon %f to %f by %f",&lo1,&lo2,&dy);
      /*printf("steve %d %d %d %s\n",la1,pos,bigdate,var);*/
   }

}
pclose(pipe);

}

/**********************************************************************/
 void grib_queryfield () {

int k;
char *var;
char fhr[20],sslev[20];
int num,rec,plev,slev;
float flev;

sprintf(cmd,"./wgrib2 -s %s",gribfile);

if ((pipe = popen(cmd, "r")) != NULL) {

   /* Read the wgrib short description of the first record */
   if (fread(buf, 1, BUFSIZE, pipe) >= BUFSIZE) {
      fprintf(stderr,"Maximum BUFSIZE exceeded\n");
      fprintf(stderr,"Increase BUFSIZE in grib_query\n");
      /*exit(1);*/
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

   /* Loop through each grib inventory record */
   for (i=0; i<nlines; i++) {

/******************************************************************/
/*** Saleeby (3-16-06) Add Units string here that is found in the
     wgrib inventory. "mb" is for pressure levels, "sfc" is for 
     surface fields, "MSL" is for mean sea level pressure fields
     and "cm down" is for soil level fields.
*******************************************************************/
      if (strstr( lines[i]," mb:")  ||
          strstr( lines[i],":surface:") ||
          strstr( lines[i],":mean sea level:") ||
          strstr( lines[i]," below ground:")) {

         k=0;
         while (var=vars[k++]) {
            if (strstr(lines[i],var)) {
               nrec++;
               sscanf(lines[i],"%f:%d:d=%d:%[^:]:%[^:]:%[^:]"
                     ,&frecs[nrec],&rec,&idates[nrec],typ,lev,fhr);

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
               if (strstr(lev,"surface")) {
                  levs[nrec]=0;
               } else if (strstr(lev,"mean sea level")) {
                  levs[nrec]=0;
               } else if (strstr(lines[i],"below ground")) {
                  sscanf(lev,"%*[^-]-%s",sslev);
                  flev=atof(sslev);
                  slev=(int)(flev*100.01);
                  levs[nrec]=-1*slev;
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
 void grib_get (float *data, float *readrec) {

FILE *input;

int res;
int i;
char dummy[250];

/* Read the wgrib short description of the first record */
if(*readrec==floor(*readrec)){
 sprintf(cmd,"./wgrib2 %s -d %5.0f -inv /dev/null -text ./file.txt",gribfile,*readrec);
 /*printf(" Record %5.0f: %i  %i\n",*readrec,nx,ny);*/
}
else if(*readrec!=floor(*readrec)){
 sprintf(cmd,"./wgrib2 %s -d %5.1f -inv /dev/null -text ./file.txt",gribfile,*readrec);
 /*printf(" Record %5.1f: %i  %i\n",*readrec,nx,ny);*/
}
system(cmd);

if((input = fopen("file.txt","r")) != NULL) {
   i=0;
   fscanf(input,"%*[^\n]\n");
   while(!feof(input)){
    fscanf(input,"%[^\n]\n",&dummy);
    sscanf(dummy,"%f",&data[i]);
    i++;
   }
   close(input);
}

}
