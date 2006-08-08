#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "mrisurf.h"
#include "utils.h"
#include "error.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "diag.h"
#include "volume_io.h"
#include "region.h"
#include "machine.h"
#include "fio.h"
#include "mri_identify.h" 
#include "mrisurf.h" 
#include "fmriutils.h" 
#include "gca.h"
#include "gcsa.h"
#include "fsgdf.h"
#include "icosahedron.h"
#include "gca.h"
#include "gcamorph.h"
#include "DICOMRead.h"
#include "fsenv.h"
#include "fsgdf.h"

/* This should be in ctype.h, but the compiler complains */
#ifndef Darwin
#ifndef isblank
int isblank (int c);
#endif
#endif

int striplessthan(char *item);
char *gdfGetNthItem(char *line, int nth);
int LoadDavidsTable(char *fname, int **pplutindex, double **pplog10p);

char *Progname = "dngtester";

/*----------------------------------------*/
int main(int argc, char **argv)
{
  int *lutindex;
  double *log10p;
  LoadDavidsTable(argv[1], &lutindex, &log10p);
  return(0);
  exit(0);
}

/*---------------------------------------------------------*/
char *gdfGetNthItem(char *line, int nth)
{
  char *item;
  int nitems,n;
  static char fmt[2000], tmpstr[2000];
  
  memset(fmt,'\0',2000);
  memset(tmpstr,'\0',2000);

  nitems = gdfCountItemsInString(line);
  if(nth < 0) nth = nitems-1;
  if(nth >= nitems){
    printf("ERROR: asking for item %d, only %d items in string\n",nth,nitems);
    printf("%s\n",line);
    return(NULL);
  }

  for(n=0; n < nth; n++) sprintf(fmt,"%s %%*s",fmt);
  sprintf(fmt,"%s %%s",fmt);
  //printf("fmt %s\n",fmt);
  sscanf(line,fmt,tmpstr);

  item = strcpyalloc(tmpstr);
  return(item);
}
/*--------------------------------------------------*/
int striplessthan(char *item)
{
  int n;

  for(n=0; n < strlen(item); n++)
    if(item[n] == '<') item[n] = '0';
  return(0);
}


/*--------------------------------------------------
  gdfCountItemsInString() returns the number of items
  in the given string, where an item is defined as
  one or more contiguous non-blank characters.
  --------------------------------------------------*/
int gdfCountItemsInString(char *str)
{
  int len, n, nhits;

  len = strlen(str);

  nhits = 0;
  n = 0;
  while(n < len){
    while(isblank(str[n])) n++;
    if(n >= len) break;
    if(str[n] == '\0' || str[n] == '\n' || str[n] == '\r') break;
    while(!isblank(str[n])) n++;
    nhits++;
  }

  //printf("nhits %d\n",nhits);

  return(nhits);
}

/*---------------------------------------------------------------------
  1. Looks for a line where the first string is "Unpaired". 
  2. Gets the last string from this line.
  3. Removes the last 4 chars from this string to get the
     segmentation name (this could be a prob?)
  4. Finds the segmentation name in the color table to get the index
     This is returned in pplutindex.
  5. Scrolls thru the file until it finds a line where the first
     string is "Mean".
  6. Goes one more line
  7. Loads the last value from this line (this is the p-value)
  8. Removes any less-than signs (ie, "<")
  9. Computes -log10 of this value (returns in pplog10p)

  Example:

   Unpaired t-test for Left-Inf-Lat-Vent_vol
   Grouping Variable: Dx
   Hypothesized Difference = 0
   Inclusion criteria: AGE > 60 from wmparc_vals_hypotest_log (imported)
           Mean Diff.      DF      t-Value P-Value
   Dementia, Nondemented   505.029 164     4.136   <.0001

  The segmentation name would be: Left-Inf-Lat-Vent, the p value
  would be .0001 (the -log10 of which would be 4.0).

  ---------------------------------------------------------------------*/

int LoadDavidsTable(char *fname, int **pplutindex, double **pplog10p)
{
  FSENV *fsenv;
  int err,tmpindex[1000];
  double tmpp[1000];
  char tmpstr[2000],tmpstr2[2000],segname[2000];
  FILE *fp;
  double p;
  int n,segindex,nitems;
  char *item;

  fsenv = FSENVgetenv();
  fp = fopen(fname,"r");
  if(fp == NULL) {
    printf("ERROR: could not open%s\n",fname);
    exit(1);
  }

  nitems = 0;
  while(1){
    if(fgets(tmpstr,2000-1,fp) == NULL) break;
    memset(tmpstr2,'\0',2000);
    sscanf(tmpstr,"%s",tmpstr2);
    if(strcmp(tmpstr2,"Unpaired") == 0){
      item = gdfGetNthItem(tmpstr,-1); // get last item
      sscanf(item,"%s",segname);
      free(item);
      // strip off _vol
      for(n=strlen(segname)-4;n<strlen(segname);n++) segname[n]='\0';
      err = CTABfindName(fsenv->ctab, segname, &segindex);
      if(segindex < 0){
	printf("ERROR: reading %s, cannot find %s in color table\n",
	       fname,segname);
	printf("%s",tmpstr);
	printf("item = %s\n",item);
	exit(1);
      }
      n=0;
      while(1){
	fgets(tmpstr,2000-1,fp);
	sscanf(tmpstr,"%s",tmpstr2);
	if(strcmp(tmpstr2,"Mean") == 0) break;
	n++;
	if(n > 1000){
	  printf("There seems to be an error finding key string 'Mean'\n");
	  exit(1);
	}
      }
      fgets(tmpstr,2000-1,fp);
      item = gdfGetNthItem(tmpstr,-1);
      striplessthan(item);
      sscanf(item,"%lf",&p);
      printf("%2d %2d %s %lf  %lf\n",nitems,segindex,segname,p,-log10(p));
      tmpindex[nitems] = segindex;
      tmpp[nitems] = p;
      nitems++;
      free(item);
    }
  }

  *pplutindex =    (int *) calloc(nitems,sizeof(int));
  *pplog10p   = (double *) calloc(nitems,sizeof(double));

  for(n=0; n < nitems; n++){
    (*pplutindex)[n] = tmpindex[n];
    (*pplog10p)[n]   = -log10(tmpp[n]);
  }


  FSENVfree(&fsenv);
  return(nitems);
}

