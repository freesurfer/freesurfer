#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/file.h>
#include <sys/types.h>
#include <time.h>
#include <sys/timeb.h>
#include <sys/time.h>
#include <ctype.h>
#include <dirent.h>
#include <malloc.h>
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri_identify.h"
#include "mri2.h"
#include "imautils.h"
#include "machine.h"

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_probe_ima.c,v 1.2 2002/01/30 22:47:30 greve Exp $";
char *Progname = NULL;

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static int  singledash(char *flag);
static int  stringmatch(char *s1, char *s2);

char* imafile = NULL; 
char* typestring = NULL;
int   type = -1;
int   typesize = 1;
int   offset = -1;
int   stringlen = 1;
char  *key = NULL;
int   keyno = -1;
int   dumpfileinfo = 0;
int debug, verbose;
FILE *fp;
char *attrname;
int  getattr = 0;

#define TMPSTRLEN 10000
static char tmpstr[TMPSTRLEN];

/*---------------------------------------------------------------*/
int main(int argc, char **argv)
{
  void *pVal;
  IMAFILEINFO *ifi;

  tmpstr[0] = 'a'; /* to stop compiler warning */

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  //if(!imaHasIMAExtension(imafile)){
  //  printf("ERROR: file %s does not have extension .ima or .IMA\n",imafile);
  //  exit(1);
  //}

  if(getattr){
    if(stringmatch(attrname,"isima")){
      printf("%d\n",imaIsSiemensIMA(imafile));
      return(0); exit(0);
    }
  }

  if(! imaIsSiemensIMA(imafile) ){
    printf("ERROR: %s is not a Siemens IMA file\n",imafile);
    exit(1);
  }

  MkImaDictionary();

  if(dumpfileinfo){
    ifi = imaLoadFileInfo(imafile);
    imaDumpFileInfo(stdout,ifi);
    return(0); exit(0);
  }

  if(getattr){
    ifi = imaLoadFileInfo(imafile);
    if(stringmatch(attrname,"voldim")){
      printf("%d %d %d\n",ifi->VolDim[0],ifi->VolDim[1],ifi->VolDim[2]);
      return(0); exit(0);
    }
    if(stringmatch(attrname,"volres")){
      printf("%g %g %g\n",ifi->VolRes[0],ifi->VolRes[1],ifi->VolRes[2]);
      return(0); exit(0);
    }
    if(stringmatch(attrname,"nframes")){
      printf("%d\n",ifi->NFrames);
      return(0); exit(0);
    }
    if(stringmatch(attrname,"tr")){
      printf("%g\n",ifi->RepetitionTime);
      return(0); exit(0);
    }
    if(stringmatch(attrname,"pulseseq")){
      printf("%s\n",ifi->PulseSequence);
      return(0); exit(0);
    }
    if(stringmatch(attrname,"patname")){
      printf("%s\n",ifi->PatientName);
      return(0); exit(0);
    }
    if(stringmatch(attrname,"patdob")){
      printf("%s\n",ifi->PatientDOB);
      return(0); exit(0);
    }
    printf("ERROR: attribute %s not recognized\n",attrname);
    return(1); exit(1);
  }

  fp = fopen(imafile,"r");
  if(fp == NULL) {
    printf("ERROR: could not open %s\n",imafile);
    exit(1);
  }

  if(keyno > -1) key = ImaDictionary[keyno].key;
  if(key != NULL){
    pVal = imaLoadValFromKey(fp,key,NULL);
    type = imaTypeFromKey(key);
    imaPrintVal(stdout, type, pVal);
    return(0); exit(0);
  }

  if(offset > -1){
    type = imaTypeFromString(typestring);
    typesize = imaTypeSize[type];
    if(debug)
      printf("type = %s (%d), offset = %d\n",typestring,type,offset);
    pVal = imaLoadVal(fp,offset,typesize,stringlen,NULL);
    imaPrintVal(stdout, type, pVal);
    printf("\n");
    return(0); exit(0);
  }

  DumpImaDictionaryVal(stdout, imafile);

  fclose(fp);
  return(0);

}
/*---------------------------------------------------------------*/
/* -----//////////////////<<<<<<<<>>>>>>>>>>\\\\\\\\\\\\\\\\---- */
/*---------------------------------------------------------------*/

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--verbose")) verbose = 1;
    else if (!strcasecmp(option, "--fileinfo")) dumpfileinfo = 1;

    /* -------- source volume inputs ------ */
    else if (!strcmp(option, "--i")){
      if(nargc < 1) argnerr(option,1);
      imafile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--attr")){
      if(nargc < 1) argnerr(option,1);
      attrname = pargv[0];
      getattr = 1;
      nargsused = 1;
    }
    else if (!strcmp(option, "--dictionary")){
      DumpImaDictionary(stdout);
      exit(0);
      nargsused = 1;
    }
    else if (!strcmp(option, "--key")){
      if(nargc < 1) argnerr(option,1);
      key = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--keyno")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&keyno);
      nargsused = 1;
    }
    else if (!strcmp(option, "--o")){
      if(nargc < 2) argnerr(option,2);
      sscanf(pargv[0],"%d",&offset);
      typestring = pargv[1];
      nargsused = 2;
      if(strcmp(typestring,"string") == 0){
  if(nargc < 3) {
    printf("ERROR: type string needs length argument\n");
    exit(1);
  }
  sscanf(pargv[2],"%d",&stringlen);
  if(stringlen < 1) {
    printf("ERROR: string length = %d, must be >= 1\n",
     stringlen);
    exit(1);
  }
  nargsused = 3;
      }
    }
    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
  fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  fprintf(stdout, "USAGE: %s \n",Progname) ;
  fprintf(stdout, "\n");
  fprintf(stdout, "   --i imafile           : path to ima file \n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --key keystring       : string from dictionary\n");
  fprintf(stdout, "   --o offset type <len> : decimal offset, data type, and len\n");
  fprintf(stdout, "   --attr attrname       : print value of file info attribute\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --fileinfo            : print interpreted file info\n");
  fprintf(stdout, "   --dictionary          : print dictionary (no file needed)\n");
#if 0
  fprintf(stdout, "   --view          : view the image  \n");
  fprintf(stdout, "   --o file        : dump binary pixel data into file\n");
  fprintf(stdout, "   --ob stem       : dump binary pixel data into bshort\n");
#endif
  fprintf(stdout, "   --help                : how to use this program \n");
  fprintf(stdout, "   --version             : print version and exit\n");

  fprintf(stdout, "\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  printf("\n");
  print_usage() ;
  printf("\n");

  printf(
"This program allows the user to query a Siemens IMA file, \n"
"and can be used to print out a single value from the IMA header \n"
"or to dump lots of info.\n");

  printf("\n"
"DESCRIPTION\n"
"\n"
"A single value can be queried in one of three ways: (1) offset/type, \n"
"(2) key string, (3) attribute name. In the offset method, the user supplies \n"
"an offset which is the number of bytes into the header; the user also supplies  \n"
"a type string indicating the data type (short, int, long, float, double, string). \n"
"If a string type is specified, the user must also supply the string length. In  \n"
"the key method, the user supplies a key string from a dictionary (see below). \n"
"In the attribute method, the user supplies an attribute name. The difference \n"
"between an attribute and a key is that an attribute may have been interpreted \n"
"in some way (eg, the raw (key-based) row resolution for mosaics is wrong, \n"
"but the attribute is correct). \n"
"\n"
"A dump of header information can be obtained in two ways. First, if only \n"
"a file name is given, then the dictionary with corresponding values will  \n"
"be dumped. If the --fileinfo flag is added, then values that have been  \n"
"interpreted are printed out. \n"
"\n"
"ARGUMENTS\n"
"\n"
"  --i imafile\n"
"\n"
"      Path to the IMA file to be probe. If this is the only option, \n"
"      the dictionary with corresponding values is printed out. See also\n"
"      --dictionary.\n"
"\n"
"  --o offset type <stringlen>\n"
"\n"
"      offset is the number of bytes from the beginning of the file. type \n"
"      is the type of data. Valid values are short, int, long, float, double, \n"
"      and string. If string is used, the length of the string must be supplied.\n"
"\n"
"  --d key \n"
"\n"
"      key is a string as found in the dictionary (see DICTIONARY below).\n"
"\n"
"  --attr attrname\n"
"\n"
"     Name of an attribute.\n"
"\n"
"  --fileinfo \n"
"\n"
"     Dump the interpreted file information.\n"
"\n"
"  --dictionary\n"
"\n"
"     Dump the dictionary (no ima file name need be supplied). Each entry has\n"
"     six columns: entry number, key string, offset, type string, number of \n"
"     bytes in the type, and the string length. If only the imafile is supplied,\n"
"     the dictionary will be printed out with the value as a seventh column.\n"
"\n"
"DICTIONARY\n"
"\n"
"The dictionary is a list of character strings (keys) that describe an \n"
"aspect of a value in the IMA header along with its offset, data type,  \n"
"and, if a string, the length of the string. The dictionary is unique to \n"
"to this program and does not represent anything official from Siemens \n"
"or anyone else. The key names, offsets, and data types were gleaned from \n"
"other programs and reverse engineering. \n"
"\n"
"ATTRIBUTES\n"
"\n"
"  isima   : returns 1 if the file is a Siemens IMA file, 0 otherwise \n"
"  voldim  : number of columns, rows, and slices in the volume. The number of \n"
"            slices will be wrong if the files does not contain a mosaic \n"
"  volres  : spacing between columns, rows, and slices in the volume\n"
"  nframes : number of frames\n"
"  tr      : repetition time (sec)\n"
"  pulseseq : pulse sequence name\n"
"  patname  : patient name \n"
"  patdob   : patient date of birth (YYYYMMDD)\n"
"\n"
"AUTHOR\n"
"\n"
"Written by Douglas N. Greve.\n"
"\n"
"BUG REPORTING\n"
"\n"
"Send bug reports to analysis-bugs@nmr.mgh.harvard.edu. Make sure to include\n"
"enough information to replicate the problem. This includes, but is not limited\n"
"to, the version and command-line.\n"
"\n"
"BUGS\n"
"\n"
"This program does not make use of any 'official' documentation from Siemens. \n"
"As such, some of the results may be wrong. Some of the header elements may or\n"
"may not be filled with legitimate values, depending upon how the pulse sequence\n"
"was programmed.\n"
"\n"
"The G28_Pre_PixelSize_Row and G28_Pre_PixelSize_Column values are incorrect for\n"
"mosaics. This is actually the way Siemens stores these values and is not an bug in \n"
"this program.\n"
"\n"
"VERSION\n"
"\n"
);
printf("   %s\n\n", vcid) ;



  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if(n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/* --------------------------------------------- */
static void check_options(void)
{

  if(imafile == NULL){
    printf("ERROR: no file name supplied\n");
    exit(1);
  }

  return;

  if(key != NULL && offset > -1){
    printf("ERROR: cannot specify key and offset\n");
    exit(1);
  }

  if(key == NULL && offset < 0){
    printf("ERROR: must specify either key or offset\n");
    exit(1);
  }

  if(offset > -1 && typestring == NULL){
    fprintf(stderr,"ERROR: no data type supplied\n");
    exit(1);
  }

  return;
}
/*-------------------------------------------------------------*/
static int stringmatch(char *s1, char *s2)
{
  if(strcmp(s1,s2) == 0) return(1);
  return(0);
}
