#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <errno.h>

#include "mri.h"
#include "error.h"
#include "diag.h"
#include "fio.h"

#ifndef lint
static char vcid[] = "$Id: mri_convert_mdh.c,v 1.2 2003/05/26 06:36:43 greve Exp $";
#endif /* lint */

#define MDH_SIZE_V15    128        //Number of bytes in the V1.5 miniheader
#define MDH_BM_ONLINE   (1 <<  4)
#define MDH_BM_PHASECOR (1 << 22)
#define MDH_BM_REFLECT  (1 << 25)

typedef struct tagMDH {
  unsigned long ulTimeStamp, BitMask1;
  unsigned short Ncols, UsedChannels, Slice, Partition, Echo, Rep,LoopCounterLine;
  float SlicePosSag,SlicePosCor,SlicePosTra;
  float TimeStamp;
  int IsPCN;
} MDH;

int PrintMiniHeader(FILE *fp, MDH *mdh);
MDH *ReadMiniHeader(FILE *fp, MDH *mdh);
float MDHtimeStamp0(char *measoutpath);
int MDHnSamplesPerLine(char *measoutpath);
int MDHbytesPerLine(char *measoutpath);
int MDHbytesPerChunk(char *measoutpath);
int MDHnPCNs(char *measoutpath);
int MDHnPhaseEncodeLines(char *measoutpath);
int MDHnEchos(char *measoutpath);
int MDHnSlices(char *measoutpath);
int MDHnFrames(char *measoutpath);
int MDHslicePosition(char *measoutpath, int Slice, float *Sag, float *Cor, float *Tra);
float MDHsliceThickness(char *measoutpath);
float MDHsliceDirCos(char *measoutpath, float *dcSag, float *dcCor, float *dcTra);
float MDHreadTR(char *measoutpath);
float MDHreadTE(char *measoutpath, int nthEcho);
int MDHfastestDim(char *measoutpath);
char *MDHparseMrProt(char *file, char *TagString);
int MDHversion(char *measoutdir);
char *MDHascPath(char *measoutdir);
int AssignFakeDCs(MRI *mri, float dcSag, float dcCor, float dcTra);

/*--------------------------------------------------------------------*/
char *Progname = NULL;
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  isflag(char *flag);
static int  singledash(char *flag);
static int  stringmatch(char *str1, char *str2);
char *srcdir, *outdir;
int rev, infoonly=1, debug=0;

/*------------------------------------------------------------------*/
int main(int argc, char **argv)
{
  FILE *fp;
  long offset;
  float adc[10000];
  MDH *mdh = NULL;
  int n,d, err;
  int nSamplesPerLine, nPCNs, nEchos, nPELs, nSlices, nFrames, FastestDim;
  float Thickness, dcSag, dcCor, dcTra, TR, TE[500];
  float PhaseEncodeFOV, ReadoutFOV, FlipAngle;
  char *tmpstr;
  MRI *pcnr, *pcni;
  MRI **echor, **echoi;
  float *rptr, *iptr, TimeStamp0;
  int s,f,e,p,l;
  char fname[1000];
  char measoutpath[2000];
  int mdhversion;
  char *mdhascfile;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  sprintf(measoutpath,"%s/meas.out",srcdir);
  if(!fio_FileExistsReadable(measoutpath)){
    printf("ERROR: %s does not exist or is not readable\n",measoutpath);
    exit(1);
  }

  mdhversion = MDHversion(srcdir);
  if(mdhversion == -1){
    printf("ERROR: cannot find MDH ascii file\n");
    exit(1);
  }
  if(mdhversion != 15 && mdhversion != 21){
    printf("ERROR: MDH version %d unrecognized\n",mdhversion);
    exit(1);
  }
  if(mdhversion == 21){
    printf("ERROR: no support for version 21 yet\n");
    exit(1);
  }

  printf("MDH Version %d\n",mdhversion);
  mdhascfile = MDHascPath(srcdir);

  nSamplesPerLine = MDHnSamplesPerLine(measoutpath);
  nPCNs   = MDHnPCNs(measoutpath);
  nEchos  = MDHnEchos(measoutpath);
  nPELs   = MDHnPhaseEncodeLines(measoutpath);
  nSlices = MDHnSlices(measoutpath);
  nFrames = MDHnFrames(measoutpath);
  TR      = MDHreadTR(measoutpath);
  Thickness  = MDHsliceDirCos(measoutpath,&dcSag,&dcCor,&dcTra);
  FastestDim = MDHfastestDim(measoutpath);
  tmpstr = MDHparseMrProt(mdhascfile, "sSliceArray.asSlice[0].dPhaseFOV");
  sscanf(tmpstr,"%f",&PhaseEncodeFOV);
  free(tmpstr);
  tmpstr = MDHparseMrProt(mdhascfile, "sSliceArray.asSlice[0].dReadoutFOV");
  sscanf(tmpstr,"%f",&ReadoutFOV);
  free(tmpstr);
  tmpstr = MDHparseMrProt(mdhascfile, "dFlipAngleDegree");
  sscanf(tmpstr,"%f",&FlipAngle);
  free(tmpstr);

  printf("TR = %g, nSPL = %d, nPCNs = %d\n",TR,nSamplesPerLine,nPCNs);
  printf("nFrames = %d, nSlices = %d, nEchos = %d, nPELs = %d\n",
	 nFrames,nSlices,nEchos,nPELs);
  printf("Thick = %g, DC: %g, %g, %g\n",Thickness,dcSag,dcCor,dcTra);
  printf("Fastest Dim: %d\n",FastestDim);
  printf("FOV:  Phase = %g, Readout = %g\n",PhaseEncodeFOV,ReadoutFOV);

  for(n=0; n<nEchos; n++){
    TE[n] = MDHreadTE(measoutpath,n);
    printf("%2d  %8.4f\n",n,TE[n]);
  }

  if(infoonly) exit(0);

  /* Create the output directory */
  err = mkdir(outdir,(mode_t)-1);
  if(err != 0 && errno != EEXIST){
    printf("ERROR: creating %s\n",outdir);
    perror(NULL);    
    exit(1);
  }

  if(rev) printf("INFO: applying reversals\n");

  echor = (MRI**) calloc(sizeof(MRI*),nEchos);
  echoi = (MRI**) calloc(sizeof(MRI*),nEchos);
  printf("Allocating\n");
  for(n=0; n<nEchos; n++){
    echor[n] = MRIallocSequence(nSamplesPerLine,nPELs,nSlices,MRI_FLOAT,nFrames);

    echor[n]->xsize = ReadoutFOV/(nSamplesPerLine/2);
    echor[n]->ysize = PhaseEncodeFOV/nPELs;
    echor[n]->zsize = Thickness;
    echor[n]->tr    = TR;
    echor[n]->te    = TE[n];
    echor[n]->flip_angle  = FlipAngle*3.1415/180;

    /* DC is correct for z, bogus for x and y */
    /* Need to change so DC is RAS */
    AssignFakeDCs(echor[n], dcSag, dcCor, dcTra);

    echoi[n] = MRIallocSequence(nSamplesPerLine,nPELs,nSlices,MRI_FLOAT,nFrames);
    echoi[n]->xsize = ReadoutFOV/(nSamplesPerLine/2);
    echoi[n]->ysize = PhaseEncodeFOV/nPELs;
    echoi[n]->zsize = Thickness;
    echoi[n]->tr    = TR;
    echoi[n]->te    = TE[n];
    echoi[n]->flip_angle  = FlipAngle*3.1415/180;
    AssignFakeDCs(echoi[n], dcSag, dcCor, dcTra); /* bogus */

  }
  pcnr = MRIallocSequence(nSamplesPerLine,nPCNs,nSlices,MRI_FLOAT,nFrames);
  pcni = MRIallocSequence(nSamplesPerLine,nPCNs,nSlices,MRI_FLOAT,nFrames);

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  printf("offset = %ld\n",offset);
  fseek(fp,offset,SEEK_SET);

  printf("Loading\n");
  d = 0;
  for(f=0; f < nFrames; f++){
    if(f==0 || f%10 ==0) printf("f=%d\n",f);
    for(s=0; s < nSlices; s++){

      /* Phase Encodes Navigators */
      for(p=0; p < nPCNs; p++){
	fseek(fp,MDH_SIZE_V15,SEEK_CUR); /* skip header */
	fread(adc,sizeof(float), 2*nSamplesPerLine, fp);
	rptr = (float*) pcnr->slices[d][p];
	iptr = (float*) pcni->slices[d][p];
	if((p+1)%2 != 0 || !rev){
	  /* odd p */
	  for(n=0; n < 2*nSamplesPerLine; n += 2){
	    (*rptr++) = adc[n];
	    (*iptr++) = adc[n+1];
	  }
	}
	else{
	  /* even p (reverse) */
	  for(n = 2*nSamplesPerLine-2; n >= 0; n -= 2){
	    (*rptr++) = adc[n];
	    (*iptr++) = adc[n+1];
	  }
	}
      }

      /* For phase encode line as fastest dim */
      /* Need to reverse lines for odd echos? */
      for(e=0; e < nEchos; e++){
	for(l=0; l < nPELs; l++){
	  fseek(fp,MDH_SIZE_V15,SEEK_CUR); /* skip header */
	  fread(adc,sizeof(float), 2*nSamplesPerLine, fp);
	  rptr = (float*) echor[e]->slices[d][l];
	  iptr = (float*) echoi[e]->slices[d][l];
	  if((l+1)%2 != 0 || !rev){
	    /* odd line */
	    for(n=0; n < 2*nSamplesPerLine; n += 2){
	      (*rptr++) = adc[n];
	      (*iptr++) = adc[n+1];
	    }
	  }
	  else{
	    /* even l (reverse) */
	    for(n = 2*nSamplesPerLine-2; n >= 0; n -= 2){
	      (*rptr++) = adc[n];
	      (*iptr++) = adc[n+1];
	    }
	  }
	} /* phase encode line */
      } /* echo */

      d++;
    } /* slice */
  } /* frame */
  fclose(fp);

  printf("Saving\n");
  for(n=0; n<nEchos; n++){
    printf("Echo %d\n",n);
    sprintf(fname,"%s/echor%d.mgh",outdir,n);
    MRIwrite(echor[n],fname);
    MRIfree(&echor[n]);
    sprintf(fname,"%s/echoi%d.mgh",outdir,n);
    MRIwrite(echoi[n],fname);
    MRIfree(&echoi[n]);
  }
  sprintf(fname,"%s/pcnr.mgh",outdir);
  MRIwrite(pcnr,fname);
  sprintf(fname,"%s/pcni.mgh",outdir);
  MRIwrite(pcni,fname);

  printf("Done\n");

  return(0);
  /*-----------------------------------------------*/

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  printf("offset = %ld\n",offset);
  fseek(fp,offset,SEEK_SET);

  n = 1;
  while(!feof(fp)){
    mdh = ReadMiniHeader(fp,mdh);
    mdh->TimeStamp -= TimeStamp0;
    if(feof(fp)) break;
    printf("n = %d ---------------------------------\n",n);
    PrintMiniHeader(stdout,mdh);
    fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
    //fread(adc,sizeof(float), 2*mdh->Ncols, fp);
    n++;
  }

  fclose(fp);

  printf("n = %d\n",n);

  return(0);
}
/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  isflag(" "); /* shuts up compiler */

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
    else if (!strcasecmp(option, "--info"))    infoonly = 1;
    else if (!strcasecmp(option, "--rev"))     rev = 1;

    else if (stringmatch(option, "--srcdir")){
      if(nargc < 1) argnerr(option,1);
      srcdir = pargv[0];
      nargsused = 1;
    }
    else if (stringmatch(option, "--outdir")){
      if(nargc < 1) argnerr(option,1);
      outdir = pargv[0];
      infoonly = 0;
      nargsused = 1;
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
  printf("\n");
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --srcdir dir : directory with meas.out\n");
  printf("   --outdir dir : output directory\n");
  printf("   --rev        : apply reversals \n");
  //printf("   --info       : no ouput, just dump info (default with no outdir)\n");
  printf("\n");
  printf("   --help : a short story about a Spanish boy named Manual\n");
  printf("\n");
  //printf("   --svol svol.img (structural volume)\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
  printf("

Probes or converts a Siemens proprietary format (MDH) file to MGH format. 

If an output directory is not specified, then the MDH file is probed,
and the various parameters are printed to stdout.

When an output is specified, the output will have two parts: phase
correction navigators and echos.  The phase correction navigators will
be stored in pcnr.mgh (real) and pcni.mgh (imaginary).  The echos will
be stored in echoNr.mgh (real) and echoNi.mgh (imaginary), where N is
the echo number (starting with 0).

--srcdir dir

Directory with the meas.out file and the ascii file. The meas.out file
contains the measured k-space data in Siemens MDH format (which includes
miniheaders). The meas.out file must be accompanied by an ascii file
called either MrProt.asc or mrprot.asc (NUMARIS 4 VA15) or meas.asc
(NUMARIS 4 VA21). The version is determined by the name of this file.

--outdir dir

Directory where the ouput will be stored. If this directory does not
exist, it will be created (only for one level). 

--rev

Reverse the readouts for even lines. Bug: should reverse lines for
even number echos, but does not do this yet.

BUGS:

Lines of even-numbered echos are not reversed with --rev.

Does not handle version V21 yet.

Does not handle 3D sequences yet.

Does not handle FLASH  sequences yet.

Does not handle multiple channels yet.

  ");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if(n==1) fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else     fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(srcdir == NULL){
    printf("ERROR: no source directory\n");
    exit(1);
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  return;
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
/*---------------------------------------------------------------*/
static int isflag(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*------------------------------------------------------------*/
static int stringmatch(char *str1, char *str2)
{
  if(! strcmp(str1,str2)) return(1);
  return(0);
}















/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
/*----------------------------------------------------------
  MDHversion() - determines the NUMARIS version number given
  the directory where the meas.out lives. The version is
  determined by the name of the ascii file.
  ----------------------------------------------------------*/
int MDHversion(char *measoutdir)
{
  char fname[2000];

  sprintf(fname,"%s/MrProt.asc",measoutdir);
  if(fio_FileExistsReadable(fname))  return(15);

  sprintf(fname,"%s/mrprot.asc",measoutdir);
  if(fio_FileExistsReadable(fname))  return(15);

  sprintf(fname,"%s/meas.asc",measoutdir);
  if(fio_FileExistsReadable(fname))  return(21);

  return(-1);
}
/*----------------------------------------------------------
  MDHascPath() - returns the path to the MDH ascii file. This
  will be either measoutdir/mrprot.asc (for NUMARIS 1.5) or
  measoutdir/MrProt.asc (also 1.5) or measoutdir/meas.asc 
  (for NUMARIS 2.1). Returns NULL if this file cannot be
  found.
  ----------------------------------------------------------*/
char *MDHascPath(char *measoutdir)
{
  char fname[2000];
  char *mrprot;
  int len;

  sprintf(fname,"%s/MrProt.asc",measoutdir);
  if(fio_FileExistsReadable(fname)){
    len = strlen(fname);
    mrprot = (char *) calloc(sizeof(char),len+1);
    memcpy(mrprot,fname,len);
    return(mrprot);
  }

  sprintf(fname,"%s/mrprot.asc",measoutdir);
  if(fio_FileExistsReadable(fname)){
    len = strlen(fname);
    mrprot = (char *) calloc(sizeof(char),len+1);
    memcpy(mrprot,fname,len);
    return(mrprot);
  }

  sprintf(fname,"%s/meas.asc",measoutdir);
  if(fio_FileExistsReadable(fname)){
    len = strlen(fname);
    mrprot = (char *) calloc(sizeof(char),len+1);
    memcpy(mrprot,fname,len);
    return(mrprot);
  }

  return(NULL);
}


/*----------------------------------------------------------*/
MDH *ReadMiniHeader(FILE *fp, MDH *mdh)
{
  unsigned long ultmp;
  unsigned short ustmp;
  float ftmp;

  if(mdh == NULL)  mdh = (MDH *) calloc(sizeof(MDH),1);
  
  fread(&ultmp,sizeof(long),1, fp); //DMALength
  fread(&ultmp,sizeof(long),1, fp); //MeasUID (long)
  fread(&ultmp,sizeof(long),1, fp); //ScanCounter
  fread(&mdh->ulTimeStamp,sizeof(long),1, fp);
  mdh->TimeStamp = 2.5*mdh->ulTimeStamp; // dont know why its 2.5
  fread(&ultmp,sizeof(long),1, fp); //PMUTimeStamp
  fread(&mdh->BitMask1,sizeof(long),1, fp); 
  //fread(&ultmp,sizeof(long),1, fp); //BitMask2 (VA21)

  fread(&mdh->Ncols,sizeof(short),1, fp); //SamplesInScan
  fread(&mdh->UsedChannels,sizeof(short),1, fp); 
  fread(&mdh->LoopCounterLine,sizeof(short),1, fp); //LoopCounter.Line
  fread(&ustmp,sizeof(short),1, fp); //Acquistion
  fread(&mdh->Slice,sizeof(short),1, fp); 
  fread(&mdh->Partition,sizeof(short),1, fp);
  fread(&mdh->Echo,sizeof(short),1, fp); 
  fread(&ustmp,sizeof(short),1, fp); // Phase
  fread(&mdh->Rep,sizeof(short),1, fp); 
  fread(&ustmp,sizeof(short),1, fp); //Set
  fread(&ustmp,sizeof(short),1, fp); //Seg
  //fread(&ustmp,sizeof(short),1, fp); //Ida (VA21)
  //fread(&ustmp,sizeof(short),1, fp); //Idb (VA21)
  //fread(&ustmp,sizeof(short),1, fp); //Idc (VA21)
  //fread(&ustmp,sizeof(short),1, fp); //Idd (VA21)
  //fread(&ustmp,sizeof(short),1, fp); //Ide (VA21)
  fread(&ustmp,sizeof(short),1, fp); //LoopCounterFree (VA15)
  fread(&ustmp,sizeof(short),1, fp); //CutOffData.Pre
  fread(&ustmp,sizeof(short),1, fp); //CutOffData.Post
  fread(&ustmp,sizeof(short),1, fp); //KSpaceCenterColumn
  fread(&ustmp,sizeof(short),1, fp); //Dummy

  fread(&ftmp,sizeof(float),1, fp); //ReadOutOffCenter
  fread(&ultmp,sizeof(long),1, fp); //TimeSinceLastRF
  
  fread(&ustmp,sizeof(short),1, fp); //KSpaceCentreLineNo
  fread(&ustmp,sizeof(short),1, fp); //KSpaceCentrePartitionNo
  //fread(&ustmp,sizeof(short),1, fp); //IceProgramPara1 (VA21)
  //fread(&ustmp,sizeof(short),1, fp); //IceProgramPara2 (VA21)
  //fread(&ustmp,sizeof(short),1, fp); //IceProgramPara3 (VA21)
  //fread(&ustmp,sizeof(short),1, fp); //IceProgramPara4 (VA21)
  fseek(fp,14*sizeof(short),SEEK_CUR); // FreePara (VA15)
  //fseek(fp,4*sizeof(short),SEEK_CUR); // FreePara (VA21)

  fread(&mdh->SlicePosSag,sizeof(float),1, fp);
  fread(&mdh->SlicePosCor,sizeof(float),1, fp);
  fread(&mdh->SlicePosTra,sizeof(float),1, fp);

  fread(&ftmp,sizeof(float),1, fp); //Quaternion1
  fread(&ftmp,sizeof(float),1, fp); //Quaternion2
  fread(&ftmp,sizeof(float),1, fp); //Quaternion3
  fread(&ftmp,sizeof(float),1, fp); //Quaternion4

  fread(&ultmp,sizeof(long),1, fp); //ChannelId

  if( (mdh->BitMask1 & MDH_BM_PHASECOR) != 0) mdh->IsPCN = 1;
  else  mdh->IsPCN = 0;

  return(mdh);
}
/*---------------------------------------------------*/
int PrintMiniHeader(FILE *fp, MDH *mdh)
{
  fprintf(fp,"TimeStamp    %f\n",mdh->TimeStamp);
  fprintf(fp,"BitMask1     %lx\n",mdh->BitMask1);
  fprintf(fp,"IsPCN        %d\n",mdh->IsPCN);
  fprintf(fp,"Reflect      %ld\n",mdh->BitMask1 & MDH_BM_REFLECT);
  fprintf(fp,"LoopCounter   %d\n",mdh->LoopCounterLine);
  fprintf(fp,"Ncols        %d\n",mdh->Ncols);
  fprintf(fp,"UsedChannels %d\n",mdh->UsedChannels);
  fprintf(fp,"Slice        %d\n",mdh->Slice);
  fprintf(fp,"Echo         %d\n",mdh->Echo);
  fprintf(fp,"Rep          %d\n",mdh->Rep);
  fprintf(fp,"SPSag        %f\n",mdh->SlicePosSag);
  fprintf(fp,"SPCor        %f\n",mdh->SlicePosCor);
  fprintf(fp,"SPTra        %f\n",mdh->SlicePosTra);
  return(0);
}
/*---------------------------------------------------*/
float MDHtimeStamp0(char *measoutpath)
{
  MDH *mdh=NULL;
  FILE *fp;
  float TimeStamp0;
  unsigned long offset;
  
  fp = fopen(measoutpath,"r");
  if(fp == NULL){
    printf("ERROR: could not open %s\n",measoutpath);
    return(-10000000);
  }
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);
  mdh = ReadMiniHeader(fp,mdh);
  TimeStamp0 = mdh->TimeStamp;
  free(mdh);
  fclose(fp);
  return(TimeStamp0);
}
/*--------------------------------------------------------------------
  MDHnSamplesPerLine() - returns number of samples collected for each
  phase encode line based on the number of samples in the
  phase correction naviagator. The number of values stored in the
  line is twice this number to account for real and imaginary.
  This may need to change for multiple channels.
  -------------------------------------------------------------------*/
int MDHnSamplesPerLine(char *measoutpath)
{
  MDH *mdh=NULL;
  FILE *fp;
  int nSamplesPerLine;
  unsigned long offset;
  
  fp = fopen(measoutpath,"r");
  if(fp == NULL){
    printf("ERROR: could not open %s\n",measoutpath);
    return(-10000000);
  }
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);
  mdh = ReadMiniHeader(fp,mdh);
  nSamplesPerLine = mdh->Ncols;
  free(mdh);
  fclose(fp);
  return(nSamplesPerLine);
}
/*--------------------------------------------------------------------
  MDHbytesPerLine() - returns number of bytes associated with
  each phase encode line based on nSamplesPerLine.  This may need to 
  change for multiple channels.
  -------------------------------------------------------------------*/
int MDHbytesPerLine(char *measoutpath)
{
  int BytesPerLine;
  // 2 = real and imaginary
  // 4 = bytes per float
  BytesPerLine = 2*4*MDHnSamplesPerLine(measoutpath);
  return(BytesPerLine);
}
/*-----------------------------------------------------------------------
  MDHbytesPerChunk() - returns the total number of bytes associated
  with each readout, including the header and data. Need to autodetect
  version.
  --------------------------------------------------------------------------*/
int MDHbytesPerChunk(char *measoutpath)
{
  int BytesPerChunk;
  BytesPerChunk = MDHbytesPerLine(measoutpath) + MDH_SIZE_V15;
  return(BytesPerChunk);
}
/*-----------------------------------------------------------------------
  MDHnPCNs() - returns the number of phase correction navigators
  collected per slice based on the number in the first slice.
  --------------------------------------------------------------------------*/
int MDHnPCNs(char *measoutpath)
{
  MDH *mdh=NULL;
  FILE *fp;
  int nPCNs;
  unsigned long offset;
  
  fp = fopen(measoutpath,"r");
  if(fp == NULL){
    printf("ERROR: could not open %s\n",measoutpath);
    return(-10000000);
  }
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  nPCNs = 0;
  while(1){
    mdh = ReadMiniHeader(fp,mdh);
    fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
    if(! mdh->IsPCN) break;
    nPCNs ++;
  }

  free(mdh);
  fclose(fp);
  return(nPCNs);
}

/*-----------------------------------------------------------------------
  MDHnSlices() - returns the number of slices based on the
  number of slices found in the first repetition.
  --------------------------------------------------------------------------*/
int MDHnSlices(char *measoutpath)
{
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nSlices, nPCNs, nEchos, nPELs, chunksize;
  int BytesPerSlice, BytesPerSliceB;
  
  chunksize = MDHbytesPerChunk(measoutpath);
  if(chunksize < 0) return(-1);
  nPCNs   = MDHnPCNs(measoutpath);
  nEchos = MDHnEchos(measoutpath);
  nPELs  = MDHnPhaseEncodeLines(measoutpath);

  // Bytes per slice, including headers and PCNs
  BytesPerSlice  = chunksize*(nPCNs + nEchos*nPELs);
  // Bytes per slice, minus one header size
  BytesPerSliceB = BytesPerSlice-MDH_SIZE_V15;

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  nSlices = 0;
  while(1){
    mdh = ReadMiniHeader(fp,mdh);
    if(feof(fp)) break;
    if(mdh->Rep != 0) break;
    if(nSlices < mdh->Slice+1) nSlices = mdh->Slice+1;
    // Jump to the beginning of the next slice
    fseek(fp,BytesPerSliceB,SEEK_CUR);
  }

  free(mdh);
  fclose(fp);
  return(nSlices);
}

/*-----------------------------------------------------------------------
  MDHslicePosition() - returns the slice position for the given zero-based
  slice number. Returns 0 is ok, otherwise returns 1.
  --------------------------------------------------------------------------*/
int MDHslicePosition(char *measoutpath, int Slice, float *Sag, float *Cor, float *Tra)
{
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nSlices, nPCNs, nEchos, nPELs, chunksize;
  int BytesPerSlice;
  
  chunksize = MDHbytesPerChunk(measoutpath);
  if(chunksize < 0) return(1);
  nPCNs   = MDHnPCNs(measoutpath);
  nEchos = MDHnEchos(measoutpath);
  nPELs  = MDHnPhaseEncodeLines(measoutpath);
  nSlices = MDHnSlices(measoutpath);

  if(Slice >= nSlices){
    printf("ERROR: MDHslicePosition: requested slice %d exceeds max %d\n",
	   Slice,nSlices-1);
    return(1);
  }

  // Bytes per slice, including headers and PCNs
  BytesPerSlice  = chunksize*(nPCNs + nEchos*nPELs);

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Jump to the beginning of the desired slice
  fseek(fp,Slice*BytesPerSlice,SEEK_CUR);
  mdh = ReadMiniHeader(fp,mdh);

  *Sag = mdh->SlicePosSag;
  *Cor = mdh->SlicePosCor;
  *Tra = mdh->SlicePosTra;

  free(mdh);
  fclose(fp);
  return(0);
}
/*-----------------------------------------------------------------------
  MDHsliceThickness - returns the distance between two adjacent slices
  --------------------------------------------------------------------------*/
float MDHsliceThickness(char *measoutpath)
{
  float dcSag, dcCor, dcTra;
  float Thickness;
  Thickness = MDHsliceDirCos(measoutpath,&dcSag,&dcCor,&dcTra);
  return(Thickness);
}
/*-----------------------------------------------------------------------
  MDHsliceDirCos - computes the direction cosines as pointing from the 0th
  slice to the 1st slice. Return value is the slice thickeness.
  --------------------------------------------------------------------------*/
float MDHsliceDirCos(char *measoutpath, float *dcSag, float *dcCor, float *dcTra)
{
  float Sag0, Cor0, Tra0;
  float Sag1, Cor1, Tra1;
  float Thickness;
  float dSag, dCor, dTra;

  MDHslicePosition(measoutpath, 0, &Sag0, &Cor0, &Tra0);
  MDHslicePosition(measoutpath, 1, &Sag1, &Cor1, &Tra1);

  dSag = (Sag1-Sag0);
  dCor = (Cor1-Cor0);
  dTra = (Tra1-Tra0);

  Thickness = sqrt( dSag*dSag + dCor*dCor + dTra*dTra);
  *dcSag = dSag/Thickness;
  *dcCor = dCor/Thickness;
  *dcTra = dTra/Thickness;
  return(Thickness);
}


/*-----------------------------------------------------------------------
  MDHnFrames() - returns the number of frames/repetitions
  --------------------------------------------------------------------------*/
int MDHnFrames(char *measoutpath)
{
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nFrames, nSlices, nPCNs, nEchos, nPELs, chunksize;
  int BytesPerFrame, BytesPerFrameB;
  
  chunksize = MDHbytesPerChunk(measoutpath);
  if(chunksize < 0) return(-1);
  nPCNs    = MDHnPCNs(measoutpath);
  nEchos  = MDHnEchos(measoutpath);
  nPELs   = MDHnPhaseEncodeLines(measoutpath);
  nSlices = MDHnSlices(measoutpath);

  // Bytes per frame, including headers and PCNs
  BytesPerFrame  = chunksize*nSlices*(nPCNs + nEchos*nPELs);
  // Bytes per frame, minus one header size
  BytesPerFrameB = BytesPerFrame-MDH_SIZE_V15;

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  nFrames = 0;
  while(1){
    mdh = ReadMiniHeader(fp,mdh);
    if(feof(fp)) break;
    if(nFrames < mdh->Rep+1) nFrames = mdh->Rep+1;
    // Jump to the beginning of the next frame/rep
    fseek(fp,BytesPerFrameB,SEEK_CUR);
  }

  free(mdh);
  fclose(fp);
  return(nFrames);
}

/*-----------------------------------------------------------------------
  MDHnEchos() - returns the number of echos as found in the first slice.
  --------------------------------------------------------------------------*/
int MDHnEchos(char *measoutpath)
{
  MDH *mdh=NULL;
  FILE *fp;
  int nEchos;
  unsigned long offset;
  int nPCNs, chunksize;
  
  chunksize = MDHbytesPerChunk(measoutpath);
  if(chunksize < 0) return(-1);
  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Get past the phase correction navigators
  nPCNs = MDHnPCNs(measoutpath);
  fseek(fp,nPCNs*chunksize,SEEK_CUR);

  nEchos = -1;
  while(1){
    mdh = ReadMiniHeader(fp,mdh);
    if(feof(fp)) break;
    if(mdh->Slice != 0) break;
    if(nEchos < mdh->Echo+1) nEchos = mdh->Echo+1;
    fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
  }

  free(mdh);
  fclose(fp);
  return(nEchos);
}
/*-----------------------------------------------------------------------
  MDHnPhaseEncodeLines() - returns the number of phase encode lines as 
  found in the first slice.
  --------------------------------------------------------------------------*/
int MDHnPhaseEncodeLines(char *measoutpath)
{
  MDH *mdh=NULL;
  FILE *fp;
  int nPhaseEncodeLines;
  unsigned long offset;
  int nPCNs, chunksize;
  
  chunksize = MDHbytesPerChunk(measoutpath);
  if(chunksize < 0) return(-1);
  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Get past the phase correction navigators
  nPCNs = MDHnPCNs(measoutpath);
  fseek(fp,nPCNs*chunksize,SEEK_CUR);

  nPhaseEncodeLines = -1;
  while(1){
    mdh = ReadMiniHeader(fp,mdh);
    if(feof(fp)) break;
    if(mdh->Slice != 0) break;
    if(nPhaseEncodeLines < mdh->LoopCounterLine+1) nPhaseEncodeLines = mdh->LoopCounterLine+1;
    fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
  }

  free(mdh);
  fclose(fp);
  return(nPhaseEncodeLines);
}
/*-------------------------------------------------------------------------
  MDHreadTR() - read the TR as the difference between the time stamps 
  two adjacent repetitions. Value is in seconds. If there is only one
  frame, a value of 0 is returned.
  -------------------------------------------------------------------------*/
float MDHreadTR(char *measoutpath)
{
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nFrames, nSlices, nPCNs, nEchos, nPELs, chunksize;
  int BytesPerFrame, BytesPerFrameB;
  float TimeStamp0, TimeStamp1, TR;
  
  chunksize = MDHbytesPerChunk(measoutpath);
  if(chunksize < 0) return(-1);
  nPCNs    = MDHnPCNs(measoutpath);
  nEchos  = MDHnEchos(measoutpath);
  nPELs   = MDHnPhaseEncodeLines(measoutpath);
  nSlices = MDHnSlices(measoutpath);
  nFrames = MDHnFrames(measoutpath);
  if(nFrames == 1) return(0);

  // Bytes per frame, including headers and PCNs
  BytesPerFrame  = chunksize*nSlices*(nPCNs + nEchos*nPELs);
  // Bytes per frame, minus one header size
  BytesPerFrameB = BytesPerFrame-MDH_SIZE_V15;

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Read time stamp from first frame/rep
  mdh = ReadMiniHeader(fp,mdh);
  TimeStamp0 = mdh->TimeStamp;

  // Jump to the beginning of the next frame/rep
  fseek(fp,BytesPerFrameB,SEEK_CUR);

  // Read time stamp from second frame/rep
  mdh = ReadMiniHeader(fp,mdh);
  TimeStamp1 = mdh->TimeStamp;

  TR = (TimeStamp1-TimeStamp0)/1000.0;

  free(mdh);
  fclose(fp);
  return(TR);
}
/*----------------------------------------------------------------------
  MDHfastestDim(char *measoutpath) - determines whether the
  fastest dim is the echo loop or the phase encode loop. For EPIs,
  it will be the phase encode loop. For FLASH, it will be the echo loop.
  This is determined by examining the first two non-PCN headers to
  see whether the LoopCounterLine changes or the Echo number changes.
  Return value = 1 if Echo is fastest
  Return value = 2 if Phase Encode is fastest
  ----------------------------------------------------------------------*/
int MDHfastestDim(char *measoutpath)
{
  MDH *mdh1, *mdh2;
  FILE *fp;
  unsigned long offset;
  int nPCNs, FastestDim, chunksize;
  
  chunksize = MDHbytesPerChunk(measoutpath);
  if(chunksize < 0) return(-1);
  nPCNs = MDHnPCNs(measoutpath);

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Get past PCNs
  fseek(fp,nPCNs*chunksize,SEEK_CUR);

  // Read the first non-PCN header
  mdh1 = ReadMiniHeader(fp,NULL);
  fseek(fp,2*mdh1->Ncols*sizeof(float),SEEK_CUR);

  // Read the second non-PCN header
  mdh2 = ReadMiniHeader(fp,NULL);

  if(mdh1->Echo != mdh2->Echo)  FastestDim = 1; // Echo is fastest
  else                          FastestDim = 2; // Phase Encode is fastest

  free(mdh1);
  free(mdh2);
  fclose(fp);

  return(FastestDim);
}

/*-----------------------------------------------------------------------
  MDHreadTE() - computes the TE of the nth Echo as the midpoint between the 
  start of the echo and the end of the echo. nthEcho is zero based. The
  TE is in milliseconds.
  ---------------------------------------------------------------------*/
float MDHreadTE(char *measoutpath, int nthEcho)
{
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nPCNs, nEchos, nPELs, chunksize, FastestDim;
  float TimeStamp0, TimeStamp1, TimeStamp1a, TimeStamp1b, TE;
  
  chunksize = MDHbytesPerChunk(measoutpath);
  if(chunksize < 0) return(-1);
  nPCNs    = MDHnPCNs(measoutpath);
  nEchos  = MDHnEchos(measoutpath);
  if(nthEcho >= nEchos){
    printf("ERROR: requested TE=%d, max is %d\n",nthEcho,nEchos);
    return(-1);
  }
  nPELs = MDHnPhaseEncodeLines(measoutpath);
  FastestDim = MDHfastestDim(measoutpath);

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Get time stamp of first PCN as that closest to RF pulse
  mdh = ReadMiniHeader(fp,mdh);
  fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
  TimeStamp0 = mdh->TimeStamp;
  // Get past remaining PCNs
  fseek(fp,(nPCNs-1)*chunksize,SEEK_CUR);

  if(FastestDim == 1){
    // Echo is the fastest (FLASH)
    // TE is time stamp at nthEcho. Is this correct?
    // This branch has never been tested!!!!!

    // Seek to the nthEcho
    fseek(fp,nthEcho*chunksize,SEEK_CUR);
    mdh = ReadMiniHeader(fp,mdh);
    TimeStamp1 = mdh->TimeStamp;

  }
  else {
    // Phase Encode is fastest (EPI)
    // TE is the average between the first and last PEL

    // Seek to first PEL for this echo
    fseek(fp,nthEcho*nPELs*chunksize,SEEK_CUR);
    mdh = ReadMiniHeader(fp,mdh);
    fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
    TimeStamp1a = mdh->TimeStamp;
    
    // Seek to last PEL for this echo
    fseek(fp,(nPELs-2)*chunksize,SEEK_CUR);
    mdh = ReadMiniHeader(fp,mdh);
    TimeStamp1b = mdh->TimeStamp;
    
    TimeStamp1 = (TimeStamp1a + TimeStamp1b)/2;
  }

  TE = TimeStamp1-TimeStamp0;

  free(mdh);
  fclose(fp);
  return(TE);
}
/*-----------------------------------------------------------------
  MDHparseMrProt() - gets the value of the given tag from the file.
  The format is assumed to be that of the Siemens MrProt. The MrProt
  is embedded as ascii text in the Siemens DICOM file, so the file may 
  be an MrProt file or a Siemens DICOM file. 

  In any case, the MrProt ascii block is begun by the line
        ### ASCCONV BEGIN ###
  and ended by the line
        ### ASCCONV END ###

  Each line inside the block has the form:
    VariableName = VariableValue

  This function searches this block for a variable named TagString
  and returns the Value as a string.

  It returns NULL if 
    1. The begining of the ASCII block cannot be found
    2. There is no match with the TagString

  Author: Douglas N. Greve, 5/24/2003
  -----------------------------------------------------------------*/
char *MDHparseMrProt(char *file, char *TagString)
{
  char linestr[1000]; 
  char tmpstr2[500];
  FILE *fp;
  int dumpline, nthchar;
  char *rt;
  char *BeginStr;
  int LenBeginStr;
  char *TestStr;
  int nTest;
  char VariableName[500];
  char *VariableValue;
  
  BeginStr = "### ASCCONV BEGIN ###";
  LenBeginStr = strlen(BeginStr);
  TestStr = (char *) calloc(LenBeginStr+1,sizeof(char));

  fp = fopen(file,"r");
  if(fp == NULL){
    printf("ERROR: could not open dicom file %s\n",file);
    exit(1);
  }

  /* This section steps through the file char-by-char until
     the BeginStr is matched */
  dumpline = 0;
  nthchar = 0;
  while(1){
    fseek(fp,nthchar, SEEK_SET);
    nTest = fread(TestStr,sizeof(char),LenBeginStr,fp);
    if(nTest != LenBeginStr) break;
    if(strcmp(TestStr,BeginStr)==0){
      fseek(fp,nthchar, SEEK_SET);
      dumpline = 1;
      break;
    }
    nthchar ++;
  }
  free(TestStr);

  
  if(! dumpline) return(NULL); /* Could not match Begin String */

  /* Once the Begin String has been matched, this section 
     searches each line until the TagString is matched
     or until the End String is matched */
  VariableValue = NULL;
  while(1){
    rt = fgets(linestr,1000,fp);
    if(rt == NULL) break;

    if(strncmp(linestr,"### ASCCONV END ###",19)==0) break;

    sscanf(linestr,"%s %*s %*s",VariableName);

    if(strlen(VariableName) != strlen(TagString)) continue;

    if( strcmp(VariableName,TagString)==0 ) { 
      /* match found */
      sscanf(linestr,"%*s %*s %s",tmpstr2);
      VariableValue = (char *) calloc(strlen(tmpstr2)+1,sizeof(char));
      memcpy(VariableValue, tmpstr2, strlen(tmpstr2));
      break;
    }
  }
  fclose(fp);

  //printf("Leaving SiemensAsciiTag() \n");fflush(stdout);fflush(stderr);
  fflush(stdout);fflush(stderr);

  return(VariableValue);
}
/*----------------------------------------------------------------
  AssignFakeDCs() - assigns direction cosines that are orthogonal.
  The z DC is correct, but the x and y are bogus.
  ---------------------------------------------------------------*/
int AssignFakeDCs(MRI *mri, float dcSag, float dcCor, float dcTra)
{
  float norm;
  /* DC is correct for z, bogus for x and y */
  /* Need to change so DC is RAS */
  mri->z_r   = dcSag; 
  mri->z_a   = dcCor;
  mri->z_s   = dcTra;

  mri->x_r   = 0;
  mri->x_a   = -dcTra;
  mri->x_s   = +dcCor;
  norm = sqrt(mri->x_r*mri->x_r + mri->x_a*mri->x_a + mri->x_s*mri->x_s);
  mri->x_r /= norm;
  mri->x_a /= norm;
  mri->x_s /= norm;

  mri->y_r   =  dcCor*dcCor + dcTra*dcTra;
  mri->y_a   = -dcSag*dcCor;
  mri->y_s   = -dcSag*dcTra;
  norm = sqrt(mri->y_r*mri->y_r + mri->y_a*mri->y_a + mri->y_s*mri->y_s);
  mri->y_r /= norm;
  mri->y_a /= norm;
  mri->y_s /= norm;

  return(0);
}
