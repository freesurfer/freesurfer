/**
 * @file  mri_convert_mdh.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:14 $
 *    $Revision: 1.27 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


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
#include "version.h"

#ifndef lint
static char vcid[] = "$Id: mri_convert_mdh.c,v 1.27 2011/03/02 00:04:14 nicks Exp $";
#endif /* lint */

#define MDH_SIZE    128        //Number of bytes in the miniheader
#define MDH_BM_ONLINE   (1 <<  4)
#define MDH_BM_PHASECOR (1 << 21)
#define MDH_BM_REFLECT  (1 << 24)
#define MDH_DIM_FRAME   1
#define MDH_DIM_SLICE   2
#define MDH_DIM_LINE    3
#define MDH_DIM_ECHO    4

typedef struct tagMDH {
  unsigned long ulTimeStamp, BitMask1, ScanCounter, ChannelId;
  unsigned short Ncols, UsedChannels, Slice, Partition, Echo, Rep,LoopCounterLine;
  unsigned short KSpaceCenterCol,KSpaceCenterLine, CutOffDataPre, CutOffDataPost ;
  float SlicePosSag,SlicePosCor,SlicePosTra;
  float TimeStamp, PED, ReadoutOffCenter;
  int IsPCN;
}
MDH;

#define MDH_ACQEND            (1 <<  0)
#define MDH_RTFEEDBACK        (1 <<  1)
#define MDH_HPFEEDBACK        (1 <<  2)
#define MDH_ONLINE            (1 <<  3)
#define MDH_OFFLINE           (1 <<  4)
#define MDH_SYNCDATA          (1 <<  5)   // readout contains synchronous data
#define MDH_6                 (1 <<  6)
#define MDH_7                 (1 <<  7)
#define MDH_LASTSCANINCONCAT  (1 <<  8)   // Flag for last scan in concatenation
#define MDH_9                 (1 <<  9)
#define MDH_RAWDATACORRECTION (1 << 10)  // Correct the rawdata with the rawdata correction factor
#define MDH_LASTSCANINMEAS    (1 << 11)  // Flag for last scan in measurement
#define MDH_SCANSCALEFACTOR   (1 << 12)  // Flag for scan specific additional scale factor
#define MDH_2NDHADAMARPULSE   (1 << 13)  // 2nd RF excitation of HADAMAR
#define MDH_REFPHASESTABSCAN  (1 << 14)  // reference phase stabilization scan
#define MDH_PHASESTABSCAN     (1 << 15)  // phase stabilization scan
#define MDH_D3FFT             (1 << 16)  // execute 3D FFT
#define MDH_SIGNREV           (1 << 17)  // sign reversal
#define MDH_PHASEFFT          (1 << 18)  // execute phase fft
#define MDH_SWAPPED           (1 << 19)  // swapped phase/readout direction
#define MDH_POSTSHAREDLINE    (1 << 20)  // shared line
#define MDH_PHASECOR          (1 << 21)  // phase correction data
#define MDH_PATREFSCAN        (1 << 22)  // additional scan for PAT reference line/partition
#define MDH_PATREFANDIMASCAN  (1 << 23)  // additional scan for PAT reference line/partition also used as image scan
#define MDH_REFLECT           (1 << 24)  // reflect line
#define MDH_NOISEADJSCAN      (1 << 25)  // noise adjust scan --> Not used in NUM4
#define MDH_SHARENOW          (1 << 26)  // all lines are acquired from the actual and previous e.g. phases
#define MDH_LASTMEASUREDLINE  (1 << 27)  // current line is the last measured line of all succeeding e.g. phases
#define MDH_FIRSTSCANINSLICE  (1 << 28)  // indicates first scan in slice (needed for time stamps)
#define MDH_LASTSCANINSLICE   (1 << 29)  // indicates last scan in slice (needed for time stamps)
#define MDH_TREFFECTIVEBEGIN  (1 << 30)  // indicates the begin time stamp for TReff (triggered measurement)
#define MDH_TREFFECTIVEEND    (1 << 31)  // indicates the end time stamp for TReff (triggered measurement)

// Not sure about these 32
#define MDH_32
#define MDH_33
#define MDH_34
#define MDH_35
#define MDH_36
#define MDH_37
#define MDH_38
#define MDH_39
#define MDH_FIRST_SCAN_IN_BLADE       (1 <<  9)  // Marks the first line of a blade
#define MDH_LAST_SCAN_IN_BLADE        (1 << 10)  // Marks the last line of a blade
#define MDH_LAST_BLADE_IN_TR          (1 << 11)  // Set for all lines of the last BLADE in each TR interval
#define MDH_43
#define MDH_44
#define MDH_RETRO_LASTPHASE           (1 << 12)  // Marks the last phase in a heartbeat
#define MDH_RETRO_ENDOFMEAS           (1 << 13)  // Marks an ADC at the end of the measurement
#define MDH_RETRO_REPEATTHISHEARTBEAT (1 << 14)  // Repeat the current heartbeat when this bit is found
#define MDH_RETRO_REPEATPREVHEARTBEAT (1 << 15)  // Repeat the previous heartbeat when this bit is found
#define MDH_RETRO_ABORTSCANNOW        (1 << 16)  // Just abort everything
#define MDH_RETRO_LASTHEARTBEAT       (1 << 17)  // This adc is from the last heartbeat (a dummy)
#define MDH_RETRO_DUMMYSCAN           (1 << 18)  // This adc is just a dummy scan, throw it away
#define MDH_RETRO_ARRDETDISABLED      (1 << 19)  // Disable all arrhythmia detection when this bit is found


// This should all add up to 128
typedef struct tagMDH_VB13 {
  unsigned int   ulFlagsAndDMALength;
  int             lMeasUID;
  unsigned int   ulScanCounter;
  unsigned int   ulTimeStamp;     // Multiply by 2.5 to get msec
  unsigned int   ulPMUTimeStamp;  // Multiply by 2.5 to get msec
  unsigned int  aulEvalInfoMask[2];
  unsigned short ushSamplesInScan ;
  unsigned short ushUsedChannels;
  unsigned short ushLine;
  unsigned short ushAcquisition; // note: acquisition is same as average (JP)
  unsigned short ushSlice;
  unsigned short ushPartition;;
  unsigned short ushEcho;
  unsigned short ushPhase;
  unsigned short ushRepetition;
  unsigned short ushSet;
  unsigned short ushSeg;
  unsigned short ushIda;
  unsigned short ushIdb;
  unsigned short ushIdc;
  unsigned short ushIdd;
  unsigned short ushIde;

  unsigned short ushPre;
  unsigned short ushPost;
  unsigned short ushKSpaceCentreColumn;
  unsigned short ushDummy;

  float            fReadOutOffcenter;
  unsigned int    ulTimeSinceLastRF;

  unsigned short  ushKSpaceCentreLineNo;
  unsigned short  ushKSpaceCentrePartitionNo;
  unsigned short aushIceProgramPara[4];
  unsigned short aushFreePara[4];

  float             fSag;
  float             fCor;
  float             fTra;
  float            afQuaternion[4];

  unsigned short  ushChannelId;
  unsigned short  ushPTABPosNeg;
}
MDH_VB13;

int PrintMDH_VB13r(FILE *fp, MDH_VB13 *mdh);
int DumpMDH_VB13(char *measfile);
int ConvertMDH_VB13(char *measfile, char *outbase);

int PrintMiniHeader(FILE *fp, MDH *mdh);
MDH *ReadMiniHeader(FILE *fp, MDH *mdh, int mdhversion);
MDH *ReadMiniHeader15(FILE *fp, MDH *mdh);
MDH *ReadMiniHeader21(FILE *fp, MDH *mdh);
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
int MDHversion(char *measoutpath);
char *MDHascPath(char *measoutdir);
int AssignFakeDCs(MRI *mri, float dcSag, float dcCor, float dcTra);
int MDHdump(char *measoutpath);
int MDHadcStats(char *measoutpath);
int *MDHind2sub(int *siz, int ind, int ndim, int *sub);
int MDHloadEchoChan(char *measout, int ChanId, int EchoId, int nReps,
                    int nSlices, int nEchoes, int nRows, int nChans,
                    int nCols, int nPCNs, int FasterDim,
                    MRI *mriReal, MRI *mriImag);
int MDHdumpADC(char *measoutpath, char *outfile, int binary);

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
char *srcdir, *outdir, *adcfile=NULL;
int rev, infoonly=1, dumpmdh=0, debug=0, adcstats=0;

int DimSpeced=0;
int nPerLine, nPCNs, nEchos, nPELs, nSlices, nFrames, FastestDim=0;
float Thickness, DistFact=0.0, dcSag, dcCor, dcTra, TR=0.0, TE[500];
float PhaseEncodeFOV, ReadoutFOV, FlipAngle;
int Strict = 1;
int nthpcn;
int BinaryADCDump=0;
char *the_basename = "meas";
int DumpPCN = 1;

/*------------------------------------------------------------------*/
int main(int argc, char **argv) {
  FILE *fp;
  long offset;
  float adc[10000];
  //MDH *mdh = NULL;
  int n,d, err;
  char *tmpstr;
  MRI *pcnr=0, *pcni=0;
  MRI **echor, **echoi;
  float *rptr, *iptr;
  //int s,f,e,p,l;
  char fname[1000];
  char measoutpath[2000];
  int mdhversion;
  char *mdhascfile;
  long nHit, nHitTot, nHitPCN, nHitTotExp, nHit10pct;
  MDH *mdh = NULL;
  int nargs;

  ConvertMDH_VB13(argv[1],argv[2]);

  //printf("MDH13 %d\n",sizeof(MDH_VB13));
  //DumpMDH_VB13(argv[1]);
  return(0);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_convert_mdh.c,v 1.27 2011/03/02 00:04:14 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  sprintf(measoutpath,"%s/%s.out",srcdir,the_basename);
  if (!fio_FileExistsReadable(measoutpath)) {
    printf("ERROR: %s does not exist or is not readable\n",measoutpath);
    exit(1);
  }

  mdhversion = MDHversion(measoutpath);
  if (mdhversion == -1) {
    printf("ERROR: cannot find MDH ascii file\n");
    exit(1);
  }
  if (mdhversion != 15 && mdhversion != 21) {
    printf("ERROR: MDH version %d unrecognized\n",mdhversion);
    exit(1);
  }

  printf("MDH Version %d\n",mdhversion);
  mdhascfile = MDHascPath(srcdir);

  if (dumpmdh) {
    MDHdump(measoutpath);
    exit(0);
  }

  if (adcfile) {
    MDHdumpADC(measoutpath,adcfile,BinaryADCDump);
    exit(0);
  }

  if (adcstats) {
    MDHadcStats(measoutpath);
    exit(0);
  }

  if (DimSpeced==0) {
    nFrames  = MDHnFrames(measoutpath);
    nSlices  = MDHnSlices(measoutpath);
    nEchos   = MDHnEchos(measoutpath);
    nPELs    = MDHnPhaseEncodeLines(measoutpath);
    nPerLine = MDHnSamplesPerLine(measoutpath);
    nPCNs    = MDHnPCNs(measoutpath);
  }
  //if(FastestDim==0) FastestDim = MDHfastestDim(measoutpath);

  //if(TR==0) TR = MDHreadTR(measoutpath);
  //if(TE[0] == 0) for(n=0; n<nEchos; n++)  TE[n] = MDHreadTE(measoutpath,n);

  //Thickness  = MDHsliceDirCos(measoutpath,&dcSag,&dcCor,&dcTra);
  tmpstr = MDHparseMrProt(mdhascfile, "sSliceArray.asSlice[0].dThickness");
  sscanf(tmpstr,"%f",&Thickness);
  free(tmpstr);
  tmpstr = MDHparseMrProt(mdhascfile, "sGroupArray.asGroup[0].dDistFact");
  if (tmpstr != NULL) {
    sscanf(tmpstr,"%f",&DistFact);
    free(tmpstr);
  }
  Thickness = Thickness * (1.0 + DistFact);
  tmpstr = MDHparseMrProt(mdhascfile, "sSliceArray.asSlice[0].dPhaseFOV");
  sscanf(tmpstr,"%f",&PhaseEncodeFOV);
  free(tmpstr);
  tmpstr = MDHparseMrProt(mdhascfile, "sSliceArray.asSlice[0].dReadoutFOV");
  sscanf(tmpstr,"%f",&ReadoutFOV);
  free(tmpstr);
  tmpstr = MDHparseMrProt(mdhascfile, "dFlipAngleDegree");
  if (tmpstr != NULL) {
    sscanf(tmpstr,"%f",&FlipAngle);
    free(tmpstr);
  }

  printf("TR = %g, nSPL = %d, nPCNs = %d\n",TR,nPerLine,nPCNs);
  printf("nFrames = %d, nSlices = %d, nEchos = %d, nPELs = %d\n",
         nFrames,nSlices,nEchos,nPELs);
  printf("Thick = %g, DC: %g, %g, %g\n",Thickness,dcSag,dcCor,dcTra);
  printf("Fastest Dim: %d\n",FastestDim);
  printf("FOV:  Phase = %g, Readout = %g\n",PhaseEncodeFOV,ReadoutFOV);
  //for(n=0; n<nEchos; n++)  printf("%2d  %8.4f\n",n,TE[n]);

  if (infoonly || dumpmdh) exit(0);

  /* Create the output directory */
  err = mkdir(outdir,0777);
  if (err != 0 && errno != EEXIST) {
    printf("ERROR: creating %s\n",outdir);
    perror(NULL);
    exit(1);
  }

#if 0
  mriReal = MRIallocSequence(nPerLine,nPELs,nSlices,MRI_FLOAT,nFrames);
  mriImag = MRIallocSequence(nPerLine,nPELs,nSlices,MRI_FLOAT,nFrames);
  for (n=0; n<nEchos; n++) {
    printf("Echo %d\n",n);
    MDHloadEchoChan(measoutpath, 0, n, nFrames, nSlices, nEchos,
                    nPELs, nChans, nPerLine, nPCNs, FastestDim,
                    mriReal,mriImag);
    printf(" Saving real\n");
    sprintf(fname,"%s/echo%03dr.mgh",outdir,n+1);
    MRIwrite(mriReal,fname);
    printf(" Saving imag\n");
    sprintf(fname,"%s/echo%03di.mgh",outdir,n+1);
    MRIwrite(mriImag,fname);
  }
  MRIfree(&mriReal);
  MRIfree(&mriImag);
#endif


  //if(rev) printf("INFO: applying reversals\n");

  echor = (MRI**) calloc(sizeof(MRI*),nEchos);
  echoi = (MRI**) calloc(sizeof(MRI*),nEchos);
  printf("Allocating\n");
  for (n=0; n<nEchos; n++) {
    echor[n] = MRIallocSequence(nPerLine,nPELs,nSlices,MRI_FLOAT,nFrames);

    echor[n]->xsize = ReadoutFOV/(nPerLine/2);
    echor[n]->ysize = PhaseEncodeFOV/nPELs;
    echor[n]->zsize = Thickness;
    echor[n]->tr    = TR;
    echor[n]->te    = TE[n];
    echor[n]->flip_angle  = FlipAngle*3.1415/180;

    /* DC is correct for z, bogus for x and y */
    /* Need to change so DC is RAS */
    AssignFakeDCs(echor[n], dcSag, dcCor, dcTra);

    echoi[n] = MRIallocSequence(nPerLine,nPELs,nSlices,MRI_FLOAT,nFrames);
    echoi[n]->xsize = ReadoutFOV/(nPerLine/2);
    echoi[n]->ysize = PhaseEncodeFOV/nPELs;
    echoi[n]->zsize = Thickness;
    echoi[n]->tr    = TR;
    echoi[n]->te    = TE[n];
    echoi[n]->flip_angle  = FlipAngle*3.1415/180;
    AssignFakeDCs(echoi[n], dcSag, dcCor, dcTra); /* bogus */

  }
  if (nPCNs > 0) {
    pcnr = MRIallocSequence(nPerLine,nPCNs,nSlices,MRI_FLOAT,nFrames);
    pcni = MRIallocSequence(nPerLine,nPCNs,nSlices,MRI_FLOAT,nFrames);
  }

  /*---------------------------------------------------------------*/
  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  printf("offset = %ld\n",offset);
  fseek(fp,offset,SEEK_SET);

  /*---------------------------------------------------------------*/
  nHitTotExp = (long)nSlices*nFrames*((long)nEchos*nPELs+nPCNs);
  printf("Loading  (lines to load = %ld)\n",nHitTotExp);
  nHitPCN = 0;
  nthpcn = 0;
  nHit = 0;
  nHitTot = 0;
  nHit10pct = nHitTotExp/10;
  while (1) {
    if (nHitTot % nHit10pct == 0) printf("%6.1f%%\n",100*(float)nHitTot/nHitTotExp);

    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    if (feof(fp)) break;
    if (debug) PrintMiniHeader(stdout,mdh);
    fread(adc,sizeof(float), 2*nPerLine, fp);
    if (feof(fp)) {
      //printf("WARNING: hit eof during data read\n");
      break;
    }
    if (mdh->Rep >= nFrames) {
      printf("ERROR: frame = %d >= nFrames = %d\n",mdh->Rep,nFrames);
      PrintMiniHeader(stdout,mdh);
      exit(1);
    }
    if (mdh->Slice >= nSlices) {
      printf("ERROR: slice = %d >= nSlices = %d\n",mdh->Slice,nSlices);
      exit(1);
    }
    d = mdh->Slice + nSlices * mdh->Rep;
    if (mdh->IsPCN) {
      /* Note: Echo becomes Line for PCNs */
      //rptr = (float*) pcnr->slices[d][mdh->Echo];
      //iptr = (float*) pcni->slices[d][mdh->Echo];
      rptr = (float*) pcnr->slices[d][nthpcn];
      iptr = (float*) pcni->slices[d][nthpcn];
      nthpcn ++;
      if (nthpcn == nPCNs) nthpcn = 0;
      nHitPCN ++;
    } else {
      rptr = (float*) echor[mdh->Echo]->slices[d][mdh->LoopCounterLine];
      iptr = (float*) echoi[mdh->Echo]->slices[d][mdh->LoopCounterLine];
      nHit ++;
    }
    /* copy the data */
    for (n=0; n < 2*nPerLine; n += 2) {
      (*rptr++) = adc[n];
      (*iptr++) = adc[n+1];
    }
    nHitTot ++;
  }
  printf("Loaded %ld lines, nHitPCN = %ld, nHit = %ld\n",nHitTot,nHitPCN,nHit);
  if (nHitTotExp != nHitTot) {
    printf("\nWARNING: number of lines loaded does "
           "not equal the number expected\n\n");
  }

  /* copy the meas.asc to the output directory */
  sprintf(fname,"cp %s/%s.asc %s/%s.asc",srcdir,the_basename,outdir,the_basename);
  system(fname);


#if 0
  printf("Loading\n");
  d = 0;
  for (f=0; f < nFrames; f++) {
    if (f==0 || f%10 ==0) printf("f=%d\n",f);
    for (s=0; s < nSlices; s++) {

      /* Phase Encode Navigators */
      for (p=0; p < nPCNs; p++) {
        fseek(fp,MDH_SIZE,SEEK_CUR); /* skip header */
        fread(adc,sizeof(float), 2*nPerLine, fp);
        rptr = (float*) pcnr->slices[d][p];
        iptr = (float*) pcni->slices[d][p];
        if ((p+1)%2 != 0 || !rev) {
          /* odd p */
          for (n=0; n < 2*nPerLine; n += 2) {
            (*rptr++) = adc[n];
            (*iptr++) = adc[n+1];
          }
        } else {
          /* even p (reverse) */
          for (n = 2*nPerLine-2; n >= 0; n -= 2) {
            (*rptr++) = adc[n];
            (*iptr++) = adc[n+1];
          }
        }
      }

      if (FastestDim == 1) {
        /* For echo as fastest dim */
        for (l=0; l < nPELs; l++) {
          for (e=0; e < nEchos; e++) {
            fseek(fp,MDH_SIZE,SEEK_CUR); /* skip header */
            fread(adc,sizeof(float), 2*nPerLine, fp);
            rptr = (float*) echor[e]->slices[d][l];
            iptr = (float*) echoi[e]->slices[d][l];
            if ((l+1)%2 != 0 || !rev) {
              /* odd line */
              for (n=0; n < 2*nPerLine; n += 2) {
                (*rptr++) = adc[n];
                (*iptr++) = adc[n+1];
              }
            } else {
              /* even l (reverse) */
              for (n = 2*nPerLine-2; n >= 0; n -= 2) {
                (*rptr++) = adc[n];
                (*iptr++) = adc[n+1];
              }
            }
          } /* phase encode line */
        } /* echo */
      }
      if (FastestDim == 2) {
        /* For phase encode line as fastest dim */
        /* Need to reverse lines for odd echos? */
        for (e=0; e < nEchos; e++) {
          for (l=0; l < nPELs; l++) {
            fseek(fp,MDH_SIZE,SEEK_CUR); /* skip header */
            fread(adc,sizeof(float), 2*nPerLine, fp);
            rptr = (float*) echor[e]->slices[d][l];
            iptr = (float*) echoi[e]->slices[d][l];
            if ((l+1)%2 != 0 || !rev) {
              /* odd line */
              for (n=0; n < 2*nPerLine; n += 2) {
                (*rptr++) = adc[n];
                (*iptr++) = adc[n+1];
              }
            } else {
              /* even l (reverse) */
              for (n = 2*nPerLine-2; n >= 0; n -= 2) {
                (*rptr++) = adc[n];
                (*iptr++) = adc[n+1];
              }
            }
          } /* phase encode line */
        } /* echo */
      }

      d++;
    } /* slice */
  } /* frame */
#endif

  fclose(fp);

  printf("Saving\n");
  for (n=0; n<nEchos; n++) {
    printf("Echo %d\n",n);
    sprintf(fname,"%s/echo%03dr.mgh",outdir,n+1);
    MRIwrite(echor[n],fname);
    MRIfree(&echor[n]);
    sprintf(fname,"%s/echo%03di.mgh",outdir,n+1);
    MRIwrite(echoi[n],fname);
    MRIfree(&echoi[n]);
  }
  if (nPCNs > 0) {
    sprintf(fname,"%s/pcnr.mgh",outdir);
    MRIwrite(pcnr,fname);
    sprintf(fname,"%s/pcni.mgh",outdir);
    MRIwrite(pcni,fname);
  }

  printf("Done\n");

  return(0);
  /*-----------------------------------------------*/

  return(0);
}
/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused, n;
  char **pargv, *option ;

  isflag(" "); /* shuts up compiler */

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--info"))    infoonly = 1;
    else if (!strcasecmp(option, "--dump"))    dumpmdh = 1;
    else if (!strcasecmp(option, "--binary"))  BinaryADCDump = 1;
    else if (!strcasecmp(option, "--adcstats")) adcstats = 1;
    else if (!strcasecmp(option, "--nopcn")) DumpPCN = 0;
    //else if (!strcasecmp(option, "--rev"))     rev = 1;

    else if (stringmatch(option, "--srcdir")) {
      if (nargc < 1) argnerr(option,1);
      srcdir = pargv[0];
      nargsused = 1;
    } else if (stringmatch(option, "--base")) {
      if (nargc < 1) argnerr(option,1);
      the_basename = pargv[0];
      nargsused = 1;
    } else if (stringmatch(option, "--outdir")) {
      if (nargc < 1) argnerr(option,1);
      outdir = pargv[0];
      infoonly = 0;
      nargsused = 1;
    } else if (stringmatch(option, "--dumpadc")) {
      if (nargc < 1) argnerr(option,1);
      adcfile = pargv[0];
      nargsused = 1;
    } else if (stringmatch(option, "--dim")) {
      if (nargc < 6) argnerr(option,6);
      sscanf(pargv[0],"%d",&nFrames);
      sscanf(pargv[1],"%d",&nSlices);
      sscanf(pargv[2],"%d",&nEchos);
      sscanf(pargv[3],"%d",&nPELs);
      sscanf(pargv[4],"%d",&nPerLine);
      sscanf(pargv[5],"%d",&nPCNs);
      DimSpeced = 1;
      nargsused = 6;
    } else if (stringmatch(option, "--fasterdim")) {
      if (nargc < 1) argnerr(option,1);
      if (strcmp(pargv[0],"echo")==0) FastestDim = 1;
      else if (strcmp(pargv[0],"line")==0) FastestDim = 2;
      else {
        printf("ERROR: fasterdim name %s unrecoginzied\n",pargv[0]);
        exit(1);
      }
      nargsused = 1;
    } else if (stringmatch(option, "--TR")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&TR);
      nargsused = 1;
    } else if (stringmatch(option, "--TE")) {
      if (nEchos == 0) {
        printf("ERROR: must spec nechos before --TE\n");
        exit(1);
      }
      if (nargc < nEchos) argnerr(option,nEchos);
      for (n=0; n<nEchos; n++) sscanf(pargv[n],"%f",&TE[n]);
      nargsused = nEchos;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("\n");
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --srcdir dir : directory with meas.out\n");
  printf("   --base basename : default is meas\n");
  printf("   --outdir dir : output directory\n");
  //printf("   --rev        : apply reversals \n");
  printf("   --dim nframes nslices nechoes nlines nperline npcns\n");
  printf("   --fasterdim dimname : line or echo\n");
  printf("   --TR TR \n");
  printf("   --TE TE1 TE2 ... \n");
  printf("\n");
  printf("   --dump : dump info about each adc to term (but not adc) \n");
  printf("   --dumpadc file : dump adc into text file (unless --binary)\n");
  printf("   --binary : dump adc into binary file\n");
  printf("   --nopcn : do not dump pcns\n");
  printf("   --adcstats : dump min, max, and avg over all adcs \n");
  printf("\n");
  printf("   --help : a short story about a Spanish boy named Manual\n");
  printf("\n");
  //printf("   --info       : no ouput, just dump info (default with no outdir)\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
  printf(

    "Probes or converts a Siemens proprietary format (MDH) file to MGH format. \n"
    "\n"
    "If an output directory is not specified, then the MDH file is probed,\n"
    "and the various parameters are printed to stdout.\n"
    "\n"
    "When an output is specified, the output will have two parts: phase\n"
    "correction navigators and echos.  The phase correction navigators will\n"
    "be stored in pcnr.mgh (real) and pcni.mgh (imaginary).  The echos will\n"
    "be stored in echoNr.mgh (real) and echoNi.mgh (imaginary), where N is\n"
    "the echo number (starting with 0).\n"
    "\n"
    "--srcdir dir\n"
    "\n"
    "Directory with the meas.out file and the ascii file. The meas.out file\n"
    "contains the measured k-space data in Siemens MDH format (which includes\n"
    "miniheaders). The meas.out file must be accompanied by an ascii file\n"
    "called either MrProt.asc or mrprot.asc (NUMARIS 4 VA15) or meas.asc\n"
    "(NUMARIS 4 VA21). The version is determined by the name of this file.\n"
    "\n"
    "--outdir dir\n"
    "\n"
    "Directory where the ouput will be stored. If this directory does not\n"
    "exist, it will be created (only for one level). \n"
    "\n"
    "--rev (DOES NOT WORK)\n"
    "\n"
    "DOES NOT WORK Reverse the readouts for even lines. Bug: should reverse \n"
    "lines for even number echos, but does not do this yet.\n"
    "\n"
    "--dump\n"
    "\n"
    "Print out the header info for each readout line (but not the ADC\n"
    "data itself). Only the source dir needs to be given for this option.\n"
    "\n"
    "--dumpadc file\n"
    "\n"
    "Print out ADC into file in simple ascii/text format. This is just the\n"
    "raw data, with no indication of what is what.\n"
    "\n"
    "-adcstats\n"
    "\n"
    "Print out min, max, and avg over all ADCs. Only the source dir \n"
    "needs to be given for this option. This is good for finding spikes.\n"
    "\n"
    "EXAMPLES:\n"
    "\n"
    "\nmri_convert_mdh --srcdir . --base spiral_32x32_1meas \n"
    "        --binary --dumpadc tmp.dump\n"
    "\n"
    "NOTES:\n"
    "\n"
    "Ncols is the number of samples where each sample has a real and\n"
    "complex component.\n"
    "\n"
    "\n"
    "BUGS:\n"
    "\n"
    "Lines of even-numbered echos are not reversed with --rev.\n"
    "\n"
    "Does not handle version V21 yet.\n"
    "\n"
    "Does not handle 3D sequences yet.\n"
    "\n"
    "Does not handle FLASH  sequences yet.\n"
    "\n"
    "Does not handle multiple channels yet.\n"
    "\n"
  );
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1) fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else     fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void) {
  if (srcdir == NULL) {
    printf("ERROR: no source directory\n");
    exit(1);
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  return;
}
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int isflag(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*------------------------------------------------------------*/
static int stringmatch(char *str1, char *str2) {
  if (! strcmp(str1,str2)) return(1);
  return(0);
}















/*----------------------------------------------------------*/
/*----------------------------------------------------------*/
/*----------------------------------------------------------
  MDHversion() - determines the NUMARIS version number given the path
  to the meas.out file. The version is determined by the name of the
  ascii
  file.
  ----------------------------------------------------------*/
int MDHversion(char *measoutpath) {
  char *measoutdir;
  char fname[2000];

  measoutdir = fio_dirname(measoutpath);

  sprintf(fname,"%s/MrProt.asc",measoutdir);
  if (fio_FileExistsReadable(fname))  return(15);

  sprintf(fname,"%s/mrprot.asc",measoutdir);
  if (fio_FileExistsReadable(fname))  return(15);

  sprintf(fname,"%s/%s.asc",measoutdir,the_basename);
  if (fio_FileExistsReadable(fname))  return(21);

  free(measoutdir);

  return(-1);
}
/*----------------------------------------------------------
  MDHascPath() - returns the path to the MDH ascii file. This
  will be either measoutdir/mrprot.asc (for NUMARIS 1.5) or
  measoutdir/MrProt.asc (also 1.5) or measoutdir/meas.asc
  (for NUMARIS 2.1). Returns NULL if this file cannot be
  found.
  ----------------------------------------------------------*/
char *MDHascPath(char *measoutdir) {
  char fname[2000];
  char *mrprot;
  int len;

  sprintf(fname,"%s/MrProt.asc",measoutdir);
  if (fio_FileExistsReadable(fname)) {
    len = strlen(fname);
    mrprot = (char *) calloc(sizeof(char),len+1);
    memmove(mrprot,fname,len);
    return(mrprot);
  }

  sprintf(fname,"%s/mrprot.asc",measoutdir);
  if (fio_FileExistsReadable(fname)) {
    len = strlen(fname);
    mrprot = (char *) calloc(sizeof(char),len+1);
    memmove(mrprot,fname,len);
    return(mrprot);
  }

  sprintf(fname,"%s/%s.asc",measoutdir,the_basename);
  if (fio_FileExistsReadable(fname)) {
    len = strlen(fname);
    mrprot = (char *) calloc(sizeof(char),len+1);
    memmove(mrprot,fname,len);
    return(mrprot);
  }

  return(NULL);
}

/*----------------------------------------------------------*/
MDH *ReadMiniHeader(FILE *fp, MDH *mdh, int mdhversion) {
  if (mdhversion == 15) return(ReadMiniHeader15(fp,mdh));
  if (mdhversion == 21) return(ReadMiniHeader21(fp,mdh));
  return(NULL);
}
/*----------------------------------------------------------*/
MDH *ReadMiniHeader15(FILE *fp, MDH *mdh) {
  unsigned long ultmp;
  unsigned short ustmp;
  float ftmp;

  if (mdh == NULL)  mdh = (MDH *) calloc(sizeof(MDH),1);

  fread(&ultmp,sizeof(long),1, fp); //DMALength
  fread(&ultmp,sizeof(long),1, fp); //MeasUID (long)
  fread(&mdh->ScanCounter,sizeof(long),1, fp); //ScanCounter
  fread(&mdh->ulTimeStamp,sizeof(long),1, fp);
  mdh->TimeStamp = 2.5*mdh->ulTimeStamp; // dont know why its 2.5
  fread(&ultmp,sizeof(long),1, fp); //PMUTimeStamp
  fread(&mdh->BitMask1,sizeof(long),1, fp);
  //fread(&ultmp,sizeof(long),1, fp); //BitMask2 (VA21)

  fread(&mdh->Ncols,sizeof(short),1, fp);          //SamplesInScan
  fread(&mdh->UsedChannels,sizeof(short),1, fp);   //Number of used channels
  fread(&mdh->LoopCounterLine,sizeof(short),1, fp);//LoopCounter.Line
  fread(&ustmp,sizeof(short),1, fp);          //Acquistion
  fread(&mdh->Slice,sizeof(short),1, fp);     //Slice number
  fread(&mdh->Partition,sizeof(short),1, fp); //Partition number
  fread(&mdh->Echo,sizeof(short),1, fp);      //Echo number
  fread(&ustmp,sizeof(short),1, fp);          //Phase
  fread(&mdh->Rep,sizeof(short),1, fp);       //Repetition number (ie,frame)
  fread(&ustmp,sizeof(short),1, fp);          //Set ??
  fread(&ustmp,sizeof(short),1, fp);          //Seg ??
  fread(&ustmp,sizeof(short),1, fp); //LoopCounterFree (VA15)
  fread(&mdh->CutOffDataPre, sizeof(short),1, fp); //CutOffData.Pre
  fread(&mdh->CutOffDataPost,sizeof(short),1, fp); //CutOffData.Post
  fread(&mdh->KSpaceCenterCol,sizeof(short),1, fp); //KSpaceCenterColumn
  fread(&ustmp,sizeof(short),1, fp); //Dummy

  fread(&mdh->ReadoutOffCenter,sizeof(float),1, fp); //ReadOutOffCenter
  fread(&ultmp,sizeof(long),1, fp); //TimeSinceLastRF
  mdh->PED = ultmp * 2.5; // convert to ms

  fread(&mdh->KSpaceCenterLine,sizeof(short),1, fp); //KSpaceCentreLineNo
  fread(&ustmp,sizeof(short),1, fp); //KSpaceCentrePartitionNo
  fseek(fp,14*sizeof(short),SEEK_CUR); // FreePara (VA15)

  fread(&mdh->SlicePosSag,sizeof(float),1, fp);
  fread(&mdh->SlicePosCor,sizeof(float),1, fp);
  fread(&mdh->SlicePosTra,sizeof(float),1, fp);

  fread(&ftmp,sizeof(float),1, fp); //Quaternion1
  fread(&ftmp,sizeof(float),1, fp); //Quaternion2
  fread(&ftmp,sizeof(float),1, fp); //Quaternion3
  fread(&ftmp,sizeof(float),1, fp); //Quaternion4

  fread(&mdh->ChannelId,sizeof(long),1, fp); //ChannelId

  if ( (mdh->BitMask1 & MDH_BM_PHASECOR) != 0) mdh->IsPCN = 1;
  else  mdh->IsPCN = 0;

  return(mdh);
}
/*----------------------------------------------------------*/
MDH *ReadMiniHeader21(FILE *fp, MDH *mdh) {
  unsigned long ultmp;
  unsigned short ustmp;
  float ftmp;

  if (mdh == NULL)  mdh = (MDH *) calloc(sizeof(MDH),1);

  fread(&ultmp,sizeof(long),1, fp); //DMALength
  fread(&ultmp,sizeof(long),1, fp); //MeasUID (long)
  fread(&mdh->ScanCounter,sizeof(long),1, fp); //ScanCounter
  fread(&mdh->ulTimeStamp,sizeof(long),1, fp);
  mdh->TimeStamp = 2.5*mdh->ulTimeStamp; // dont know why its 2.5
  fread(&ultmp,sizeof(long),1, fp); //PMUTimeStamp
  fread(&mdh->BitMask1,sizeof(long),1, fp);
  fread(&ultmp,sizeof(long),1, fp); //BitMask2 (VA21)

  fread(&mdh->Ncols,sizeof(short),1, fp);          //SamplesInScan
  fread(&mdh->UsedChannels,sizeof(short),1, fp);   //Number of used channels
  fread(&mdh->LoopCounterLine,sizeof(short),1, fp);//LoopCounter.Line
  fread(&ustmp,sizeof(short),1, fp);               //Acquistion
  fread(&mdh->Slice,sizeof(short),1, fp);          //Slice number
  fread(&mdh->Partition,sizeof(short),1, fp);      //Partition number
  fread(&mdh->Echo,sizeof(short),1, fp);           //Echo number
  fread(&ustmp,sizeof(short),1, fp);               //Phase
  fread(&mdh->Rep,sizeof(short),1, fp);            //Repetition number (ie,frame)
  fread(&ustmp,sizeof(short),1, fp);               //Set ??
  fread(&ustmp,sizeof(short),1, fp);               //Seg ??
  fread(&ustmp,sizeof(short),1, fp); //Ida (VA21)
  fread(&ustmp,sizeof(short),1, fp); //Idb (VA21)
  fread(&ustmp,sizeof(short),1, fp); //Idc (VA21)
  fread(&ustmp,sizeof(short),1, fp); //Idd (VA21)
  fread(&ustmp,sizeof(short),1, fp); //Ide (VA21)

  fread(&ustmp,sizeof(short),1, fp); //CutOffData.Pre
  fread(&ustmp,sizeof(short),1, fp); //CutOffData.Post
  fread(&ustmp,sizeof(short),1, fp); //KSpaceCenterColumn
  fread(&ustmp,sizeof(short),1, fp); //Dummy

  fread(&ftmp,sizeof(float),1, fp); //ReadOutOffCenter
  fread(&ultmp,sizeof(long),1, fp); //TimeSinceLastRF
  mdh->PED = ultmp * 2.5; // convert to ms

  fread(&ustmp,sizeof(short),1, fp); //KSpaceCentreLineNo
  fread(&ustmp,sizeof(short),1, fp); //KSpaceCentrePartitionNo
  fread(&ustmp,sizeof(short),1, fp); //IceProgramPara1 (VA21)
  fread(&ustmp,sizeof(short),1, fp); //IceProgramPara2 (VA21)
  fread(&ustmp,sizeof(short),1, fp); //IceProgramPara3 (VA21)
  fread(&ustmp,sizeof(short),1, fp); //IceProgramPara4 (VA21)
  fseek(fp,4*sizeof(short),SEEK_CUR); // FreePara (VA21)

  fread(&mdh->SlicePosSag,sizeof(float),1, fp);
  fread(&mdh->SlicePosCor,sizeof(float),1, fp);
  fread(&mdh->SlicePosTra,sizeof(float),1, fp);

  fread(&ftmp,sizeof(float),1, fp); //Quaternion1
  fread(&ftmp,sizeof(float),1, fp); //Quaternion2
  fread(&ftmp,sizeof(float),1, fp); //Quaternion3
  fread(&ftmp,sizeof(float),1, fp); //Quaternion4

  fread(&ultmp,sizeof(long),1, fp); //ChannelId

  if ( (mdh->BitMask1 & MDH_BM_PHASECOR) != 0) mdh->IsPCN = 1;
  else  mdh->IsPCN = 0;

  return(mdh);
}
/*---------------------------------------------------*/
int PrintMiniHeader(FILE *fp, MDH *mdh) {
  int n;
  fprintf(fp,"TimeStamp      %f, PED = %f\n",mdh->TimeStamp,mdh->PED);
  fprintf(fp,"IsPCN          %d\n",mdh->IsPCN);
  fprintf(fp,"ScanCounter    %ld\n",mdh->ScanCounter);
  fprintf(fp,"Line           %d\n",mdh->LoopCounterLine);
  fprintf(fp,"Echo           %d\n",mdh->Echo);
  fprintf(fp,"Slice          %d\n",mdh->Slice);
  fprintf(fp,"Partition      %d\n",mdh->Partition);
  fprintf(fp,"Rep            %d\n",mdh->Rep);
  fprintf(fp,"KSCenterCol    %d\n",mdh->KSpaceCenterCol);
  fprintf(fp,"KSCenterLine   %d\n",mdh->KSpaceCenterLine);
  fprintf(fp,"CutOffDataPre  %d\n",mdh->CutOffDataPre);
  fprintf(fp,"CutOffDataPost %d\n",mdh->CutOffDataPost);
  fprintf(fp,"ROOffCenter    %f\n",mdh->ReadoutOffCenter);
  fprintf(fp,"NoChannels     %d\n",mdh->UsedChannels);
  fprintf(fp,"ChannelId      %ld\n",mdh->ChannelId);
  fprintf(fp,"SP  %7.2f %7.2f %7.2f \n",
          mdh->SlicePosSag,mdh->SlicePosCor,mdh->SlicePosTra);
  fprintf(fp,"Ncols        %d\n",mdh->Ncols);
  fprintf(fp,"BitMask1     %lx\n",mdh->BitMask1);
  for (n = 0; n < 32; n++) printf("%d ",n%10);
  printf("\n");
  for (n = 0; n < 32; n++) printf("%ld ",(mdh->BitMask1 >> n)&1);
  printf("\n");

  return(0);
  //fprintf(fp,"Reflect      %ld\n",mdh->BitMask1 & MDH_BM_REFLECT);
  //fprintf(fp,"UsedChannels %d\n",mdh->UsedChannels);
}
/*---------------------------------------------------*/
float MDHtimeStamp0(char *measoutpath) {
  MDH *mdh=NULL;
  FILE *fp;
  float TimeStamp0;
  unsigned long offset;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  fp = fopen(measoutpath,"r");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n",measoutpath);
    return(-10000000);
  }
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);
  mdh = ReadMiniHeader(fp,mdh,mdhversion);
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
int MDHnSamplesPerLine(char *measoutpath) {
  MDH *mdh=NULL;
  FILE *fp;
  int nSamplesPerLine;
  unsigned long offset;
  int mdhversion;

  mdhversion = MDHversion(measoutpath);
  fp = fopen(measoutpath,"r");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n",measoutpath);
    return(-10000000);
  }
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);
  mdh = ReadMiniHeader(fp,mdh,mdhversion);
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
int MDHbytesPerLine(char *measoutpath) {
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
int MDHbytesPerChunk(char *measoutpath) {
  int BytesPerChunk;
  BytesPerChunk = MDHbytesPerLine(measoutpath) + MDH_SIZE;
  return(BytesPerChunk);
}
/*-----------------------------------------------------------------------
  MDHnPCNs() - returns the number of phase correction navigators
  collected per slice based on the number in the first slice.
  --------------------------------------------------------------------------*/
int MDHnPCNs(char *measoutpath) {
  MDH *mdh=NULL;
  FILE *fp;
  int nPCNs;
  unsigned long offset;
  int mdhversion;

  mdhversion = MDHversion(measoutpath);
  fp = fopen(measoutpath,"r");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n",measoutpath);
    return(-10000000);
  }
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  nPCNs = 0;
  while (1) {
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
    if (! mdh->IsPCN) break;
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
int MDHnSlices(char *measoutpath) {
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nSlices, nPCNs, nEchos, nPELs, chunksize;
  int BytesPerSlice, BytesPerSliceB;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  chunksize = MDHbytesPerChunk(measoutpath);
  if (chunksize < 0) return(-1);
  nPCNs   = MDHnPCNs(measoutpath);
  nEchos = MDHnEchos(measoutpath);
  nPELs  = MDHnPhaseEncodeLines(measoutpath);

  // Bytes per slice, including headers and PCNs
  BytesPerSlice  = chunksize*(nPCNs + nEchos*nPELs);
  // Bytes per slice, minus one header size
  BytesPerSliceB = BytesPerSlice-MDH_SIZE;

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  nSlices = 0;
  while (1) {
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    if (feof(fp)) break;
    if (mdh->Rep != 0) break;
    if (nSlices < mdh->Slice+1) nSlices = mdh->Slice+1;
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
int MDHslicePosition(char *measoutpath, int Slice, float *Sag, float *Cor, float *Tra) {
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nSlices, nPCNs, nEchos, nPELs, chunksize;
  int BytesPerSlice;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  chunksize = MDHbytesPerChunk(measoutpath);
  if (chunksize < 0) return(1);
  nPCNs   = MDHnPCNs(measoutpath);
  nEchos = MDHnEchos(measoutpath);
  nPELs  = MDHnPhaseEncodeLines(measoutpath);
  nSlices = MDHnSlices(measoutpath);

  if (Slice >= nSlices) {
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
  mdh = ReadMiniHeader(fp,mdh,mdhversion);

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
float MDHsliceThickness(char *measoutpath) {
  float dcSag, dcCor, dcTra;
  float Thickness;
  Thickness = MDHsliceDirCos(measoutpath,&dcSag,&dcCor,&dcTra);
  return(Thickness);
}
/*-----------------------------------------------------------------------
  MDHsliceDirCos - computes the direction cosines as pointing from the 0th
  slice to the 1st slice. Return value is the slice thickeness.
  --------------------------------------------------------------------------*/
float MDHsliceDirCos(char *measoutpath, float *dcSag, float *dcCor, float *dcTra) {
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
int MDHnFrames(char *measoutpath) {
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nFrames, nSlices, nPCNs, nEchos, nPELs, chunksize;
  int BytesPerFrame, BytesPerFrameB;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  chunksize = MDHbytesPerChunk(measoutpath);
  if (chunksize < 0) return(-1);
  nPCNs    = MDHnPCNs(measoutpath);
  nEchos  = MDHnEchos(measoutpath);
  nPELs   = MDHnPhaseEncodeLines(measoutpath);
  nSlices = MDHnSlices(measoutpath);

  // Bytes per frame, including headers and PCNs
  BytesPerFrame  = chunksize*nSlices*(nPCNs + nEchos*nPELs);
  // Bytes per frame, minus one header size
  BytesPerFrameB = BytesPerFrame-MDH_SIZE;

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  nFrames = 0;
  while (1) {
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    if (feof(fp)) break;
    if (nFrames < mdh->Rep+1) nFrames = mdh->Rep+1;
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
int MDHnEchos(char *measoutpath) {
  MDH *mdh=NULL;
  FILE *fp;
  int nEchos;
  unsigned long offset;
  int nPCNs, chunksize;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  chunksize = MDHbytesPerChunk(measoutpath);
  if (chunksize < 0) return(-1);
  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Get past the phase correction navigators
  nPCNs = MDHnPCNs(measoutpath);
  fseek(fp,nPCNs*chunksize,SEEK_CUR);

  nEchos = -1;
  while (1) {
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    if (feof(fp)) break;
    if (mdh->Slice != 0) break;
    if (nEchos < mdh->Echo+1) nEchos = mdh->Echo+1;
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
int MDHnPhaseEncodeLines(char *measoutpath) {
  MDH *mdh=NULL;
  FILE *fp;
  int nPhaseEncodeLines;
  unsigned long offset;
  int nPCNs, chunksize;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  chunksize = MDHbytesPerChunk(measoutpath);
  if (chunksize < 0) return(-1);
  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Get past the phase correction navigators
  nPCNs = MDHnPCNs(measoutpath);
  fseek(fp,nPCNs*chunksize,SEEK_CUR);

  nPhaseEncodeLines = -1;
  while (1) {
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    if (feof(fp)) break;
    if (mdh->Slice != 0) break;
    if (nPhaseEncodeLines < mdh->LoopCounterLine+1)
      nPhaseEncodeLines = mdh->LoopCounterLine+1;
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
float MDHreadTR(char *measoutpath) {
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nFrames, nSlices, nPCNs, nEchos, nPELs, chunksize;
  int BytesPerFrame, BytesPerFrameB;
  float TimeStamp0, TimeStamp1, TR;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  chunksize = MDHbytesPerChunk(measoutpath);
  if (chunksize < 0) return(-1);
  nPCNs    = MDHnPCNs(measoutpath);
  nEchos  = MDHnEchos(measoutpath);
  nPELs   = MDHnPhaseEncodeLines(measoutpath);
  nSlices = MDHnSlices(measoutpath);
  nFrames = MDHnFrames(measoutpath);
  if (nFrames == 1) return(0);

  // Bytes per frame, including headers and PCNs
  BytesPerFrame  = chunksize*nSlices*(nPCNs + nEchos*nPELs);
  // Bytes per frame, minus one header size
  BytesPerFrameB = BytesPerFrame-MDH_SIZE;

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Read time stamp from first frame/rep
  mdh = ReadMiniHeader(fp,mdh,mdhversion);
  TimeStamp0 = mdh->TimeStamp;

  // Jump to the beginning of the next frame/rep
  fseek(fp,BytesPerFrameB,SEEK_CUR);

  // Read time stamp from second frame/rep
  mdh = ReadMiniHeader(fp,mdh,mdhversion);
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
int MDHfastestDim(char *measoutpath) {
  MDH *mdh1, *mdh2;
  FILE *fp;
  unsigned long offset;
  int nPCNs, FastestDim, chunksize;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  chunksize = MDHbytesPerChunk(measoutpath);
  if (chunksize < 0) return(-1);
  nPCNs = MDHnPCNs(measoutpath);

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Get past PCNs
  fseek(fp,nPCNs*chunksize,SEEK_CUR);

  // Read the first non-PCN header
  mdh1 = ReadMiniHeader(fp,NULL,mdhversion);
  fseek(fp,2*mdh1->Ncols*sizeof(float),SEEK_CUR);

  // Read the second non-PCN header
  mdh2 = ReadMiniHeader(fp,NULL,mdhversion);

  if (mdh1->Echo != mdh2->Echo)  FastestDim = 1; // Echo is fastest
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
float MDHreadTE(char *measoutpath, int nthEcho) {
  MDH *mdh=NULL;
  FILE *fp;
  unsigned long offset;
  int nPCNs, nEchos, nPELs, chunksize, FastestDim;
  float TimeStamp0, TimeStamp1, TimeStamp1a, TimeStamp1b, TE;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  chunksize = MDHbytesPerChunk(measoutpath);
  if (chunksize < 0) return(-1);
  nPCNs    = MDHnPCNs(measoutpath);
  nEchos  = MDHnEchos(measoutpath);
  if (nthEcho >= nEchos) {
    printf("ERROR: requested TE=%d, max is %d\n",nthEcho,nEchos);
    return(-1);
  }
  nPELs = MDHnPhaseEncodeLines(measoutpath);
  FastestDim = MDHfastestDim(measoutpath);

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  // Get time stamp of first PCN as that closest to RF pulse
  mdh = ReadMiniHeader(fp,mdh,mdhversion);
  fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
  TimeStamp0 = mdh->TimeStamp;
  // Get past remaining PCNs
  fseek(fp,(nPCNs-1)*chunksize,SEEK_CUR);

  if (FastestDim == 1) {
    // Echo is the fastest (FLASH)
    // TE is time stamp at nthEcho. Is this correct?
    // This branch has never been tested!!!!!

    // Seek to the nthEcho
    fseek(fp,nthEcho*chunksize,SEEK_CUR);
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    TimeStamp1 = mdh->TimeStamp;

  } else {
    // Phase Encode is fastest (EPI)
    // TE is the average between the first and last PEL

    // Seek to first PEL for this echo
    fseek(fp,nthEcho*nPELs*chunksize,SEEK_CUR);
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
    TimeStamp1a = mdh->TimeStamp;

    // Seek to last PEL for this echo
    fseek(fp,(nPELs-2)*chunksize,SEEK_CUR);
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
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
char *MDHparseMrProt(char *file, char *TagString) {
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
  if (fp == NULL) {
    printf("ERROR: could not open dicom file %s\n",file);
    exit(1);
  }

  /* This section steps through the file char-by-char until
     the BeginStr is matched */
  dumpline = 0;
  nthchar = 0;
  while (1) {
    fseek(fp,nthchar, SEEK_SET);
    nTest = fread(TestStr,sizeof(char),LenBeginStr,fp);
    if (nTest != LenBeginStr) break;
    if (strcmp(TestStr,BeginStr)==0) {
      fseek(fp,nthchar, SEEK_SET);
      dumpline = 1;
      break;
    }
    nthchar ++;
  }
  free(TestStr);


  if (! dumpline) return(NULL); /* Could not match Begin String */

  /* Once the Begin String has been matched, this section
     searches each line until the TagString is matched
     or until the End String is matched */
  VariableValue = NULL;
  while (1) {
    rt = fgets(linestr,1000,fp);
    if (rt == NULL) break;

    if (strncmp(linestr,"### ASCCONV END ###",19)==0) break;

    sscanf(linestr,"%s %*s %*s",VariableName);

    if (strlen(VariableName) != strlen(TagString)) continue;

    if ( strcmp(VariableName,TagString)==0 ) {
      /* match found */
      sscanf(linestr,"%*s %*s %s",tmpstr2);
      VariableValue = (char *) calloc(strlen(tmpstr2)+1,sizeof(char));
      memmove(VariableValue, tmpstr2, strlen(tmpstr2));
      break;
    }
  }
  fclose(fp);

  //printf("Leaving SiemensAsciiTag() \n");fflush(stdout);fflush(stderr);
  fflush(stdout);
  fflush(stderr);

  return(VariableValue);
}
/*--------------------------------------------------------------
  MDHdump - dump the miniheader contents of the given meas.out
  ---------------------------------------------------------------*/
int MDHdump(char *measoutpath) {
  FILE *fp;
  unsigned long offset;
  int n;
  float TimeStamp0=0.0;
  MDH *mdh=NULL;
  int mdhversion;
  mdhversion = MDHversion(measoutpath);

  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  printf("offset = %ld\n",offset);
  fseek(fp,offset,SEEK_SET);

  n = 1;
  while (!feof(fp)) {
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    if (n==1) TimeStamp0 = mdh->TimeStamp;
    mdh->TimeStamp -= TimeStamp0;
    if (feof(fp)) break;
    printf("n = %d ---------------------------------\n",n);
    PrintMiniHeader(stdout,mdh);
    fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
    //fread(adc,sizeof(float), 2*mdh->Ncols, fp);
    n++;
  }

  fclose(fp);

  printf("n = %d\n",n);
  return(n);
}
/*--------------------------------------------------------------
  MDHdumpADC - dump the ADC contents of the given meas.out to outfile
  ---------------------------------------------------------------*/
int MDHdumpADC(char *measoutpath, char *outfile, int binary) {
  FILE *fp, *outfp;
  unsigned long offset;
  int m;
  MDH *mdh=NULL;
  int mdhversion,n;
  float adc[10000];

  mdhversion = MDHversion(measoutpath);
  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  outfp = fopen(outfile,"w");
  if (outfp == NULL) {
    printf("ERROR: opening %s\n",outfile);
    exit(1);
  }

  n = 0;
  while (!feof(fp)) {
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    if (feof(fp)) break;
    if (mdh->BitMask1 == 1) break;
    fread(adc,sizeof(float), 2*mdh->Ncols, fp);
    if (mdh->IsPCN && ! DumpPCN) continue;
    for (m=0; m < 2*mdh->Ncols; m++) {
      if (! binary)
        fprintf(outfp,"%20.40f\n",adc[m]);
      else
        fwrite(&adc[m],sizeof(float),1,outfp);
    }
    n++;
  }
  fclose(fp);
  fclose(outfp);

  printf("dumped %d lines\n",n);

  return(0);
}
/*--------------------------------------------------------------
  MDHadcStats - compute stats for data
  ---------------------------------------------------------------*/
int MDHadcStats(char *measoutpath) {
  FILE *fp;
  unsigned long offset;
  int n,m;
  //float TimeStamp0=0.0;
  MDH *mdh=NULL;
  int mdhversion;
  float adc[10000];
  double min, max, avg, sum;

  mdhversion = MDHversion(measoutpath);
  fp = fopen(measoutpath,"r");
  fread(&offset,sizeof(long),1, fp);
  //printf("offset = %ld\n",offset);
  fseek(fp,offset,SEEK_SET);

  min =  10e12;
  max = -10e12;
  sum = 0;
  n = 1;
  while (!feof(fp)) {
    mdh = ReadMiniHeader(fp,mdh,mdhversion);
    if (feof(fp)) break;
    //if(n==1) TimeStamp0 = mdh->TimeStamp;
    //mdh->TimeStamp -= TimeStamp0;
    //printf("n = %d ---------------------------------\n",n);
    //PrintMiniHeader(stdout,mdh);
    //fseek(fp,2*mdh->Ncols*sizeof(float),SEEK_CUR);
    fread(adc,sizeof(float), 2*mdh->Ncols, fp);
    for (m=0; m < 2*mdh->Ncols; m++) {
      if (min > adc[m]) min = adc[m];
      if (max < adc[m]) max = adc[m];
      sum = sum + adc[m];
    }
    n++;
  }
  fclose(fp);

  avg = sum/(2*n*mdh->Ncols);
  printf("n = %d, min = %g, max = %g, avg = %g\n",n,min,max,avg);
  return(n);
}


/*----------------------------------------------------------------
  AssignFakeDCs() - assigns direction cosines that are orthogonal.
  The z DC is correct, but the x and y are bogus.
  ---------------------------------------------------------------*/
int AssignFakeDCs(MRI *mri, float dcSag, float dcCor, float dcTra) {
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
/*-------------------------------------------------------------------
  MDHind2sub() - converts a zero-based index into a multidimensional
   array into a set of zero-based subscripts into the array.
  -------------------------------------------------------------------*/
int *MDHind2sub(int *siz, int ind, int ndim, int *sub) {
  int i,n,cps[100];

  if (sub == NULL) sub = (int *) calloc(sizeof(int),ndim);

  cps[0] = 1;
  for (n=0; n<ndim-1; n++) cps[n+1] = cps[n]*siz[n];

  //for(n=0; n<ndim; n++) printf("n=%d, cps[n]=%d\n",n,cps[n]);

  //printf("ind = %d\n",ind);
  ind--;
  for (i=ndim-1; i >=0;  i--) {
    sub[i] = floor((float)ind/cps[i]);
    ind = ind%cps[i];
  }

  //for(n=0; n<ndim; n++) printf("n=%d, siz[n]=%d, sub[n]=%d\n",n,siz[n],sub[n]);

  return(sub);
}

/*---------------------------------------------------------------
  MDHloadEchoChan() - loads the data for a given echo and channel
  without trying to figure anything out about the size of any
  of the dimensions. No reversals are applied.
  ----------------------------------------------------------------*/
int MDHloadEchoChan(char *measout, int ChanId, int EchoId, int nReps,
                    int nSlices, int nEchoes, int nRows, int nChans,
                    int nCols, int nPCNs, int FasterDim,
                    MRI *mriReal, MRI *mriImag) {
  extern int Strict;
  MDH *mdh;
  int mdhversion;
  FILE *fp;
  long offset;
  float *adc, *rptr, *iptr;
  int rep, slice, n1, n2, n1max=0, n2max=0, row=0, echo=0, chan, col, d;
  int ADCSizeBytes;

  fp = fopen(measout,"r");
  if (fp == NULL) {
    printf("ERROR: could not open %s\n",measout);
    exit(1);
  }

  printf("EchoId = %d, ChanId = %d\n",EchoId,ChanId);
  if (FasterDim == 1) { // Echo is faster
    printf("Echo is faster\n");
    n1max = nRows;
    n2max = nEchoes;
  }
  if (FasterDim == 2) { // Row is faster
    printf("Row is faster\n");
    n1max = nEchoes;
    n2max = nRows;
  }

  fread(&offset,sizeof(long),1, fp);
  fseek(fp,offset,SEEK_SET);

  ADCSizeBytes = 2*nCols * sizeof(float); // 2 = real/imag, 4 = float
  adc = (float *) malloc(ADCSizeBytes);
  mdh = NULL;
  mdhversion = MDHversion(measout);

  for (rep = 0; rep < nReps; rep++) {
    //printf("\n");
    //printf("rep = %d/%d: Slice:  ",rep+1,nReps);
    printf("%d ",rep+1);
    fflush(stdout);

    for (slice = 0; slice < nSlices; slice++) {
      //printf(" %d",slice);
      d = slice + nSlices * rep;

      mdh = ReadMiniHeader(fp,mdh,mdhversion);
      if (mdh->Slice != slice && Strict) {
        printf("ERROR: mdh->Slice = %d, slice = %d\n",
               mdh->Slice, slice);
        exit(1);
      }
      if (mdh->Rep != rep && Strict) {
        printf("ERROR: mdh->Rep = %d, rep = %d\n",
               mdh->Rep, rep);
        exit(1);
      }
      fseek(fp,ADCSizeBytes,SEEK_CUR);

      /* skip the PCNs */
      fseek(fp,(nPCNs-1)*(MDH_SIZE+ADCSizeBytes),SEEK_CUR);

      for (n1 = 0; n1 < n1max; n1++) {

        for (n2 = 0; n2 < n2max; n2++) {

          if (FasterDim == 1) { // Echo is faster
            row  = n1;
            echo = n2;
          }
          if (FasterDim == 2) { // Row is faster
            echo = n1;
            row  = n2;
          }

          for (chan = 0; chan < nChans; chan++) {

            //printf("rep=%d, slc=%d, n1=%d, n2=%d, chan = %d\n",
            //   rep,slice,n1,n2,chan);

            /* Skip the mdh header */
            //fseek(fp,MDH_SIZE,SEEK_CUR);
            mdh = ReadMiniHeader(fp,mdh,mdhversion);
            if (feof(fp)) {
              printf("ERROR: EOF\n");
              exit(1);
            }
            if (echo != mdh->Echo  && Strict) {
              printf("ERROR: mdh->Echo = %d, echo = %d \n",
                     mdh->Echo,echo);
              exit(1);
            }
            if (echo == 0 && mdh->LoopCounterLine != row && Strict) {
              printf("ERROR: mdh->LCL = %d, row = %d (echo = %d)\n",
                     mdh->LoopCounterLine,row,echo);
              exit(1);
            }
            if (echo == 1 && mdh->LoopCounterLine != 63-row && Strict) {
              printf("ERROR: mdh->LCL = %d, row = %d (echo = %d)\n",
                     mdh->LoopCounterLine,row,echo);
              exit(1);
            }

            /* Skip the adc if its not the one we want */
            if (chan != ChanId || echo != EchoId) {
              fseek(fp,ADCSizeBytes,SEEK_CUR);
              if (feof(fp)) {
                printf("ERROR: EOF\n");
                exit(1);
              }
              continue;
            }

            fread(adc,sizeof(float), 2*nCols, fp);
            if (feof(fp)) {
              printf("ERROR: EOF\n");
              exit(1);
            }

            rptr = (float*) mriReal->slices[d][row];
            iptr = (float*) mriImag->slices[d][row];
            for (col=0; col < 2*nCols; col += 2) {
              (*rptr++) = adc[col];
              (*iptr++) = adc[col+1];
            }


          } /* channel */

        } /* n2 (echo or row) */
      } /* n1 (row or echo) */
    } /* slice */
  } /* rep */
  printf("\n");

  fclose(fp);
  return(0);
}
int PrintMDH_VB13(FILE *fp, MDH_VB13 *mdh)
{
  fprintf(fp,"ulFlagsAndDMALength %d \n",mdh->ulFlagsAndDMALength);
  fprintf(fp,"lMeasUID            %d\n",mdh->lMeasUID);
  fprintf(fp,"ulScanCounter       %d\n",mdh->ulScanCounter);
  fprintf(fp,"ulTimeStamp         %d\n",mdh->ulTimeStamp);
  fprintf(fp,"ulPMUTimeStamp      %d\n",mdh->ulPMUTimeStamp);
  fprintf(fp,"aulEvalInfoMask1    %x\n",mdh->aulEvalInfoMask[0]);
  fprintf(fp,"aulEvalInfoMask2    %x\n",mdh->aulEvalInfoMask[1]);
  fprintf(fp,"ushSamplesInScan    %d\n",mdh->ushSamplesInScan);
  fprintf(fp,"ushUsedChannels     %d\n",mdh->ushUsedChannels);
  fprintf(fp,"ushLine             %d\n",mdh->ushLine);
  fprintf(fp,"ushAcquisition      %d\n",mdh->ushAcquisition);
  fprintf(fp,"ushSlice            %d\n",mdh->ushSlice);
  fprintf(fp,"ushPartition        %d\n",mdh->ushPartition);
  fprintf(fp,"ushEcho             %d\n",mdh->ushEcho);
  fprintf(fp,"ushPhase            %d\n",mdh->ushPhase);
  fprintf(fp,"ushRepetition       %d\n",mdh->ushRepetition);
  fprintf(fp,"ushSet              %d\n",mdh->ushSet);
  fprintf(fp,"ushSeg              %d\n",mdh->ushSeg);
  fprintf(fp,"ushIda              %d\n",mdh->ushIda);
  fprintf(fp,"ushIdb              %d\n",mdh->ushIdb);
  fprintf(fp,"ushIdc              %d\n",mdh->ushIdc);
  fprintf(fp,"ushIdd              %d\n",mdh->ushIdd);
  fprintf(fp,"ushIde              %d\n",mdh->ushIde);

  fprintf(fp,"ushPre              %d\n",mdh->ushPre);
  fprintf(fp,"ushPost             %d\n",mdh->ushPost);
  fprintf(fp,"ushKSpaceCentreColumn %d\n",mdh->ushKSpaceCentreColumn);
  fprintf(fp,"ushDummy            %d\n",mdh->ushDummy);

  fprintf(fp,"fReadOutOffcenter   %f\n",mdh->fReadOutOffcenter);
  fprintf(fp,"ulTimeSinceLastRF   %d\n",mdh->ulTimeSinceLastRF);

  fprintf(fp,"ushKSpaceCentreLineNo %d\n",mdh->ushKSpaceCentreLineNo);
  fprintf(fp,"ushKSpaceCentrePartitionNo %d\n",mdh->ushKSpaceCentrePartitionNo);
  fprintf(fp,"aushIceProgramPara1 %d\n",mdh->aushIceProgramPara[0]);
  fprintf(fp,"aushIceProgramPara2 %d\n",mdh->aushIceProgramPara[1]);
  fprintf(fp,"aushIceProgramPara3 %d\n",mdh->aushIceProgramPara[2]);
  fprintf(fp,"aushIceProgramPara4 %d\n",mdh->aushIceProgramPara[3]);
  fprintf(fp,"aushFreePara1       %d\n",mdh->aushFreePara[0]);
  fprintf(fp,"aushFreePara2       %d\n",mdh->aushFreePara[1]);
  fprintf(fp,"aushFreePara3       %d\n",mdh->aushFreePara[2]);
  fprintf(fp,"aushFreePara4       %d\n",mdh->aushFreePara[3]);

  fprintf(fp,"fSag                %f\n",mdh->fSag);
  fprintf(fp,"fCor                %f\n",mdh->fCor);
  fprintf(fp,"fTra                %f\n",mdh->fTra);
  fprintf(fp,"afQuaternion1       %f\n",mdh->afQuaternion[0]);
  fprintf(fp,"afQuaternion2       %f\n",mdh->afQuaternion[1]);
  fprintf(fp,"afQuaternion3       %f\n",mdh->afQuaternion[2]);
  fprintf(fp,"afQuaternion4       %f\n",mdh->afQuaternion[3]);

  fprintf(fp,"ushChannelId        %d\n",mdh->ushChannelId);
  fprintf(fp,"ushPTABPosNeg       %d\n",mdh->ushPTABPosNeg);

  return(0);

}

int DumpMDH_VB13(char *measfile)
{
  FILE *fp;
  int offset;
  MDH_VB13 mdh;
  int nread,sz,ADCSizeBytes,n, isPCN;
  long nth;

  sz = sizeof(MDH_VB13);

  fp = fopen(measfile,"r");
  if(fp == NULL){
    printf("ERROR: opening %s\n",measfile);
    return(1);
  }
  fread(&offset,sizeof(int),1, fp);
  printf("offset = %d\n",offset);
  fseek(fp,offset,SEEK_SET);

  nth = 0;
  while(1){
    nread = fread(&mdh,sz,1,fp);
    if(nread != 1) break;
    printf("%6ld -----------------------\n",nth);
    PrintMDH_VB13(stdout, &mdh);

    printf("BM  ");
    for (n = 0; n < 32; n++) printf("%d ",n%10);
    printf("\n");
    printf("BM1 ");
    for (n = 0; n < 32; n++) printf("%d ",(mdh.aulEvalInfoMask[0] >> n)&1);
    printf("\n");
    printf("BM2 ");
    for (n = 0; n < 32; n++) printf("%d ",(mdh.aulEvalInfoMask[1] >> n)&1);
    printf("\n");

    if(mdh.aulEvalInfoMask[0] & MDH_BM_PHASECOR){
      if(mdh.aulEvalInfoMask[0] & MDH_ONLINE) isPCN = 2;
      else isPCN = 1;
    }
    else isPCN = 0;
    printf("PCN %d\n",isPCN);

    ADCSizeBytes = mdh.ushSamplesInScan*2*sizeof(float); // 2 = real/imag, 4 = float
    fseek(fp,ADCSizeBytes,SEEK_CUR);
    nth = nth + 1;
  }
  fclose(fp);

  return(0);
}

int ConvertMDH_VB13(char *measfile, char *outbase)
{
  FILE *fp, *fpkd, *fpki, *fppd, *fppi, *fpi=NULL, *fpd=NULL;
  FILE *fppd2, *fppi2;
  unsigned int offset;
  MDH_VB13 mdh;
  int nread,sz,ADCSizeBytes, isPCN, err, isReflected;
  long nth;
  char *outdir;
  char tmpstr[1000];
  float ADC[1000];

  sz = sizeof(MDH_VB13);
  if(sz != 128){
    printf("ERROR: sizeof(MDH_VB13) = %d, != 128\n",(int)sizeof(MDH_VB13));
    return(1);
  }

  // Open the input meas.dat file
  fp = fopen(measfile,"r");
  if(fp == NULL){
    printf("ERROR: opening %s\n",measfile);
    return(1);
  }
  fread(&offset,sizeof(unsigned int),1, fp);
  printf("offset = %d\n",offset);
  fseek(fp,offset,SEEK_SET);

  // Create the output dir
  outdir = fio_dirname(outbase);
  err = mkdir(outdir,0777);
  if(err != 0 && errno != EEXIST) {
    printf("ERROR: creating directory %s\n",outdir);
    perror(NULL);
    return(1);
  }
  free(outdir);

  // Need four output files:
  sprintf(tmpstr,"%s-kdata.dat",outbase);
  fpkd = fopen(tmpstr,"w");
  if(fpkd == NULL){
    printf("ERROR: creating opening %s\n",tmpstr);
    perror(NULL);
    return(1);
  }
  sprintf(tmpstr,"%s-kdata.info",outbase);
  fpki = fopen(tmpstr,"w");
  if(fpki == NULL){
    printf("ERROR: creating opening %s\n",tmpstr);
    perror(NULL);
    return(1);
  }
  sprintf(tmpstr,"%s-kpcn.dat",outbase);
  fppd = fopen(tmpstr,"w");
  if(fppd == NULL){
    printf("ERROR: creating opening %s\n",tmpstr);
    perror(NULL);
    return(1);
  }
  sprintf(tmpstr,"%s-kpcn.info",outbase);
  fppi = fopen(tmpstr,"w");
  if(fppi == NULL){
    printf("ERROR: creating opening %s\n",tmpstr);
    perror(NULL);
    return(1);
  }

  sprintf(tmpstr,"%s-kpcn2.dat",outbase);
  fppd2 = fopen(tmpstr,"w");
  if(fppd2 == NULL){
    printf("ERROR: creating opening %s\n",tmpstr);
    perror(NULL);
    return(1);
  }
  sprintf(tmpstr,"%s-kpcn2.info",outbase);
  fppi2 = fopen(tmpstr,"w");
  if(fppi2 == NULL){
    printf("ERROR: creating opening %s\n",tmpstr);
    perror(NULL);
    return(1);
  }

  nth = -1;
  while(1){
    nth = nth + 1;

    nread = fread(&mdh,sz,1,fp);
    if(nread != 1) {
      printf("WARNING: MDH: unexpected end of file, nth = %ld\n",nth);
      break;
    }

    if(mdh.aulEvalInfoMask[0] & MDH_ACQEND) break;

    ADCSizeBytes = mdh.ushSamplesInScan*2*sizeof(float);
    nread = fread(ADC,ADCSizeBytes,1,fp);
    if(nread != 1) {
      printf("WARNING: DATA: unexpected end of file, nth = %ld\n",nth);
      break;
    }

    isPCN = 0;
    if(mdh.aulEvalInfoMask[0] & MDH_PHASECOR) isPCN = 1;
    if( !(mdh.aulEvalInfoMask[0] & MDH_PHASECOR) &&
	!(mdh.aulEvalInfoMask[0] & MDH_ONLINE)  ){
      isPCN = 2;
      // 2D PhaseCor
      //continue;
    }

    isReflected = 0;
    if(mdh.aulEvalInfoMask[0] & MDH_REFLECT) isReflected = 1;


    if(isPCN == 0) {
      fpi = fpki;
      fpd = fpkd;
    }
    if(isPCN == 1) {
      fpi = fppi;
      fpd = fppd;
    }
    if(isPCN == 2) {
      fpi = fppi2;
      fpd = fppd2;
    }

    // 1.Line 2.Chan 3.Set  4.Echo 5.Phase 6.Seg  7.Rep 8.Part 9.Slice  10.Time 11.isRefl
    fprintf(fpi,"%3d %3d %3d   %3d %3d %3d   %3d %3d %3d  %d  %6d\n",
	    mdh.ushLine, mdh.ushChannelId, 
	    mdh.ushSet, mdh.ushEcho, mdh.ushPhase, 
	    mdh.ushSeg,mdh.ushRepetition,
	    mdh.ushPartition, mdh.ushSlice,
	    isReflected, mdh.ulTimeStamp); 
    fwrite(ADC,ADCSizeBytes,1,fpd);
  }

  printf("INFO: found %ld MDHs\n",nth);

  fclose(fp);
  fclose(fpkd);
  fclose(fpki);
  fclose(fppd);
  fclose(fppi);
  fclose(fppd2);
  fclose(fppi2);

  return(0);
}
