/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


// mideface.h - include file for minimally invasive defacing tool

#ifndef MIDEFACE_H
#define MIDEFACE_H

#undef X

class MiDeface {
public:
  MRI *invol=NULL, *headmask=NULL, *xmask=NULL;
  MRIS *tempsurf=NULL;
  LABEL *templabel=NULL;
  MRI *outvol=NULL, *faceseg=NULL;
  MRIS *minsurf=NULL,*maxsurf=NULL;
  int nface1vox=0,nface2vox=0;
  double gmean1=0, gstddev1=0,gmean2=0, gstddev2=0;
  double min1=0, max1=0, min2=0, max2=0, mode2=0;
  double DistIn=0,DistOut=0,DeltaDist=-1,dL=-1;
  double DistInRaw=0,DistOutRaw=0;
  float *DistInList=NULL,*DistOutList=NULL;
  double DistInMin=2,DistInMax=100; // these will be set in the caller
  double DistOutMin=2,DistOutMax=100; 
  double DistInFrac=0.9, DistOutFrac=1.0;
  int FillType=1;
  double FillConstIn=0,FillConstOut=0;
  int cPad=5, rPad=5, sPad=5;
  int nxmask = 0;
  int DoRipple=0,rippleaxis=1;
  double ripplecenter[3] = {0,0,0},rippleperiod=20,rippleamp=1;
  int SegFace(void);
  int FaceIntensityStats(void);
  int SetDeltaDist(void);
  int DistanceBounds(void);
  int Deface(void);
  int VoxOutOfBounds(int c, int r, int s);
  MRI *Surf2VolProjFill(MRI *vol, double Dist, double FillVal);
  int PrintParams(FILE *fp);
  int PrintStats(FILE *fp);
  int ripple(LABEL *label);
  int watermark(LABEL *watermark, double dwatermark);
  const char *embedded_code = "mideface-freesurfer";
  const char *version = "v001";
  MRI *EmbedCode(MRI *input, MRI *output);
  int EmbedCodeCheck(const MRI *vol);
  int embeddebug=0;
};

#endif
