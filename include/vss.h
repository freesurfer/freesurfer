/**
 * @file  vss.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.3 $
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


/* $Id */
#ifndef VTXSAMPSPC_H
#define VTXSAMPSPC_H

#ifdef __cplusplus
extern "C"
{
#endif
  /* ----------------------------------------------------- */

#include <stdio.h>
#include "mrisurf.h"
#include "matrix.h"
#include "filecode.h"

#define VSS_NBRID_METHOD_DIST     1
#define VSS_NBRID_METHOD_KNEAREST 2
#define VSS_NBRID_METHOD_MLI      3
#define VSS_NBRID_METHOD_MRISURF  4

#define VSS_WEIGHT_METHOD_NBRAVG 1
#define VSS_WEIGHT_METHOD_MLI    2

#define VSS_ISAVG_METHOD_RANDEFF  1
#define VSS_ISAVG_METHOD_FIXEDEFF 2

#define VSS_USERDATAID_NONE     0
#define VSS_USERDATAID_KNEAREST 1
#define VSS_USERDATAID_SELAVG   2
#define VSS_USERDATAID_MLI      3

#define MVSS_CODE_STAT 1

  /* vertex-sampled space structure */
  typedef struct
  {
    int    version;
    char   sampname[1024];
    char   spacename[1024];
    int    nvertices;
    int    spacedim;
    float **vect;
    int   *dof;
    int   sznotes;
    char  *notes;
    int    udid;
    void  *ud;
  }
  VTX_SAMP_SPC, VSS;

  /* VSS representation of inter-space average */
  typedef struct
  {
    int    version;
    char   name[2048];
    int    nspaces;
    float *wspace;
    int    nstats;
    int   *dofstat;
    VSS   *avg;
    VSS   *var;
    char   model[2048];
    int    avgholes;
    int   sznotes;
    char  *notes;
    int    udid;
    void  *ud;
  }
  INTER_SPACE_STAT, ISS;

  /* vertex projection matrix */
  typedef struct
  {
    char vpmname[2048];
    int nrows, ncols;
    int affine;
    float **m;
  }
  VSS_VTXPROJMTX, VPM;

  /** vertex sampled stat (with covariance mtx) **/
  typedef struct
  {
    int   version;
    char  name[2048];
    int   avgholes;
    char  model[2048];
    int   dof;
    VSS  *avg;    /* Nv X Nh */
    VPM  *gcvm;   /* Nh X Nh, global (ie, inter-vertex) covariance matrix */
    VSS  *vtxvar; /* Nv X 1 */
    int   udid;
    void *ud;
    int   sznotes;
    char *notes;
  }
  VTX_SAMP_STAT, VSSTAT, VST;

  /* vertex neighborhood structure */
  typedef struct
  {
    int    nnbrs;
    int   *id;
    float *weight;
    float *dist;
  }
  VTXNBRHD;

  /* vertex interpolation table */
  typedef struct
  {
    int      version;
    char     interpspacename[1024];
    char     targsampname[1024];
    char     srcsampname[1024];
    int      nvertices;
    VTXNBRHD *vtxnbrhd;
    int      nbridmethod;
    void    *nmdata;
    int      weightmethod;
    void    *wmdata;
    int      sznotes;
    char    *notes;
  }
  VTX_INTERP_TBL, VIT;

  /* KNearest-Neighbor Structure */
  typedef struct
  {
    int   knearest;
    float dmax;
    int   noself;
  }
  VSS_KNN;

  /* Selavg */
  typedef struct
  {
    int    version;
    int    Nc, Nh;
    float  TR,TW,TPS;
    float  TER;
    int    Ntp;
    int    DTOrder;
    int    *dof;
    int    GammaFit;
    float *gfDelta;
    float *gfTau;
    int    NullCondId;
  }
  VSS_SELAVG;

  /* spatial interpolation table - projection*/
  typedef struct
  {
    int version;
    char sampname[1024];
    char inspacename[1024];
    int  inspacedim;
    char outspacename[1024];
    int  outspacedim;
    float **w;
    int  sznotes;
    char  *notes;
  }
  SPC_INTERP_TBL, SIT;

  /* Multi-linear interpolation */
  typedef struct
  {
    int   version;
    char  sampname[1024];
    int   spacedim; /* number of axes */
    int   nvertices; /* number of vertices in the uniform quantization */
    int   *nsamp; /* number of samples per axis */
    float *dsamp; /* distance between samples */
    float *val0;  /* axis value at sample 0 */
  }
  MULTI_LIN_INTERP, MLI;

  void *vssAllocVPM(int rows, int cols);
  int   vssFreeVPM(void **pm);
  int   vssDumpVPM(FILE *fp, char *fmt, void *pm);
  int  vssWriteVPM(FILE *fp, void *v);
  int   vssWriteVPMFile(char *filename, void *pm);
  void *vssReadVPM(FILE *fp);
  void *vssReadVPMFile(char *filename);
  int vssVPMRand(VPM *vpm, float min, float max, int seed);
  int vssVPMIdentity(VPM *vpm);

  VSS *vssProject(VPM *vpm, VSS *src);
  VPM *vssMatrix2VPM(MATRIX *m);
  MATRIX *vssVPM2Matrix(VPM *vpm);
  VPM *vssReadRegisterAsVPM(char *regfile);

  void *vssAllocVSSHeader(int version);
  int   vssAllocVSSData(void *vss);
  int   vssFreeVSS(void **vss);
  int   vssDumpVSSHeader(FILE *fp, void *vss);
  int   vssDumpVSSData(FILE *fp, char *fmt, void *vss);
  void *vssCopyVSSHeader(void *psrc);
  int   vssCopyVSSData(void *pdest, void *psrc);
  void *vssCopyVSS(void *src);
  int   vssDumpVtxNbrhd(FILE *fp, char *fmt, void *pvn);

  void *vssAllocVITHeader(int version);
  int   vssAllocVITData(void *vit);
  int   vssFreeVIT(void **vit);
  int   vssDumpVITHeader(FILE *fp, void *vit);
  int   vssDumpVITData(FILE *fp, char *fmt, void *vit);

  void *vssAllocSITHeader(int version);
  int   vssAllocSITData(void *sit);
  int   vssFreeSIT(void **sit);
  int   vssDumpSITHeader(FILE *fp, void *sit);
  int   vssDumpSITData(FILE *fp, char *fmt, void *sit);


  void *vssAllocMLIHeader(int version);
  int vssAllocMLIData(void *pmli);
  int vssDumpMLIHeader(FILE *fp, void *pmli);
  int vssDumpMLIData(FILE *fp, char *fmt, void *pmli);
  int vssSubscript2Vertex(int *ss, MLI *mli);
  int vssVertex2Subscript(int *ss, MLI *mli, int vtx);
  int vssFreeMLI(void **ppmli);
  int vssCheckSubscr2Vtx(FILE *fp, MLI *mli, int roe);
  int vssBuildMLINbrhd(VTXNBRHD *tvn, float *v, void *pmli);
  int vssMLINVertices(MLI *mli);


  int vssGetVtxNbrNum(VTXNBRHD *vn, int id);
  int vssAddVtxNbr(VTXNBRHD *vn, int id);

  VIT *vssBuildInterpTable(VSS *targ,        VSS *src,
                           int nbridmethod,  void *nmdata,
                           int weightmethod, void *wmdata);
  VIT *vssCompositeInterpTable(VIT *a2b, VIT *b2c);
  VIT *vssBuildMLITable(VSS *targ, void *pmli);

  VSS *vssResample(VIT *vit, VSS *src);
  VSS *vssResampleNearestNbr(VIT *vit, VSS *src);
  int vssCountHoles(VSS *vss);


  int vssNormalizeVector(float *v, int vdim, int NormOrder);
  int vssNormalizeSITWeight(SIT *sit, int NormOrder);
  int vssSITWeightIdentity(SIT *sit, float scale);
  int vssSITWeightURand(SIT *sit, float min, float max, long seed);
  VSS *vssInterpSpace(SIT *sit, VSS *inspc);


  int vssNbrIdMethod_Dist(VTXNBRHD *nbrhd, int tvtx, VSS *targ,
                          int candvtx, VSS *src, void *wm);
  int vssNbrIdMethod_KNearest(VTXNBRHD *nbrhd, int tvtx, VSS *targ,
                              int candvtx, VSS *src, void *wm);
  int vssNbrIdMethod_QNearest(VTXNBRHD *nbrhd, int tvtx, VSS *targ,
                              int candvtx, VSS *src, void *wm);

  int vssWeightMethod_NbrAvg(int tvtx, VSS *targ, VTXNBRHD *nbrhd,
                             VSS *src, int nbidmethod, void *nmdata,
                             void *wmdata);

  int vssRandFill(VSS *vss, float *min, float *max);
  VSS *vssRandomQuantization (MLI *mli, char *spacename);
  VSS *vssUniformQuantization(MLI *mli, char *spacename);


  VSS *vssSamplePoly(VSS *coord, int pdim);
  VSS *vssSampleGaussianNoise(char *sampname, int nvertices, int spcdim,
                              float mean, float var, long seed);
  VSS *vssSampleSphere(char *sampname, int nvertices, int spcdim,
                       float radius, long seed);

  VSS *vssVertexStat(VSS *inspc);
  float **vssSpaceStat(VSS *spc);
  float *vssRMSDataDiff(VSS *vss1, VSS *vss2);
  float *vssVtxStat(VSS *vss, int order, int avgholes);


  ISS *vssInterSpaceRandomEffects(VSS **spc, int nobservations,
                                  float *weight, int IncludeHoles);
  ISS *vssISAvgRE(VSS **spc, int nspaces, float *wspace, int AvgHoles);
  ISS *vssISAvgFE(ISS **spc, int nspaces, float *wspace, int AvgHoles);

  void *vssAllocISSHeader(int version);
  void *vssCopyISSHeader(void *dest, void *src);
  void *vssCopyISS(void *dest, void *src);
  int  vssFreeISS(void **ppiss);
  int vssAllocISSData(void *piss);
  int vssDumpISSData(FILE *fp, char *fmt, void *piss);
  int vssDumpISSHeader(FILE *fp, void *piss);
  int vssWriteISSFile(char *filename, void *piss);
  void *vssReadISSFile(char *filename);

  int   vssDumpVSTHeader(FILE *fp, void *pvst);
  int   vssFreeVST(void **ppvst);
  void *vssAllocVST(int version);
  int vssWriteVSTFile(char *filename, void *pvst);
  void *vssReadVSTFile(char *filename);

  VSS *vssWFormat2VSS(int num, int *vtxnum, float *wval);
  int vssWriteWFormat(char *wfile, VSS *vss, int component);
  VSS *vssSurf2VSS(MRI_SURFACE *surf, char *vtxname, char *spacename, char *sampname);

  VSS *vssLoadAsVSS(char *fname, char *format);
  VSS *vssReadWFormat(char *wfile);
  VSS *vssReadIcoFormat(char *fname);
  VSS *vssReadSurfFormat(char *surffilename);
  VSS *vssReadTalairachFormat(char *stem);
  VSS *vssReadBVolumeFormat(char *stem);
  int vssWriteBVolumeFormat(char *stem, VSS *vss,
                            int Nrows, int Ncols, int Nslices);
  int vssWriteSelavgFormat(char *stem, ISS *iss,
                           int Nrows, int Ncols, int Nslices);

  float *vssReadBFile(char *stem, char *ext, int slice,
                      int Nrows, int Ncols, int Ntp, int Endian);
  ISS *vssReadSelavgVolume(char *stem);
  VST *vssReadSelXAvgAsVST(char *stem);
  VSS *vssReadBVolumeFormat(char *stem);
  MLI *vssReadBVolumeAsMLI(char *stem);

  char *vssSurfaceFromPath(char *fname);
  char *vssSubjectFromPath(char *fname);
  char *vssHemifieldFromPath(char *fname);

  /* ------------ VSS File IO --------------------*/
  int   vssWriteVSSFile(char *fname, void *vss);
  void *vssReadVSSFile(char *fname);
  int   vssWriteVSSHeader(FILE *fp, void *v);
  void *vssReadVSSHeader(FILE *fp);
  int   vssWriteVSSData(FILE *fp, void *v);
  int   vssReadVSSData(FILE *fp, void *v);

  /* ------------ VIT File IO --------------------*/
  int   vssWriteVITFile(char *fname, void *vit);
  void *vssReadVITFile(char *fname);
  int   vssWriteVITHeader(FILE *fp, void *v);
  void *vssReadVITHeader(FILE *fp);
  int   vssWriteVITData(FILE *fp, void *v);
  int   vssReadVITData(FILE *fp, void *v);

  /* -------------------------------------------- */
  int byteswapbuf(void *buf, long int nbufbytes);
  long int getfilesize(char *fname);
  char *TimeStamp(void);

  /* -------------------------------------------- */
  void *vssAllocUserData(int udid);
  int   vssFreeUserData(void **ppud, int udid);
  void *vssCopyUserData(void *uddest, void *udsrc, int udid);
  int   vssWriteUserData(FILE *fp, void *pud, int udid);
  void *vssReadUserData(FILE *fp, int udid);
  char const *vssUserDataName(int udid);
  int vssDumpUserData(FILE *fp, char *fmt, void *pud, int udid);


  /* -------------------------------------------- */
  void *vssAllocKNearest(int version);
  int vssFreeKNearest(void **knn);
  int vssWriteKNearest(FILE *fp, void *pknn);
  void *vssReadKNearest(FILE *fp);
  void *vssCopyKNearest(void *dest, void *src);
  int vssDumpKNearest(FILE *fp, char *fmt, void *pknn);

  /* -------------------------------------------- */
  void *vssAllocSelavg(int version);
  int vssFreeSelavg(void **sa);
  int vssWriteSelavg(FILE *fp, void *psa);
  void *vssReadSelavg(FILE *fp);
  void *vssCopySelavg(void *dest, void *src);
  int vssDumpSelavg(FILE *fp, char *fmt, void *psa);

  int GetFileCode(char *fname);
  int CountBVolumeSlices(char *stem);
  float * vssReadBSlice(char *stem, char *ext, int slice,
                        int *Nrows, int *Ncols, int *Nframes);
  VPM *vssReadHCovAsVPM(char *stem, int slice);

  /* ----------------------------------------------------- */
#ifdef __cplusplus
}
#endif /* #ifdef __cplusplus */

#endif /* #ifndef VTXSAMPSPC_H */

