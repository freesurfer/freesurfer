#ifndef VENTFIX_H
#define VENTFIX_H

#include "mri.h"
#include "volcluster.h"

class VentFix
{
public:
  static MRI *fixasegps(MRI *asegps, MRI *brainmask, char *segids, float threshmin, int niters, int nmax, int topo);
 
#if 0
  VentFix();
  ~VentFix();

  void init(const char* outdir, const char* subject, const char *asegps, float threshmin0, float threshmax0 = 0);
  void binarize();
  void apply_mask(); 

  void create_volcluster(const char *sumfile, const char *outcnid);
  void expandSegIndices(int segid, int niters, int nmax, int topo, const char *ventfilldat, const char *newseg, const char* ventfilljson);

  void reset();

  static void fixasegps(MRI *asegps, MRI *brainmask, MRI *newseg, const char *segids, float threshmin, int niters, int nmax, int topo);
 
private:
  //int nthreads;

  MRI *asegVol;
  MRI *binVol;
  MRI *maskVol;
  MRI *ocnVol;
  MRI *newsegVol;

  int nClusters;
  VOLCLUSTER **ClusterList;

  char* sumfile;
  char* outocn;
  float threshmin;
  float threshmax;

  char mdir[256];          // subject mri directory
  char asegpspath[256];    // full path to aseg in subject mri directory
  char brainmask[256];     // brainmask.mgz in subject mri directory
  char newseg[256];        // 
  char outbin[256];        // binarized aseg bin.mgh

  // mri_volcluster parameters might be useful in the future
  CSD *csd;
  double fwhm;
  int FixMNI;
  int BonferroniMax;
  int sig2pmax;
  char *segvolfile;
  MRI *segvol;
  COLOR_TABLE *segctab;
  char *regfile;

private:
  void _outClusterVol(const char *outcnid);   // write clusters numbers to a volume, include color LUT
  void _outsumfile(const char *sumfile);      // output summary file ocn.sum.dat
#endif
};

#endif
