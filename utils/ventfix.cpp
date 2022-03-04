#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "ventfix.h"
#include "colortab.h"
#include "mri_identify.h"
#include "diag.h"
#include "mri2.h"

MRI* VentFix::fixasegps(MRI *asegps, MRI *brainmask, char *segids, float threshmin, int niters, int nmax, int topo)
{
  // Prepare the output binary volume
  MRI *binVol = MRIallocSequence(asegps->width, asegps->height, asegps->depth, MRI_INT, 1);
  if (binVol == NULL) exit(1);
  MRIcopyHeader(asegps, binVol);

  // binarize input asegps
  MRIbinarize(asegps, binVol, threshmin, 1, 0);

  // apply brainmask to binary volume
  MRI *newbinVol = MRImask(binVol, brainmask, NULL, 0, 0);
  MRIfree(&binVol);
  binVol = newbinVol;

  int nClusters;
  VOLCLUSTER **ClusterList = clustGetClusters(binVol, 0, threshmin, 0, 0, 0, NULL, &nClusters, NULL);

  MRI *ocnVol = clustClusterList2Vol(ClusterList, nClusters, binVol, 0, 0);

  // mgz doesn't have cluster information.
  // ??? do we need to do this if we are not outputing .lut ???
#if 1
  ocnVol->ct = CTABalloc(nClusters+1);
  strcpy(ocnVol->ct->entries[0]->name, "Unknown");

  int n;
  for (n = 0; n < nClusters; n++)
    sprintf(ocnVol->ct->entries[n+1]->name, "Cluster-%03d", n+1);
#endif

#if 0  // testing
  CTABwriteFileASCII(ocnVol->ct, "ocn.lut");
#endif

  // free binary volume
  MRIfree(&binVol);

  // free ClusterList
  clustFreeClusterList(&ClusterList, nClusters);


  int nc;
  PointSet centroid;
  MRI *segVol = NULL, *newsegVol = NULL;

  char *restsegids = (char*)segids;
  char *nextsegid = strtok_r(segids, ",", &restsegids);
  while (nextsegid != NULL)
  {
    if (segVol == NULL)
      segVol = asegps;
    else
      segVol = newsegVol;

    int segid = atoi(nextsegid);
    nextsegid = strtok_r(NULL, ",", &restsegids);
    printf("\nexpandSegIndices for segid %d\n", segid);
    newsegVol = ExpandSegIndices(segVol, segid, ocnVol, niters, nmax, topo, &nc, centroid, NULL);
    printf("segid %d niters %d nmax %d topo %d   nc %5d \n", segid, niters, nmax, topo, nc);

    // free MRI structures returned by previous ExpandSegIndices() call
    if (segVol != asegps)
      MRIfree(&segVol);
  }

  MRIfree(&ocnVol);

  return newsegVol;
}

#if 0
VentFix::VentFix()
{
  // these are for future use, set these in init(), so each run can be different ???
  csd = NULL;
  fwhm = -1;
  FixMNI = 1;
  BonferroniMax = 0;
  sig2pmax = 0; // convert max value from -log10(p) to p
  segvolfile = NULL;
  segvol = NULL;
  segctab = NULL;
  regfile = NULL;

  memset(mdir, '\0', sizeof(mdir));
  memset(asegpspath, '\0', sizeof(asegpspath));
  memset(brainmask, '\0', sizeof(brainmask));
  memset(newseg, '\0', sizeof(newseg));
  memset(outbin, '\0', sizeof(outbin));

  threshmin = 0;
  threshmax = 0;

  asegVol   = NULL; 
  binVol    = NULL;
  maskVol   = NULL;
  ocnVol    = NULL;
  newsegVol = NULL;
}

VentFix::~VentFix()
{
  if (maskVol != NULL)
    MRIfree(&maskVol);

  if (binVol != NULL)
    MRIfree(&binVol);

  if (asegVol != NULL)
    MRIfree(&asegVol);

  if (ocnVol != NULL)
    MRIfree(&ocnVol);

  if (newsegVol != NULL)
    MRIfree(&newsegVol);

  // free ClusterList
  clustFreeClusterList(&ClusterList, nClusters);
}

void VentFix::reset()
{
  if (maskVol != NULL)
    MRIfree(&maskVol);

  if (binVol != NULL)
    MRIfree(&binVol);

  if (asegVol != NULL)
    MRIfree(&asegVol);

  if (ocnVol != NULL)
    MRIfree(&ocnVol);

  if (newsegVol != NULL)
    MRIfree(&newsegVol);
  newsegVol = NULL;

  memset(mdir, '\0', sizeof(mdir));
  memset(asegpspath, '\0', sizeof(asegpspath));
  memset(brainmask, '\0', sizeof(brainmask));
  memset(newseg, '\0', sizeof(newseg));
  memset(outbin, '\0', sizeof(outbin));

  threshmin = 0;
  threshmax = 0;
}

void VentFix::init(const char* outdir, const char* subject, const char* asegps, float threshmin0, float threshmax0)
{
  threshmin = threshmin0;
  threshmax = threshmax0;

  //set mdir = $SUBJECTS_DIR/$subject/mri/
  //set asegpspath = $mdir/$asegps
  //set brainmask = $mdir/brainmask.mgz
  //set newseg = $outdir/newseg.mgz
  //set outbin = $tmpdir/bin.mgh

  char *subjdir = getenv("SUBJECTS_DIR");
  sprintf(mdir, "%s/%s/mri", subjdir, subject);
  sprintf(asegpspath, "%s/%s", mdir, asegps);
  sprintf(brainmask, "%s/brainmask.mgz", mdir);

  //char *pwd = get_current_dir_name();
  //sprintf(outbin, "%s/tmp/bin.mgh", pwd);
}

void VentFix::binarize()
{
  //sprintf(asegpspath, "%s", "/autofs/space/sulc_001/users/ventfix/yujing/0619_m24/mri/aseg.presurf.mgz");
  sprintf(outbin, "bin.0.mgh");

  // initialize asegVol
  asegVol = MRIread(asegpspath);

  // Prepare the output binary volume
  binVol = MRIallocSequence(asegVol->width, asegVol->height, asegVol->depth, MRI_INT, 1);
  if (binVol == NULL) exit(1);
  MRIcopyHeader(asegVol, binVol);

  MRIbinarize(asegVol, binVol, threshmin, 1, 0);

  printf("Writing output to %s\n", outbin);
  int err = MRIwrite(binVol, outbin);
  if(err) 
    return;

  /* print report at the end: mri_binarize.cpp::dump_options()
   * input      /autofs/space/sulc_001/users/ventfix/yujing/0619_m24/mri//aseg.presurf.mgz
frame      0
nErode3d   0
nErode2d   0
output     tm/bin.mgh
Binarizing based on threshold
min        0.5
max        +infinity
binval        0
binvalnot     1
fstart = 0, fend = 0, nframes = 1
Starting parallel 1
Found 1268197 values in range
Counting number of voxels in first frame
Found 1268196 voxels in final mask
Writing output to tm/bin.mgh
Count: 1268196 1268196.000000 16777216 7.559037
   */
}

void VentFix::apply_mask()
{
  //sprintf(brainmask, "/autofs/space/sulc_001/users/ventfix/yujing/0619_m24/mri/brainmask.mgz");

  // read brainmask.mgz
  maskVol = MRIread(brainmask);

  MRI *newbinVol = MRImask(binVol, maskVol, NULL, 0, 0);

  sprintf(outbin, "bin.mgh");
  printf("Writing output to %s\n", outbin);
  int err = MRIwrite(newbinVol, outbin);
  if(err) 
    return;

  MRIfree(&binVol);
  binVol = newbinVol;
}

// we don't need binVol, maskVol, ClusterList after this call.
// free them at the end??? 
void VentFix::create_volcluster(const char *sumfile, const char *outcnid)
{
  ClusterList = clustGetClusters(binVol, 0, threshmin, threshmax, 0, 0, NULL, &nClusters, NULL);

  printf("INFO: Found %d clusters\n", nClusters);

  if (outcnid != NULL)
    _outClusterVol(outcnid);

  if (sumfile != NULL)
    _outsumfile(sumfile);
}


// write clusters numbers to a volume, include color LUT
// extracted from mri_volcluster.cpp - if (outcnid != 0)
void VentFix::_outClusterVol(const char *outcnid)
{
  ocnVol = clustClusterList2Vol(ClusterList, nClusters, binVol, 0, 0);
  ocnVol->ct = CTABalloc(nClusters+1);
  strcpy(ocnVol->ct->entries[0]->name, "Unknown");

  int n;
  for (n = 0; n < nClusters; n++)
    sprintf(ocnVol->ct->entries[n+1]->name, "Cluster-%03d", n+1);

  char tmpstr[2000] = {'\0'};
  char *stem = IDstemFromName(outcnid);
  sprintf(tmpstr, "%s.lut", stem);

  CTABwriteFileASCII(ocnVol->ct, tmpstr);
  printf("INFO: writing OCN to %s\n", outcnid);
  MRIwrite(ocnVol, outcnid);

  free(stem);
}

// output summary file ocn.sum.dat
// this is extracted from mri_volcluster.cpp
// the summary is commented out for now
void VentFix::_outsumfile(const char *sumfile)
{
  FILE *fpsum = NULL;

  /* Open the Summary File (or set its pointer to stdout) */
  if (sumfile != NULL) {
    fpsum = fopen(sumfile,"w");
    if (fpsum == NULL) {
      printf("ERROR: could not open %s for writing\n",sumfile);
      exit(1);
    }
    printf("INFO: writing summary to %s\n",sumfile);
  } else fpsum = stdout;

  float colres   = binVol->xsize;
  float rowres   = binVol->ysize;
  float sliceres = binVol->zsize;
  float voxsize  = colres * rowres * sliceres;

#if 0
  /* Dump summary to file or stdout */
  fprintf(fpsum,"# Cluster Growing Summary (mri_volcluster)\n");
  fprintf(fpsum,"# %s\n",getVersion().c_str());
  fprintf(fpsum,"# cwd %s\n",cwd);
  fprintf(fpsum,"# cmdline %s\n",cmdline);
  if(SUBJECTS_DIR) fprintf(fpsum,"# SUBJECTS_DIR  %s\n",SUBJECTS_DIR);
  fprintf(fpsum,"# sysname  %s\n",uts.sysname);
  fprintf(fpsum,"# hostname %s\n",uts.nodename);
  fprintf(fpsum,"# machine  %s\n",uts.machine);
  fprintf(fpsum,"# user     %s\n",VERuser());

  fprintf(fpsum,"# Input Volume:      %s\n",volid);
  fprintf(fpsum,"# Frame Number:      %d\n",frame);
  fprintf(fpsum,"# VoxSize_mm3 %g\n",voxsize);
  fprintf(fpsum,"# SearchSpace_mm3 %g\n",searchspace);
  fprintf(fpsum,"# SearchSpace_vox %d\n",nmask);
  fprintf(fpsum,"# Minimum Threshold: %g\n",threshmin);
  fprintf(fpsum,"# Bonferroni %d\n",Bonferroni);
  if (threshmax < 0)
    fprintf(fpsum,"# Maximum Threshold: inifinity\n");
  else
    fprintf(fpsum,"# Maximum Threshold: %g\n",threshmax);
  fprintf(fpsum,"# Threshold Sign:    %s\n",signstring);
  fprintf(fpsum,"# AdjustThreshWhenOneTail %d\n",AdjustThreshWhenOneTail);

  if (distthresh > 0)
    fprintf(fpsum,"# Distance Threshold: %g (mm)\n",distthresh);

  if(cwpvalthresh > 0)
    fprintf(fpsum,"# CW PValue Threshold: %g \n",cwpvalthresh);

  fprintf(fpsum,"# Size Threshold:    %g mm^3\n",sizethresh);
  fprintf(fpsum,"# Size Threshold:    %g voxels\n",sizethresh/voxsize);
  fprintf(fpsum,"# Voxel Size:        %g mm^3\n",voxsize);
  if (regfile) fprintf(fpsum,"# Registration:      %s\n",regfile);
  else fprintf(fpsum,"# Registration:      None : Tal Coords invalid\n");
  if (synthfunction != NULL)
    fprintf(fpsum,"# Synthesize:        %s\n",synthfunction);
  if (maskid != NULL) {
    fprintf(fpsum,"# Mask Vol:          %s\n",maskid);
    fprintf(fpsum,"# Mask Thresh:       %f\n",maskthresh);
    fprintf(fpsum,"# Mask Sign:         %s\n",masksignstring);
    fprintf(fpsum,"# Mask Invert:       %d\n",maskinvert);
  }
  fprintf(fpsum,"# AllowDiag:         %d\n",allowdiag);
  fprintf(fpsum,"# NClusters          %d\n",nclusters);
  if (csd != NULL) {
    fprintf(fpsum,"# CSD thresh  %lf\n",csd->thresh);
    fprintf(fpsum,"# CSD nreps    %d\n",csd->nreps);
    fprintf(fpsum,"# CSD simtype  %s\n",csd->simtype);
    fprintf(fpsum,"# CSD contrast %s\n",csd->contrast);
    fprintf(fpsum,"# CSD confint  %lf\n",ciPct);
  }
  if (fwhm > 0) {
    fprintf(fpsum,"# FWHM        %lf\n",fwhm);
    // dLh is an FSL parameter. Include here for comparison
    fprintf(fpsum,"# dLh         %lf\n", pow((fwhm/colres)/sqrt(4.0*log(2.0)),-3) );
  }

  fprintf(fpsum,"# \n");
#endif
  if (regfile) {
    if (FixMNI) {
      fprintf(fpsum,"# Reporting Coordinates in Talairach Space\n");
      fprintf(fpsum,"# Cluster   Size(n)   Size(mm^3)     "
              "TalX   TalY    TalZ              Max");
    } else {
      fprintf(fpsum,"# Reporting Coordinates in MNI305 Space\n");
      fprintf(fpsum,"# Cluster   Size(n)   Size(mm^3)     "
              "MNIX   MNIY    MNIZ              Max");
    }
  } else {
    fprintf(fpsum,"# Reporting Coordinates in Voxel Indices\n");
    fprintf(fpsum,"# Cluster   Size(n)   Size(mm^3)     "
            "VoxX    VoxY    VoxZ             Max");
  }

  if (csd != NULL)  fprintf(fpsum,"    CWP    CWPLow    CWPHi\n");
  else if (fwhm > 0) fprintf(fpsum,"     GRFCWP\n");
  else fprintf(fpsum,"\n");

  int n, col, row, slc;
  float x, y, z;
  MATRIX *CRS2MNI = MatrixIdentity(4, NULL);
  for (n = 0; n < nClusters; n++) {
    double maxval = ClusterList[n]->maxval;
    if(sig2pmax) {
      maxval = pow(10.0,-fabs(ClusterList[n]->maxval));
      if(BonferroniMax > 1) maxval *= BonferroniMax;
    }
    clustComputeTal(ClusterList[n], CRS2MNI); /* for "true" Tal coords */
    //clustComputeXYZ(ClusterList[n],CRS2FSA); /* for FSA coords */
    //clustComputeXYZ(ClusterList[n],CRS2MNI); /* for MNI coords */
    col = ClusterList[n]->col[ClusterList[n]->maxmember];
    row = ClusterList[n]->row[ClusterList[n]->maxmember];
    slc = ClusterList[n]->slc[ClusterList[n]->maxmember];
    if (regfile) {
      x = ClusterList[n]->x[ClusterList[n]->maxmember];
      y = ClusterList[n]->y[ClusterList[n]->maxmember];
      z = ClusterList[n]->z[ClusterList[n]->maxmember];
    } else {
      x=col;
      y=row;
      z=slc;
    }
    fprintf(fpsum,"%3d        %5d      %8.1f    %7.2f %7.2f %7.2f   %15.5f",
            n+1, ClusterList[n]->nmembers, voxsize*ClusterList[n]->nmembers,
            x, y, z, maxval);
    if (Gdiag_no > 0) fprintf(fpsum,"  %3d %3d %3d \n", col, row, slc);
    if (csd != NULL)
      fprintf(fpsum,"  %7.5lf  %7.5lf  %7.5lf",
              ClusterList[n]->pval_clusterwise,
              ClusterList[n]->pval_clusterwise_low,
              ClusterList[n]->pval_clusterwise_hi);
    else if (fwhm > 0)
      fprintf(fpsum, "  %9.7lf", ClusterList[n]->pval_clusterwise);
    if(segvolfile){
      int ctabindex = MRIgetVoxVal(segvol, col, row, slc, 0);
      fprintf(fpsum,"  %s", segctab->entries[ctabindex]->name);
    }
    fprintf(fpsum,"\n");
  }

  if (sumfile != NULL)
    fclose(fpsum);
}

// we only need asegVol and ocnVol for this call
// the rest MRI structures can be freed???
void VentFix::expandSegIndices(int segid, int niters, int nmax, int topo, const char *ventfilldat, const char *newseg, const char* ventfilljson)
{
  int nc;
  PointSet centroid;

  MRI *segVol = asegVol;
  if (newsegVol != NULL)
    segVol = newsegVol;

  newsegVol = ExpandSegIndices(segVol, segid, ocnVol, niters, nmax, topo, &nc, centroid, NULL);
  printf("segid %d niters %d nmax %d topo %d   nc %5d \n", segid, niters, nmax, topo, nc);

  if (ventfilldat != NULL)
  {
    FILE *fp = fopen(ventfilldat, "w");
    fprintf(fp,"segid %d niters %d nmax %d topo %d   nc %5d \n", segid, niters, nmax, topo, nc);
    fclose(fp);
  }

  if (newseg != NULL)
  {
    printf("Wrting to %s\n", newseg);
    MRIwrite(newsegVol, newseg);
  }

  centroid.save(ventfilljson);
}
#endif
