/**
 * @brief fuse a set of segmentations (asegs)
 *
 * program to fuse a group of cross sectional segmentations into
 * an initial estimate of a longitudinal one
 * See Sabuncu et al., MICCA 2009 (SNIP paper).
 *
 *
 * Original Author: Bruce Fischl
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

#include <string>
#include <vector>
#include <iostream>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "cma.h"
#include "transform.h"

struct Parameters
{
  std::string prog_name;
  std::vector<std::string> norm_paths;
  std::vector<std::string> aseg_paths;
  std::vector<std::string> aseg_nocc_paths;
  std::vector<std::string> trx_paths;
  double cross_time_sigma = 3.0;
  std::string in_path;
  std::string out_path;
};


static void parseCommand(int argc, char *argv[], Parameters &par);
static void parseList(int &argc, char **&argv, std::vector<std::string> &list);
static MRI *readMriOrError(const std::string path);
static MRI *warpReplace(TRANSFORM *trx, MRI* mri_src, const MRI* mri_dst,
                        int interp_method);
static MRI *MRIfuseSegmentations(MRI *mri_in,
                                 std::vector<MRI*> mri_asegs,
                                 std::vector<MRI*> mri_nocc_asegs,
                                 std::vector<MRI*> mri_norms,
                                 double sigma);
static void forward(int &argc, char **&argv);
static void printUsage(void);


int main(int argc, char *argv[])
{
  int nargs = handleVersionOption(argc, argv, "mri_fuse_segmentations");
  if (nargs && argc - nargs == 1) exit(EXIT_SUCCESS);
  argc -= nargs;
  
  Timer start;
  start.reset();
  Parameters par;
  parseCommand(argc, argv, par);
  if (par.trx_paths.empty()) {
    std::cout << "No transforms: assuming images are in same space\n";
  }
  
  // Prepare data: try to be efficient by not keeping GCAMs in memory.
  MRI* mri_in = readMriOrError(par.in_path);
  const size_t num_tp = par.aseg_paths.size();
  std::vector<MRI*> mri_asegs;
  std::vector<MRI*> mri_nocc_asegs;
  std::vector<MRI*> mri_norms;
  for (size_t i=0; i<num_tp; i++) {
    std::cout << "Loading data for TP " << i+1 << std::endl;
    mri_asegs.push_back( readMriOrError(par.aseg_paths[i]) );
    mri_nocc_asegs.push_back( readMriOrError(par.aseg_nocc_paths[i]) );
    mri_norms.push_back( readMriOrError(par.norm_paths[i]) );
    if (par.trx_paths.empty()) {
      continue;
    }
    else if (!par.trx_paths[i].compare("identity.nofile")) {
      std::cout << "Passing on identity.nofile for TP " << i+1 << std::endl;
      continue;
    }
    std::cout << "Resampling data for TP " << i+1 << std::endl;
    TRANSFORM *trx = TransformRead( par.trx_paths[i].c_str() );
    if (!trx) {
      ErrorExit(ERROR_NOFILE, "ERROR: could not read transform %s",
                par.trx_paths[i].c_str());
    }
    mri_asegs[i] = warpReplace(trx, mri_asegs[i], mri_in, SAMPLE_NEAREST);
    mri_nocc_asegs[i] = warpReplace(trx, mri_nocc_asegs[i], mri_in,
                                    SAMPLE_NEAREST);
    // Prevent rounding of interpolated voxel intensities:
    MRI *mri_tmp = mri_norms[i];
    mri_norms[i] = MRIchangeType(mri_norms[i], MRI_FLOAT, 0, 0, 0);
    MRIfree(&mri_tmp);
    mri_norms[i] = warpReplace(trx, mri_norms[i], mri_in, SAMPLE_TRILINEAR);
    TransformFree(&trx);
  }
  
  std::cout << "Fusing...\n";
  MRI *mri_fused = MRIfuseSegmentations(mri_in, mri_asegs, mri_nocc_asegs,
                                        mri_norms, par.cross_time_sigma);
  std::cout << "Writing fused segmentation to " << par.out_path << std::endl;
  MRIwrite(mri_fused, par.out_path.c_str());
  int msec = start.milliseconds();
  int sec = nint((float)msec/1000.0f);
  int min = sec / 60;
  sec = sec % 60;
  std::cout << "Fusing took " << min << " minutes and " << sec << " seconds\n";
  
  MRIfree(&mri_in);
  MRIfree(&mri_fused);
  for (size_t i = 0; i < num_tp; i++) {
    MRIfree(&mri_asegs[i]);
    MRIfree(&mri_nocc_asegs[i]);
    MRIfree(&mri_norms[i]);
  }

  printf("#VMPC# mri_fuse_segmentations VmPeak  %d\n",GetVmPeak());
  printf("mri_fuse_segmentations done\n");
  return (NO_ERROR);
}


static void parseCommand(int argc, char *argv[], Parameters &par)
{
  par.prog_name = std::string(argv[0]);
  forward(argc, argv);

  if (argc == 0) {
    printUsage();
    exit(EXIT_FAILURE);
  }
  if (argc < 2 || ISOPTION(*argv[argc-1]) || ISOPTION(*argv[argc-2])) {
    ErrorExit(ERROR_BADPARM, "ERROR: positional arguments not specified");
  }
  par.out_path = std::string(argv[--argc]);
  par.in_path = std::string(argv[--argc]);
  std::cout << "Input volume: " << par.in_path << std::endl;
  std::cout << "Output aseg: " << par.out_path << std::endl;
  
  while (argc > 0) {
    if (!strcmp(*argv, "--aseg") || !strcmp(*argv, "-a")) {
      parseList(argc, argv, par.aseg_paths);
    }
    else if (!strcmp(*argv, "--nocc") || !strcmp(*argv, "-c")) {
      parseList(argc, argv, par.aseg_nocc_paths);
    }
    else if (!strcmp(*argv, "--norm") || !strcmp(*argv, "-n")) {
      parseList(argc, argv, par.norm_paths);
    }
    else if (!strcmp(*argv, "--trx") || !strcmp(*argv, "-t")) {
      parseList(argc, argv, par.trx_paths);
    }
    else if (!strcmp(*argv, "--sigma") || !strcmp(*argv, "-s")) {
      forward(argc, argv);
      if (argc==0 || ISOPTION(*argv[0])) {
        ErrorExit(ERROR_BADPARM, "ERROR: no cross-time sigma specified");
      }
      par.cross_time_sigma = atof(*argv);
      forward(argc, argv);
    }
    else if (!strcmp(*argv, "--debug") || !strcmp(*argv, "-d")) {
      forward(argc, argv);
      size_t coord[3];
      for (int i=0; i<3; i++) {
        if (argc==0 || ISOPTION(*argv[0])) {
          ErrorExit(ERROR_BADPARM, "ERROR: no debug voxel specified");
        }
        coord[i] = atoi(*argv);
        forward(argc, argv);
      }
      Gx = coord[0];
      Gy = coord[1];
      Gz = coord[2];
    }
    else {
      ErrorExit(ERROR_BADPARM, "ERROR: unknown argument %s", *argv);
    }
  } // End of while.
  
  const size_t num_tp = par.aseg_paths.size();
  if (num_tp == 0) {
    ErrorExit(ERROR_BADPARM, "ERROR: no asegs specified");
  }
  if (par.norm_paths.size() != num_tp) {
    ErrorExit(ERROR_BADPARM,
              "ERROR: number of asegs and norm volumes do not match");
  }
  if (par.aseg_nocc_paths.size() != num_tp) {
    ErrorExit(ERROR_BADPARM,
              "ERROR: number of asegs and no-CC-label asegs do not match");
  }
  if (!par.trx_paths.empty() && par.trx_paths.size() != num_tp) {
    ErrorExit(ERROR_BADPARM,
              "ERROR: number of asegs and transforms do not match");
  }
}


#include "mri_fuse_segmentations.help.xml.h"
static void printUsage(void)
{
  outputHelpXml(mri_fuse_segmentations_help_xml,
                mri_fuse_segmentations_help_xml_len);
}


static void forward(int &argc, char **&argv)
{
  argc--;
  argv++;
}


static void parseList(int &argc, char **&argv, std::vector<std::string> &list)
{
  std::string option(*argv);
  forward(argc, argv);
  while (argc>0 && !ISOPTION(*argv[0])) {
    list.push_back( std::string(*argv) );
    forward(argc, argv);
  }
}


static MRI *readMriOrError(const std::string path)
{
  MRI *mri = MRIread( path.c_str() );
  if (!mri) {
    exit(EXIT_FAILURE); // Underlying functions print info.
  }
  return (mri);
}


static MRI *warpReplace(TRANSFORM *trx, MRI *mri_src, const MRI *mri_dst,
                        int interp_method)
{
  if (trx->type != MORPH_3D_TYPE) {
    LTA *lta = (LTA *)trx->xform;
    trx->xform = (void *)LTAreduce(lta); // Apply full transform.
    LTAfree(&lta);
  }
  vg_isEqual_Threshold = 10e-4; // Override, include/transform.h.
  VOL_GEOM vg_src, vg_dst, vg_trx_src, vg_trx_dst;
  TransformGetSrcVolGeom(trx, &vg_trx_src);
  TransformGetDstVolGeom(trx, &vg_trx_dst);
  getVolGeom(mri_src, &vg_src);
  getVolGeom(mri_dst, &vg_dst);
  if (!vg_isEqual(&vg_src, &vg_trx_src) || !vg_isEqual(&vg_dst, &vg_trx_dst)) {
    if (vg_isEqual(&vg_src, &vg_trx_dst) && vg_isEqual(&vg_dst, &vg_trx_src)) {
      // Pass MRI in case path changed (see GCAMfillInverse):
      std::cout << "Inverting transform to match MRI geometries\n";
      TransformInvertReplace(trx, mri_dst);
    }
    else {
      ErrorExit(ERROR_BADPARM, "ERROR: transform geometry does not match");
    }
  }
  MRI *mri = mri_src;
  mri = TransformApplyType(trx, mri, NULL, interp_method);
  MRIfree(&mri_src);
  return (mri);
}


static MRI *MRIfuseSegmentations(MRI *mri_in,
                                 std::vector<MRI*> mri_asegs,
                                 std::vector<MRI*> mri_nocc_asegs,
                                 std::vector<MRI*> mri_norms,
                                 double sigma)
{
  std::vector<double> label_pvals(MAX_CMA_LABELS); // Initialized to zero.
  int cc_relabel_count = 0;
  const double s = -0.5 / (sigma*sigma);
  const int frame = 0;
  MRI *mri_fused = MRIclone(mri_in, NULL);
  printf("Smoothing temporal bias field with sigma = %2.1f\n", sigma);

  for (int x = 0 ; x < mri_fused->width ; x++) {
    for (int y = 0 ; y < mri_fused->height ; y++) {
      for (int z = 0 ; z < mri_fused->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        double val = MRIgetVoxVal(mri_in, x, y, z, frame) ;
        int min_label = MAX_CMA_LABELS ;
        int max_label = 0 ;
        for (size_t i = 0 ; i < mri_asegs.size() ; i++) {
          double oval = MRIgetVoxVal(mri_norms[i], x, y, z, frame) ;
          int label = MRIgetVoxVal(mri_asegs[i], x, y, z, frame) ;
          switch (label) {
            case CC_Posterior:
            case CC_Mid_Posterior:
            case CC_Central:
            case CC_Mid_Anterior:
            case CC_Anterior:
              // need to relabel these back to white matter, from noCC aseg
              // (downstream mri_cc will re-seg the CC)
              int label_nocc = MRIgetVoxVal(mri_nocc_asegs[i], x, y, z, frame);
              if (label_nocc == Left_Cerebral_White_Matter ||
                  label_nocc == Right_Cerebral_White_Matter) {
                label = label_nocc;
                cc_relabel_count++;
              }
            break;
          }
          if (label < min_label) {
            min_label = label;
          }
          if (label > max_label) {
            max_label = label;
          }
          double dif = oval - val ;
          // if sigma is zero only trust images with identical intensity
          double p = 0;
          const double threshold = 0.0000001;
          if (fabs(sigma) < threshold)
          {
            p = fabs(dif) < threshold ? 1.0 : 0.0;
          }
          else {
            p = exp((dif*dif)*s);
          }
          label_pvals[label] += p ;
        }
        double p = 0;
        // select label with highest probability
        for (int label = min_label ; label <= max_label ; label++) {
          if (label_pvals[label] > p) {
            p = label_pvals[label] ;
            MRIsetVoxVal(mri_fused, x, y, z, frame, label) ;
          }
          label_pvals[label] = 0 ; // reset for next time to avoid big memset
        }
      }
    }
  }
  
  // confirm that no CC labels still exist in the fused volume
  for (int x = 0 ; x < mri_fused->width ; x++) {
    for (int y = 0 ; y < mri_fused->height ; y++) {
      for (int z = 0 ; z < mri_fused->depth ; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        int label = MRIgetVoxVal(mri_fused, x, y, z, frame) ;
        switch (label) {
          case CC_Posterior:
          case CC_Mid_Posterior:
          case CC_Central:
          case CC_Mid_Anterior:
          case CC_Anterior:
            // look at neighbors to determine who is more predominant: left or
            // right white matter voxels; predominance determines relabeling
            int leftwm=0;
            int rightwm=0;
            int xk, yk, zk, xi, yi, zi;
            // brute force 3x3 search:
            for (zk = -1 ; zk <= 1 ; zk++) {
              zi = mri_fused->zi[z+zk] ;
              for (yk = -1 ; yk <= 1 ; yk++) {
                yi = mri_fused->yi[y+yk] ;
                for (xk = -1 ; xk <= 1 ; xk++) {
                  xi = mri_fused->xi[x+xk] ;
                  if (!zk && !yk && !xk) continue ;
                  int lbl = MRIgetVoxVal(mri_fused, xi, yi, zi, frame);
                  if (lbl == Left_Cerebral_White_Matter) leftwm++;
                  else if (lbl == Right_Cerebral_White_Matter) rightwm++;
                }
              }
            }
            //printf("Neighbors: leftwm=%d, rightwm=%d\n", leftwm, rightwm);
            if (leftwm > rightwm) {
              MRIsetVoxVal(mri_fused, x, y, z, frame, Left_Cerebral_White_Matter);
              cc_relabel_count++;
            }
            else if (rightwm > 0) {
              MRIsetVoxVal(mri_fused, x, y, z, frame, Right_Cerebral_White_Matter);
              cc_relabel_count++;
            }
            else {
              printf("ERROR: could not determine a label for relabeling\n");
              printf("CC label (%d) at xyz %d %d %d\n", label, x, y, z);
            }
          break;
        }
      }
    }
  }
  printf("Relabeled %d CC labels\n", cc_relabel_count);
  return (mri_fused);
}

