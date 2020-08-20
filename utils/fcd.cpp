/**
 * @brief I/O and analysis algorithms for FCDs (focal cortical dysplasias)
 *
 */
/*
 * Original Author: Bruce Fischl
 *
 * Copyright © 2013-2014 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "fcd.h"
#include "cma.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "label.h"
#include "macros.h"
#include "mrisegment.h"
#include "utils.h"

#include "romp_support.h"

static int sort_labels(FCD_DATA *fcd);
static int most_frequent_label(MRI *mri_seg, MRI_SEGMENT *mseg);
static int compare_labels(const void *v1, const void *v2);
static int fcdFreeLabels(FCD_DATA *fcd);
static MRI *build_distance_by_intensity_histo(
    MRI *mri_norm, MRI *mri_dist, MRI *mri_aseg, double dist_spacing, double max_dist);
static int augment_thicknesses(FCD_DATA *fcd, MRI *mri_pvals, double min_dist, double max_dist, double thresh);

static int most_frequent_label(MRI *mri_seg, MRI_SEGMENT *mseg)
{
  int label_counts[MAX_CMA_LABELS], i, max_count, max_label, label;

  memset(label_counts, 0, sizeof(label_counts));
  for (i = max_count = max_label = 0; i < mseg->nvoxels; i++) {
    label = MRIgetVoxVal(mri_seg, mseg->voxels[i].x, mseg->voxels[i].y, mseg->voxels[i].z, 0);
    if (IS_WM(label) == 0 && IS_UNKNOWN(label) == 0) {
      label_counts[label]++;
    }
  }
  for (i = max_count = max_label = 0; i < MAX_CMA_LABELS; i++) {
    if (label_counts[i] > max_count) {
      max_count = label_counts[i];
      max_label = i;
    }
  }
  return (max_label);
}

#define MAX_DIST 10
#define DIST_SPACING .5

FCD_DATA *FCDloadData(char *sdir, char *subject, char *suffix_in)
{
  FCD_DATA *fcd;
  MRI *mri_interior, *mri_dist, *mri_int_lh, *mri_int_rh, *mri_pvals;
  std::string fname, suffix;
  if( suffix_in != nullptr ) {
    suffix = std::string(suffix_in);
    suffix.insert(suffix.begin(), '.');
  }

  fcd = (FCD_DATA *)calloc(1, sizeof(FCD_DATA));

  std::string baseName = std::string(sdir) + "/" + std::string(subject) + "/surf/";
  fname = baseName + "lh.white" + suffix;
  if (!FileExists(fname.c_str())) {
    fname = baseName + "lh.white";
  }
  fcd->mris_lh = MRISread(fname.c_str());
  if (fcd->mris_lh == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }
  MRISsaveVertexPositions(fcd->mris_lh, WHITE_VERTICES);

  std::string pialname("pial");
  if (suffix_in) {
    pialname += std::string(".") + suffix_in;
  }
  fname = baseName + "lh.pial" + suffix;
  if (!FileExists(fname.c_str())) {
    fname = baseName + "lh.pial";
    pialname = "pial";
  }
  if (MRISreadPialCoordinates(fcd->mris_lh, pialname.c_str()) != NO_ERROR) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load lh pial vertices");
  }

  fcd->mris_lh_pial = MRISread(fname.c_str());
  if (fcd->mris_lh_pial == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  fname = baseName + "lh.sphere.d1.left_right" + suffix;
  if (!FileExists(fname.c_str())) {
    fname = baseName + "lh.sphere.d1.left_right";
  }
  fcd->mris_lh_sphere_d1 = MRISread(fname.c_str());
  if (fcd->mris_lh_sphere_d1 == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  exec_progress_callback(1, 12, 0, 1);
  fname = baseName + "rh.white" + suffix;
  if (!FileExists(fname.c_str())) {
    fname = baseName + "rh.white";
  }
  fcd->mris_rh = MRISread(fname.c_str());
  if (fcd->mris_rh == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }
  MRISsaveVertexPositions(fcd->mris_rh, WHITE_VERTICES);

  if (suffix_in) {
    pialname = std::string("pial.") + suffix_in;
  }
  fname = baseName + "rh.pial" + suffix;
  if (!FileExists(fname.c_str())) {
    fname = baseName + "rh.pial";
    pialname = "pial";
  }
  if (MRISreadPialCoordinates(fcd->mris_rh, pialname.c_str()) != NO_ERROR) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load rh pial vertices");
  }

  fcd->mris_rh_pial = MRISread(fname.c_str());
  if (fcd->mris_rh_pial == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  fname = baseName + "rh.sphere.d1.left_right" + suffix;
  if (!FileExists(fname.c_str())) {
    fname = baseName + "rh.sphere.d1.left_right";
  }
  fcd->mris_rh_sphere_d1 = MRISread(fname.c_str());
  if (fcd->mris_rh_sphere_d1 == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  exec_progress_callback(2, 12, 0, 1);
  std::string mriBaseName = std::string(sdir) + "/" + std::string(subject) + "/mri/";
  fname = mriBaseName + "aseg" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = mriBaseName + "aseg.mgz";
  }
  fcd->mri_aseg = MRIread(fname.c_str());
  if (fcd->mri_aseg == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  exec_progress_callback(3, 12, 0, 1);
  fname = mriBaseName + "aparc+aseg" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = mriBaseName + "aparc+aseg.mgz";
  }
  fcd->mri_aparc = MRIread(fname.c_str());
  if (fcd->mri_aparc == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  exec_progress_callback(4, 12, 0, 1);
  fcd->mri_flair = NULL;
  fname = mriBaseName + "flair.reg.norm" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = mriBaseName + "FLAIR" + suffix + ".mgz";
    if (!FileExists(fname.c_str())) {
      fname = mriBaseName + "FLAIRax" + suffix + ".mgz";
      if (!FileExists(fname.c_str())) {
	fname = mriBaseName + "FLAIRcor" + suffix + ".mgz";
        if (!FileExists(fname.c_str())) {
	  fname = " ";
        }
      }
    }
  }
  if ( fname.size() <= 1) {
    fname = mriBaseName + "flair.reg.norm.mgz";
    if (!FileExists(fname.c_str())) {
      fname = mriBaseName + "FLAIR.mgz";
      if (!FileExists(fname.c_str())) {
	fname = mriBaseName + "FLAIRax.mgz";
        if (!FileExists(fname.c_str())) {
	  fname = mriBaseName + "FLAIRcor.mgz";
          if (!FileExists(fname.c_str())) {
	    fname = " ";
          }
        }
      }
    }
  }
  if ( fname.size() > 1) {
    fcd->mri_flair = MRIread(fname.c_str());
    if (fcd->mri_flair == NULL) {
      ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load $s", fname.c_str());
    }
  }

  fcd->mri_t2 = NULL;
  fname = mriBaseName + "T2" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = mriBaseName + "T2ax" + suffix + ".mgz";
    if (!FileExists(fname.c_str())) {
      fname = mriBaseName + "T2cor" + suffix + ".mgz";
      if (!FileExists(fname.c_str())) {
	fname = " ";
      }
    }
  }
  if ( fname.size() <= 1) {
    fname = mriBaseName + "T2.mgz";
    if (!FileExists(fname.c_str())) {
      fname = mriBaseName + "T2ax.mgz";
      if (!FileExists(fname.c_str())) {
	fname = mriBaseName + "T2cor.mgz";
        if (!FileExists(fname.c_str())) {
	  fname = " ";
        }
      }
    }
  }
  if ( fname.size() > 1) {
    fcd->mri_t2 = MRIread(fname.c_str());
    if (fcd->mri_t2 == NULL) {
      ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load $s", fname.c_str());
    }
  }

  exec_progress_callback(5, 12, 0, 1);
  fname = mriBaseName + "norm" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = mriBaseName + "norm.mgz";
  }
  fcd->mri_norm = MRIread(fname.c_str());
  if (fcd->mri_norm == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  fcd->mri_thickness_increase = MRIcloneDifferentType(fcd->mri_aseg, MRI_FLOAT);
  fcd->mri_thickness_decrease = MRIcloneDifferentType(fcd->mri_aseg, MRI_FLOAT);
  fcd->mri_thickness_difference = MRIadd(fcd->mri_thickness_increase, fcd->mri_thickness_decrease, NULL);

  exec_progress_callback(6, 12, 0, 1);
  fname = baseName + "lh.rh.thickness.smooth0" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = baseName + "lh.rh.thickness.smooth0.mgz";
  }
  fcd->rh_thickness_on_lh = MRIread(fname.c_str());
  if (fcd->rh_thickness_on_lh == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  exec_progress_callback(7, 12, 0, 1);
  fname = baseName + "rh.thickness" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = baseName + "rh.thickness.mgz";
  }
  fcd->rh_thickness_on_rh = MRIread(fname.c_str());
  if (fcd->rh_thickness_on_rh == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  exec_progress_callback(8, 12, 0, 1);
  fname = baseName + "lh.thickness" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = baseName + "lh.thickness.mgz";
  }
  fcd->lh_thickness_on_lh = MRIread(fname.c_str());
  if (fcd->lh_thickness_on_lh == NULL) {
    ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
  }

  exec_progress_callback(9, 12, 0, 1);
  fname = baseName + "rh.lh.thickness.smooth0" + suffix + ".mgz";
  if (!FileExists(fname.c_str())) {
    fname = baseName + "rh.lh.thickness.smooth0.mgz";
  }
  if (!FileExists(fname.c_str())) {
    fname = baseName + "rh.lh.thickness" + suffix + ".mgz";
    if (!FileExists(fname.c_str())) {
      fname = baseName + "rh.lh.thickness.mgz";
    }
    if (fcd->lh_thickness_on_rh == NULL) {
      ErrorExit(ERROR_NOFILE, "FCDloadData: couldn't load %s", fname.c_str());
    }
  }
  fcd->lh_thickness_on_rh = MRIread(fname.c_str());

  exec_progress_callback(10, 12, 0, 1);
  mri_int_lh = MRIclone(fcd->mri_norm, NULL);
  mri_int_rh = MRIclone(fcd->mri_norm, NULL);
  mri_interior = MRIclone(fcd->mri_norm, NULL);
  mri_dist = MRIcloneDifferentType(mri_interior, MRI_FLOAT);
  MRISrestoreVertexPositions(fcd->mris_lh, PIAL_VERTICES);
  MRISrestoreVertexPositions(fcd->mris_rh, PIAL_VERTICES);
  MRISfillInterior(fcd->mris_lh, mri_interior->xsize, mri_int_lh);
  MRISfillInterior(fcd->mris_rh, mri_interior->xsize, mri_int_rh);

  exec_progress_callback(11, 12, 0, 1);
  MRIor(mri_int_lh, mri_int_rh, mri_interior, 0);
  MRIfree(&mri_int_lh);
  MRIfree(&mri_int_rh);
  MRIbinarize(mri_interior, mri_interior, 1, 0, 1);
  if (Gdiag & DIAG_WRITE) {
    MRIwrite(mri_interior, "int.mgz");
  }
  MRIdistanceTransform(mri_interior, mri_dist, 1, 2 * MAX_DIST, DTRANS_MODE_SIGNED, NULL);
  if (Gdiag & DIAG_WRITE) {
    MRIwrite(mri_dist, "dist.mgz");
  }
  mri_pvals = build_distance_by_intensity_histo(fcd->mri_norm, mri_dist, fcd->mri_aseg, DIST_SPACING, MAX_DIST);
  exec_progress_callback(12, 12, 0, 1);
  augment_thicknesses(fcd, mri_pvals, 1.5, 5, 1);

  MRISrestoreVertexPositions(fcd->mris_lh, WHITE_VERTICES);
  MRISrestoreVertexPositions(fcd->mris_rh, WHITE_VERTICES);
  MRIfree(&mri_dist);
  MRIfree(&mri_interior);
  MRIfree(&mri_pvals);

  return (fcd);
}

static int sort_labels(FCD_DATA *fcd)
{
  int i;

  qsort(fcd->labels, fcd->nlabels, sizeof(LABEL *), compare_labels);
  for (i = 0; i < fcd->nlabels; i++) {
    strcpy(fcd->label_names[i], fcd->labels[i]->name);
  }
  return (NO_ERROR);
}

static int compare_labels(const void *v1, const void *v2)
{
  LABEL *l1, *l2;

  l1 = *((LABEL **)v1);
  l2 = *((LABEL **)v2);
  if (l1->avg_stat > l2->avg_stat) {
    return (-1);
  }
  if (l1->avg_stat < l2->avg_stat) {
    return (+1);
  }
  return (0);  // equal
}

int FCDcomputeThicknessLabels(FCD_DATA *fcd, double thickness_thresh, double sigma, int size_thresh)
{
  MRI *mri_lh, *mri_rh, *mri_lh_diff, *mri_rh_diff;
  int niter, vno, s;
  MRI_SEGMENTATION *mriseg;

  fcdFreeLabels(fcd);  // free old ones if they exist
  niter = SIGMA_TO_SURFACE_SMOOTH_STEPS(sigma);

  // do LH
  mri_lh = MRIclone(fcd->lh_thickness_on_lh, NULL);
  mri_rh = MRIclone(fcd->lh_thickness_on_lh, NULL);

  exec_progress_callback(1, 8, 0, 1);
  MRISwriteFrameToValues(fcd->mris_lh, fcd->lh_thickness_on_lh, 0);
  MRISaverageVals(fcd->mris_lh, niter);
  MRISreadFrameFromValues(fcd->mris_lh, mri_lh, 0);

  exec_progress_callback(2, 8, 0, 1);
  MRISwriteFrameToValues(fcd->mris_lh, fcd->rh_thickness_on_lh, 0);
  MRISaverageVals(fcd->mris_lh, niter);
  MRISreadFrameFromValues(fcd->mris_lh, mri_rh, 0);
  mri_lh_diff = MRIsubtract(mri_lh, mri_rh, NULL);  // lh minus rh on lh
  MRIfree(&mri_lh);
  MRIfree(&mri_rh);

  // do RH
  mri_lh = MRIclone(fcd->lh_thickness_on_rh, NULL);
  mri_rh = MRIclone(fcd->lh_thickness_on_rh, NULL);

  exec_progress_callback(3, 8, 0, 1);
  MRISwriteFrameToValues(fcd->mris_rh, fcd->lh_thickness_on_rh, 0);
  MRISaverageVals(fcd->mris_rh, niter);
  MRISreadFrameFromValues(fcd->mris_rh, mri_lh, 0);

  exec_progress_callback(4, 8, 0, 1);
  MRISwriteFrameToValues(fcd->mris_rh, fcd->rh_thickness_on_rh, 0);
  MRISaverageVals(fcd->mris_rh, niter);
  MRISreadFrameFromValues(fcd->mris_rh, mri_rh, 0);
  mri_rh_diff = MRIsubtract(mri_rh, mri_lh, NULL);  // lh minus rh on rh
  MRIfree(&mri_lh);
  MRIfree(&mri_rh);

  MRIclear(fcd->mri_thickness_increase);
  MRIclear(fcd->mri_thickness_decrease);
  exec_progress_callback(5, 8, 0, 1);

// process left hemisphere
  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) shared(fcd, mri_lh_diff, Gdiag_no, thickness_thresh) schedule(static, 1)
#endif
  for (vno = 0; vno < fcd->mris_lh->nvertices; vno++) {
    ROMP_PFLB_begin
    double d;
    float val, val2, thickness;
    int base_label;
    VERTEX *v;

    v = &fcd->mris_lh->vertices[vno];
    if (v->ripflag) {
      ROMP_PFLB_continue;
    }
    thickness = MRIgetVoxVal(fcd->lh_thickness_on_lh, vno, 0, 0, 0);
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    val = MRIgetVoxVal(mri_lh_diff, vno, 0, 0, 0);
    if (fabs(val) < thickness_thresh) {
      ROMP_PFLB_continue;
    }

    for (d = 0, base_label = 0; d < thickness; d += 0.25) {
      double xv, yv, zv;
      double xs = v->x + d * v->nx;
      double ys = v->y + d * v->ny;
      double zs = v->z + d * v->nz;
      MRISsurfaceRASToVoxel(fcd->mris_lh, fcd->mri_thickness_increase, xs, ys, zs, &xv, &yv, &zv);
      int xvi = nint(xv);
      int yvi = nint(yv);
      int zvi = nint(zv);
      int label = MRIgetVoxVal(fcd->mri_aparc, xvi, yvi, zvi, 0);
      if (IS_WM(label) == 0 && label >= MIN_CORTICAL_PARCELLATION && label != ctx_lh_unknown) {
        if (label != base_label) {
          if (base_label) {
            break;
          }
        }
        else {
          base_label = label;
        }
        if (val >= 0) {
          val2 = MRIgetVoxVal(fcd->mri_thickness_increase, xvi, yvi, zvi, 0);
          // check another thread already populated this voxel
          if (val > val2) {
            MRIsetVoxVal(fcd->mri_thickness_increase, xvi, yvi, zvi, 0, val);
          }
        }
        else {
          val2 = MRIgetVoxVal(fcd->mri_thickness_decrease, xvi, yvi, zvi, 0);
          // check if another thread already populated this voxel
          if (val < val2) {
            MRIsetVoxVal(fcd->mri_thickness_decrease, xvi, yvi, zvi, 0, val);
          }
        }
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end

  exec_progress_callback(6, 8, 0, 1);

// now do right hemisphere

  ROMP_PF_begin
#ifdef HAVE_OPENMP
  #pragma omp parallel for if_ROMP(experimental) shared(fcd, mri_rh_diff, Gdiag_no, thickness_thresh) schedule(static, 1)
#endif
  for (vno = 0; vno < fcd->mris_rh->nvertices; vno++) {
    ROMP_PFLB_begin
    double d;
    float val, val2, thickness;
    int base_label;
    VERTEX *v;

    v = &fcd->mris_rh->vertices[vno];
    if (v->ripflag) {
      ROMP_PFLB_continue;
    }
    if (vno == Gdiag_no) {
      DiagBreak();
    }
    val = MRIgetVoxVal(mri_rh_diff, vno, 0, 0, 0);
    if (fabs(val) < thickness_thresh) {
      ROMP_PFLB_continue;
    }
    thickness = MRIgetVoxVal(fcd->rh_thickness_on_rh, vno, 0, 0, 0);

    for (d = 0, base_label = 0; d < thickness; d += 0.25) {
      double xv, yv, zv;
      double xs = v->x + d * v->nx;
      double ys = v->y + d * v->ny;
      double zs = v->z + d * v->nz;
      MRISsurfaceRASToVoxel(fcd->mris_rh, fcd->mri_thickness_increase, xs, ys, zs, &xv, &yv, &zv);
      int xvi = nint(xv);
      int yvi = nint(yv);
      int zvi = nint(zv);
      int label = MRIgetVoxVal(fcd->mri_aparc, xvi, yvi, zvi, 0);
      if (IS_WM(label) == 0 && label >= MIN_CORTICAL_PARCELLATION && label != ctx_rh_unknown) {
        if (label != base_label) {
          if (base_label) {
            break;
          }
        }
        else {
          base_label = label;
        }

        if (val >= 0) {
          val2 = MRIgetVoxVal(fcd->mri_thickness_increase, xvi, yvi, zvi, 0);
          if (val > val2) {
            MRIsetVoxVal(fcd->mri_thickness_increase, xvi, yvi, zvi, 0, val);
          }
        }
        else {
          val2 = MRIgetVoxVal(fcd->mri_thickness_decrease, xvi, yvi, zvi, 0);
          if (val < val2) {
            MRIsetVoxVal(fcd->mri_thickness_decrease, xvi, yvi, zvi, 0, val);
          }
        }
      }
    }
    ROMP_PFLB_end
  }
  ROMP_PF_end
  
  exec_progress_callback(7, 8, 0, 1);
  mriseg = MRIsegment(fcd->mri_thickness_increase, thickness_thresh, 1e10);
  MRIeraseSmallSegments(mriseg, fcd->mri_thickness_increase, thickness_thresh);
  MRIsegmentFree(&mriseg);
  MRIclose(fcd->mri_thickness_increase, fcd->mri_thickness_increase);
  mriseg = MRIsegment(fcd->mri_thickness_increase, thickness_thresh, 1e10);
  MRIremoveSmallSegments(mriseg, size_thresh);
  printf(
      "segmenting volume at threshold %2.1f with %d "
      "smoothing iters yields %d segments\n",
      thickness_thresh,
      niter,
      mriseg->nsegments);
  fflush(stdout);

  exec_progress_callback(8, 8, 0, 1);
  fcd->nlabels = mriseg->nsegments;
  for (s = 0; s < mriseg->nsegments; s++) {
    int label;

    fcd->labels[s] = MRIsegmentToLabel(mriseg, fcd->mri_thickness_increase, s);
    label = most_frequent_label(fcd->mri_aparc, &mriseg->segments[s]);
    strcpy(fcd->labels[s]->name, cma_label_to_name(label));
  }
  sort_labels(fcd);

  MRIadd(fcd->mri_thickness_increase, fcd->mri_thickness_decrease, fcd->mri_thickness_difference);

  for (s = 0; s < mriseg->nsegments; s++) {
    printf("%s: %2.3fmm\n", fcd->label_names[s], fcd->labels[s]->avg_stat);
    fflush(stdout);
  }
  MRIfree(&mri_lh_diff);
  MRIfree(&mri_rh_diff);
  MRIsegmentFree(&mriseg);

  return (fcd->nlabels);
}

int FCDwriteLabels(FCD_DATA *fcd, char *dir)
{
  int s;
  char label_name[STRLEN];

  for (s = 0; s < fcd->nlabels; s++) {
    if (fcd->labels[s]) {
      sprintf(label_name, "%s/fcd_%02d_%s", dir, s, fcd->labels[s]->name);
      LabelWrite(fcd->labels[s], label_name);
    }
  }

  printf("wrote FCD labels to %s\n", dir);

  return (NO_ERROR);
}

static int fcdFreeLabels(FCD_DATA *fcd)
{
  int s;

  for (s = 0; s < fcd->nlabels; s++)
    if (fcd->labels[s]) {
      LabelFree(&fcd->labels[s]);
    }
  fcd->nlabels = 0;
  return (NO_ERROR);
}

int FCDfree(FCD_DATA **pfcd)
{
  FCD_DATA *fcd;

  fcd = *pfcd;
  *pfcd = NULL;
  if (fcd->mris_lh) {
    MRISfree(&fcd->mris_lh);
  }
  if (fcd->mris_rh) {
    MRISfree(&fcd->mris_rh);
  }
  if (fcd->mris_lh_pial) {
    MRISfree(&fcd->mris_lh_pial);
  }
  if (fcd->mris_rh_pial) {
    MRISfree(&fcd->mris_rh_pial);
  }

  if (fcd->mris_lh_sphere_d1) {
    MRISfree(&fcd->mris_lh_sphere_d1);
  }
  if (fcd->mris_rh_sphere_d1) {
    MRISfree(&fcd->mris_rh_sphere_d1);
  }

  if (fcd->mri_aseg) {
    MRIfree(&fcd->mri_aseg);
  }
  if (fcd->mri_aparc) {
    MRIfree(&fcd->mri_aparc);
  }
  if (fcd->mri_norm) {
    MRIfree(&fcd->mri_norm);
  }
  if (fcd->mri_flair) {
    MRIfree(&fcd->mri_flair);
  }
  if (fcd->mri_thickness_increase) {
    MRIfree(&fcd->mri_thickness_increase);
  }
  if (fcd->mri_thickness_decrease) {
    MRIfree(&fcd->mri_thickness_decrease);
  }
  if (fcd->mri_thickness_difference) {
    MRIfree(&fcd->mri_thickness_difference);
  }
  if (fcd->lh_thickness_on_lh) {
    MRIfree(&fcd->lh_thickness_on_lh);
  }
  if (fcd->lh_thickness_on_rh) {
    MRIfree(&fcd->lh_thickness_on_rh);
  }
  if (fcd->rh_thickness_on_lh) {
    MRIfree(&fcd->rh_thickness_on_lh);
  }
  if (fcd->rh_thickness_on_rh) {
    MRIfree(&fcd->rh_thickness_on_rh);
  }

  fcdFreeLabels(fcd);
  free(fcd);
  return (NO_ERROR);
}

static MRI *build_distance_by_intensity_histo(
    MRI *mri_norm, MRI *mri_dist, MRI *mri_aseg, double dist_spacing, double max_dist)
{
  HISTOGRAM2D *h_dist_by_int;
  int x, y, z, b1, b2, label;
  float val, dist;
  double total, unlikely, pval;
  MRI *mri_pvals;

  h_dist_by_int = HISTO2Dalloc((int)ceil(max_dist / dist_spacing), 256);
  HISTO2Dinit(h_dist_by_int, h_dist_by_int->nbins1, h_dist_by_int->nbins2, 0, MAX_DIST, 0, 255);
  for (x = 0; x < mri_dist->width; x++)
    for (y = 0; y < mri_dist->height; y++)
      for (z = 0; z < mri_dist->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        dist = MRIgetVoxVal(mri_dist, x, y, z, 0);
        if (dist < 0 || dist > MAX_DIST)  // in interior or too far away
        {
          continue;
        }
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0);
        if (IS_WHITE_CLASS(label) == 0 && IS_CORTEX(label) == 0 && label < MIN_CORTICAL_PARCELLATION) {
          continue;
        }
        val = MRIgetVoxVal(mri_norm, x, y, z, 0);
        HISTO2DaddSample(h_dist_by_int, dist, val, 0, max_dist, 0, 255);
      }

  // normalize the counts for each distance
  for (b1 = 0; b1 < h_dist_by_int->nbins1; b1++) {
    for (total = 0.0, b2 = 0; b2 < h_dist_by_int->nbins2; b2++) {
      total += h_dist_by_int->counts[b1][b2];
    }

    if (total > 0) {
      unlikely = 1.0 / (10 * total);
      for (b2 = 0; b2 < h_dist_by_int->nbins2; b2++) {
        h_dist_by_int->counts[b1][b2] /= total;
        if (DZERO(h_dist_by_int->counts[b1][b2])) {
          h_dist_by_int->counts[b1][b2] = unlikely;
        }
      }
    }
  }

  mri_pvals = MRIclone(mri_dist, NULL);
  for (x = 0; x < mri_dist->width; x++)
    for (y = 0; y < mri_dist->height; y++)
      for (z = 0; z < mri_dist->depth; z++) {
        if (x == Gx && y == Gy && z == Gz) {
          DiagBreak();
        }
        dist = MRIgetVoxVal(mri_dist, x, y, z, 0);
        if (dist < 0 || dist > MAX_DIST)  // in interior
        {
          continue;
        }
        label = MRIgetVoxVal(mri_aseg, x, y, z, 0);
        if (IS_WHITE_CLASS(label) == 0 && IS_CORTEX(label) == 0 && label < MIN_CORTICAL_PARCELLATION) {
          continue;
        }
        val = MRIgetVoxVal(mri_norm, x, y, z, 0);
        pval = HISTO2DgetCount(h_dist_by_int, dist, val);
        if (pval > 0) {
          pval = -log10(pval);
        }
        else {
          pval = -10000;
        }
        MRIsetVoxVal(mri_pvals, x, y, z, 0, pval);
      }

  HISTO2Dfree(&h_dist_by_int);
  return (mri_pvals);
}

static int augment_thicknesses(FCD_DATA *fcd, MRI *mri_pvals, double min_dist, double max_dist, double thresh)
{
  int h, vno;
  VERTEX *v;
  MRI_SURFACE *mris;
  MRI *mri_thickness;
  double nx, ny, nz, x0, y0, z0, d, x, y, z, val;
  MRI_SEGMENTATION *mriseg;

  mriseg = MRIsegment(mri_pvals, thresh, 1e10);
  MRIeraseSmallSegments(mriseg, mri_pvals, 20);
  MRIremoveSmallSegments(mriseg, 100);
  if (Gdiag & DIAG_WRITE) {
    MRIwrite(mri_pvals, "pvals.mgz");
  }
  MRIsegmentFree(&mriseg);

  for (h = 0; h <= 1; h++)  // do each hemi
  {
    if (h == 0)  // left hemi
    {
      mri_thickness = fcd->lh_thickness_on_lh;
      mris = fcd->mris_lh;
    }
    else  // right hemi
    {
      mri_thickness = fcd->rh_thickness_on_rh;
      mris = fcd->mris_rh;
    }
    for (vno = 0; vno < mris->nvertices; vno++) {
      if (vno == Gdiag_no) {
        DiagBreak();
      }
      v = &mris->vertices[vno];
      if (v->ripflag) {
        continue;
      }
      MRISvertexNormalInVoxelCoords(mris, mri_pvals, vno, &nx, &ny, &nz);
      MRISvertexToVoxel(mris, v, mri_pvals, &x0, &y0, &z0);

      for (d = 0; d <= max_dist; d += 0.5) {
        x = x0 + d * nx;
        y = y0 + d * ny;
        z = z0 + d * nz;
        MRIsampleVolume(mri_pvals, x, y, z, &val);
        if (val < thresh) {
          break;
        }
      }

      if (d > min_dist)  // a string of unlikely values
      {
        val = MRIgetVoxVal(mri_thickness, vno, 0, 0, 0);
        MRIsetVoxVal(mri_thickness, vno, 0, 0, 0, val + d);
      }
    }
  }
  return (NO_ERROR);
}
