#include <string>
#include <iostream>
#include <fstream>

#include "argparse.hpp"
#include "pointset.hpp"

extern "C" {
#include "utils.h"
#include "mri.h"
#include "label.h"
#include "mrisurf.h"
#include "mri_circulars.h"
}

static PointSet::Point sras2ras(MRIS* surf, PointSet::Point point);

int main(int argc, const char **argv) 
{
  LABEL *label;
  int v;

  // --- setup ----

  ArgumentParser parser;
  parser.addArgument("-s", "--surf",    1, String, true);
  parser.addArgument("-d", "--defects", 1, String, true);
  parser.addArgument("-o", "--out",     1, String, true);
  parser.addArgument("-l", "--label",     1, String, true);
  if(argc < 3){
    printf("mris_defects_pointset \n");
    printf("  -s surface (eg, lh.orig.nofix) \n");
    printf("  -d defectseg (eg, lh.defect_labels) \n");
    printf("  -l label : restrict pointset to this label (optional, must be in eg, orig.nofix space) \n");
    printf("\n");
    exit(1);
  }
  parser.parse(argc, argv);

  // load surface
  std::string surfpath = parser.retrieve<std::string>("surf");
  MRIS *surf = MRISread(surfpath.c_str());
  printf("Reading in surface %s\n",surfpath.c_str());
  if (!surf) exit(1);

  // load defect overlay
  std::string defectpath = parser.retrieve<std::string>("defects");
  printf("Reading in defect segmentation %s\n",defectpath.c_str());
  MRI *overlay = MRIread(defectpath.c_str());
  if (!overlay) exit(1);

  if (overlay->width != surf->nvertices) {
    std::cout << "error: defect overlay (" << overlay->width << " points) "
              << "does not match surface (" << surf->nvertices << " vertices)"
              << std::endl;
    exit(1);
  }

  //mris_apply_reg --src-label ../label/lh.cortex.label --streg lh.orig lh.orig.nofix --trg ../label/lh.cortex.nofix.label
  std::string labelpath = parser.retrieve<std::string>("label");
  if(labelpath.length() > 0){
    printf("Reading in label %s\n",labelpath.c_str());
    label = LabelRead(NULL, labelpath.c_str());
    if(label==NULL){
      std::cout << "ERROR: loading " << labelpath << std::endl;
      exit(1);
    }
    // Set values outside of the label to 0
    MRI *tmpmri;
    tmpmri = MRISlabel2Mask(surf, label, NULL);
    for(v=0; v < surf->nvertices; v++){
      if(MRIgetVoxVal(tmpmri,v,0,0,0)) continue;
      MRIsetVoxVal(overlay, v, 0, 0, 0, 0);      
    }
    MRIfree(&tmpmri);
    LabelFree(&label);
  }

  // ---- get number of defects (max value in overlay) ----

  int ndefects = 0;
  for (v = 0 ; v < surf->nvertices ; v++) {
    int value = MRIgetVoxVal(overlay, v, 0, 0, 0);
    if (value > ndefects) ndefects = value;
  }

  // ---- compute centroids of each defect ----

  // create a vector of vectors that contain the vertices for each defect
  std::vector<std::vector<int>> defects(ndefects);
  for (int v = 0 ; v < surf->nvertices ; v++) {
    int defect_idx = MRIgetVoxVal(overlay, v, 0, 0, 0) - 1;
    if (defect_idx >= 0) defects[defect_idx].push_back(v);
  }

  // cycle through each defect and find its center position
  PointSet centroids;
  for (auto &defect : defects) {
    if (!defect.empty()) {
      // average the points in the defect
      PointSet::Point centroid;
      int npoints = 0;
      for (auto &vnum : defect) {
        centroid.x += surf->vertices[vnum].x;
        centroid.y += surf->vertices[vnum].y;
        centroid.z += surf->vertices[vnum].z;
        npoints++;
      }
      centroid.x /= npoints;
      centroid.y /= npoints;
      centroid.z /= npoints;
      // if the surface contains volume geometry, convert the point to
      // scanner RAS coordinates
      if (surf->vg.valid) centroid = sras2ras(surf, centroid);
      centroids.add(centroid);
    }
  }


  // ---- write the pointset file ----

  if (surf->vg.valid) {
    centroids.vox2ras = "scanner_ras";
  } else {
    centroids.vox2ras = "tkreg";
  }

  centroids.save(parser.retrieve<std::string>("out"));


  // ---- cleanup ----

  MRIfree(&overlay);
  MRISfree(&surf);

  printf("#VMPC# mris_defects_pointset %d\n",GetVmPeak());
  printf("mris_defects_pointset done\n");

  exit(0);

}


// returns a point transformed from surface RAS to scanner RAS coordinates - ideally
// this function should be moved to mrisurf.c or somewhere else more appropriate 
static PointSet::Point sras2ras(MRIS* surf, PointSet::Point point)
{
  static MATRIX *sras2ras_matrix;
  static VECTOR *v1, *v2;

  // allocate the matrix and vectors (this will only be done once in the entire program)
  if (!sras2ras_matrix) {
    // quick way to get the sras2ras matrix by using RASFromSurfaceRAS_ with a tmp volume
    MRI* tmp = MRIallocHeader(surf->vg.width, surf->vg.height, surf->vg.depth, MRI_UCHAR, 1);
    MRIcopyVolGeomToMRI(tmp, &surf->vg);
    sras2ras_matrix = RASFromSurfaceRAS_(tmp);
    MRIfree(&tmp);
    // allocate input and output vectors
    v1 = VectorAlloc(4, MATRIX_REAL);
    v2 = VectorAlloc(4, MATRIX_REAL);
    VECTOR_ELT(v1, 4) = 1.0;
    VECTOR_ELT(v2, 4) = 1.0;
  }

  V3_X(v1) = point.x; V3_Y(v1) = point.y; V3_Z(v1) = point.z;
  MatrixMultiply(sras2ras_matrix, v1, v2);

  return PointSet::Point(V3_X(v2), V3_Y(v2), V3_Z(v2), point.value);
}
