#include <string>
#include <iostream>
#include <fstream>

#include "mris_defects_pointset.help.xml.h"
#include "argparse.h"
#include "pointset.h"
#include "log.h"

 
#include "utils.h"
#include "mri.h"
#include "label.h"
#include "mrisurf.h"
#include "mri_circulars.h"



static fsPointSet::Point sras2ras(MRIS* surf, fsPointSet::Point point);

int main(int argc, char **argv) 
{

  // --- setup ----

  ArgumentParser parser;
  // required
  parser.addArgument("-s", "--surf",    1, String, true);
  parser.addArgument("-d", "--defects", 1, String, true);
  parser.addArgument("-o", "--out",     1, String, true);
  // optional
  parser.addArgument("-l", "--label",   1, String);
  parser.addArgument("-c", "--control", 0);
  // help text
  parser.addHelp(mris_defects_pointset_help_xml, mris_defects_pointset_help_xml_len);
  parser.parse(argc, argv);

  // load surface
  std::string surfpath = parser.retrieve<std::string>("surf");
  std::cout << "Reading in surface " << surfpath << std::endl;
  MRIS *surf = MRISread(surfpath.c_str());
  if (!surf) fs::fatal() << "could not read surface";

  // load defect overlay
  std::string defectpath = parser.retrieve<std::string>("defects");
  std::cout << "Reading in defect segmentation " << defectpath << std::endl;
  MRI *overlay = MRIread(defectpath.c_str());
  if (!overlay) fs::fatal() << "could not read defect segmentation";

  if (overlay->width != surf->nvertices) {
    fs::fatal() << "error: defect overlay (" << overlay->width << " points) "
                << "does not match surface (" << surf->nvertices << " vertices)";
  }

  // ---- optionally mask the defect segmentation to a label ----

  // the provided label must be in the correct space (eg. orig.nofix space):
  // mris_apply_reg --src-label ../label/lh.cortex.label --streg lh.orig lh.orig.nofix --trg ../label/lh.cortex.nofix.label

  if (parser.exists("label")) {
    std::string labelpath = parser.retrieve<std::string>("label");
    std::cout << "Reading in label " << labelpath << std::endl;
    LABEL *label = LabelRead(NULL, labelpath.c_str());
    if (!label) fs::fatal() << "could not read label";
    // set values outside of the label to 0
    MRI *tmp = MRISlabel2Mask(surf, label, NULL);
    for (int v = 0 ; v < surf->nvertices ; v++) {
      if (!MRIgetVoxVal(tmp, v, 0, 0, 0)) MRIsetVoxVal(overlay, v, 0, 0, 0, 0);
    }
    MRIfree(&tmp);
    LabelFree(&label);
  }

  // ---- get number of defects (max value in overlay) ----

  int ndefects = 0;
  for (int v = 0 ; v < surf->nvertices ; v++) {
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
  fsPointSet centroids;
  for (auto &defect : defects) {
    if (!defect.empty()) {
      // average the points in the defect
      fsPointSet::Point centroid;
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
      // if the surface contains volume geometry, convert the point to scanner RAS
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

  std::string outfile = parser.retrieve<std::string>("out");
  if (parser.exists("control")) {
    centroids.save_as_ctrlpoint(outfile);
  } else {
    centroids.save(outfile);
  }

  // ---- cleanup ----

  MRIfree(&overlay);
  MRISfree(&surf);

  std::cout << "#VMPC# mris_defects_pointset " << GetVmPeak() << std::endl;
  std::cout << "mris_defects_pointset done" << std::endl;

  exit(0);

}


// returns a point transformed from surface RAS to scanner RAS coordinates - ideally
// this function should be moved to mrisurf.c or somewhere else more appropriate 
static fsPointSet::Point sras2ras(MRIS* surf, fsPointSet::Point point)
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

  return fsPointSet::Point(V3_X(v2), V3_Y(v2), V3_Z(v2), point.value);
}
