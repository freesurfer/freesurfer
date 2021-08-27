#include <string>
#include <iostream>
#include <fstream>

#include "mris_reposition_surface.help.xml.h"
#include "argparse.h"
#include "pointset.h"
#include "log.h"

 
#include "utils.h"
#include "mri.h"
#include "label.h"
#include "mrisurf.h"
#include "mri_circulars.h"
#include "json.h"

using json = nlohmann::json;

int main(int argc, char **argv) 
{

  // --- setup ----

  ArgumentParser parser;
  // required
  parser.addArgument("-s", "--surf",    1, String, true);
  parser.addArgument("-v", "--volume",    1, String, true);
  parser.addArgument("-p", "--points", 1, String, true);
  parser.addArgument("-o", "--out",     1, String, true);
  // optional
  parser.addArgument("-z", "--size",   1, Int);
  parser.addArgument("-g", "--sigma", 1, Float);
  parser.addArgument("-i", "--iterations", 1, Int);

  // help text
  parser.addHelp(mris_reposition_surface_help_xml, mris_reposition_surface_help_xml_len);
  parser.parse(argc, argv);

  // load surface
  std::string surfpath = parser.retrieve<std::string>("surf");
  std::cout << "Reading in surface " << surfpath << std::endl;
  MRIS *surf = MRISread(surfpath.c_str());
  if (!surf) fs::fatal() << "could not read surface";

  // load volume
  std::string mripath = parser.retrieve<std::string>("volume");
  std::cout << "Reading in volume " << mripath << std::endl;
  MRI* mri = MRIread(mripath.c_str());
  if (!mri) fs::fatal() << "could not read volume";

  // load pointset
  std::string ptspath = parser.retrieve<std::string>("points");
  std::cout << "Reading in control points " << ptspath << std::endl;
  std::ifstream istream(ptspath);
  json js;
  istream >> js;
  json js_pts = js["points"];
  if (js_pts.size() < 1) fs::fatal() << "could not read control points";

  bool bTkReg = (js["vox2ras"].get<std::string>() == std::string("tkreg"));
  PointSet points;
  for (int i = 0; i < js_pts.size(); i++) {
      PointSet::Point p;
      p.x = js_pts[i]["coordinates"]["x"].get<float>();
      p.y = js_pts[i]["coordinates"]["y"].get<float>();
      p.z = js_pts[i]["coordinates"]["z"].get<float>();
      // if the surface contains volume geometry, convert the point to scanner RAS
      if (!bTkReg) {
        double sx, sy, sz;
        MRIRASToSurfaceRAS(mri, p.x, p.y, p.z, &sx, &sy, &sz);
        p.x = sx; p.y = sy; p.z = sz;
      }
      points.add(p);
  }

  // size input
  int nsize = 1;
  if (parser.exists("size")) {
    nsize = parser.retrieve<int>("size");
    if (nsize < 1) fs::fatal() << "size must be greater than 0";
  }

  // sigma input
  float sigma = 2;
  if (parser.exists("sigma")) {
    sigma = parser.retrieve<float>("sigma");
    if (sigma <= 0) fs::fatal() << "sigma must be greater than 0";
  }

  // number of iterations
  float niterations = 1;
  if (parser.exists("iterations")) {
    niterations = parser.retrieve<int>("iterations");
    if (niterations < 1) niterations = 1;
  }

  for (int i = 0; i < niterations; i++) {
    double max_spacing;
    int max_vno;
    float distance;
    MRIScomputeVertexSpacingStats(surf, NULL, NULL, &max_spacing, NULL, &max_vno, CURRENT_VERTICES);
    MRIS_HASH_TABLE* hash = MHTcreateVertexTable_Resolution(surf, CURRENT_VERTICES, max_spacing);
    // cycle through each points
    for (auto it = points.begin(); it != points.end(); it++) {
      PointSet::Point p = *it;
      int vno = MHTfindClosestVertexNoXYZ(hash, surf, p.x, p.y, p.z, &distance);
      if (vno < 0){
	printf("Failed to find closest vertex in hash, using brute force\n");
	vno = MRISfindClosestVertex(surf, p.x, p.y, p.z, &distance, CURRENT_VERTICES);
	//fs::fatal() << "failed to find closest vertex";
      }
      MRISrepositionSurfaceToCoordinate(surf, mri, vno, p.x, p.y, p.z, nsize, sigma, 0);
    }
    MRIScomputeNormals(surf);
    MHTfree(&hash);
  }

  // ---- output
  if (MRISwrite(surf, parser.retrieve<std::string>("out").c_str()) != 0)
    fs::fatal() << "failed to write surface";

  // ---- cleanup ----
  MRISfree(&surf);
  MRIfree(&mri);

  std::cout << "#VMPC# mris_reposition_surface " << GetVmPeak() << std::endl;
  std::cout << "mris_reposition_surface done" << std::endl;

  exit(0);
}
