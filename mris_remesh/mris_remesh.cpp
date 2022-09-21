#include "argparse.h"
#include "mrisurf.h"
#include "remesher.h"
#include "surfgrad.h"
#include "mris_remesh.help.xml.h"


int main(int argc, char **argv) 
{
  ArgumentParser parser;
  parser.addHelp(mris_remesh_help_xml, mris_remesh_help_xml_len);
  // required
  parser.addArgument("-i", "--input",  1, String, true);
  parser.addArgument("-o", "--output", 1, String, true);
  // one of these is required
  parser.addArgument("--remesh");
  parser.addArgument("--nvert", 1, Int);
  parser.addArgument("--edge-len", 1, Float);
  parser.addArgument("--desired-face-area", 1, Float);
  // optional
  parser.addArgument("--iters", 1, Int);
  parser.parse(argc, argv);

  // read input surface
  std::string inputname = parser.retrieve<std::string>("input");
  MRIS *surf = MRISread(inputname.c_str());
  if (!surf) fs::fatal() << "could not read input surface " << inputname;

  // number of iterations
  int iters = parser.exists("iters") ? parser.retrieve<int>("iters") : 5;
  std::cout << "iters = " << iters << std::endl;

  MRIScomputeMetricProperties(surf);

  // init the remesher
  Remesher remesher = Remesher(surf);

  // quick sanity check to make sure only one remesh method was specified
  int numtargets = int(parser.exists("remesh")) +
                   int(parser.exists("nvert")) +
                   int(parser.exists("edge-len")) +
                   int(parser.exists("desired-face-area"));
  if (numtargets > 1) fs::fatal() << "must only specify one remeshing target";

  if (parser.exists("nvert")) {
    // remesh to number of vertices
    int nverts = parser.retrieve<int>("nvert");
    std::cout << "target vertices = " << nverts << std::endl;
    remesher.remeshBKV(iters, nverts);
  }
  else if (parser.exists("remesh")) {
    // remesh without changing nvertices - just modify the quality
    std::cout << "standard remeshing without target" << std::endl;
    remesher.remeshBK(iters);
  }
  else if (parser.exists("edge-len")) {
    // remesh to target edge length
    float edgelength = parser.retrieve<float>("edge-len");
    std::cout << "target edge length = " << edgelength << std::endl;
    remesher.remeshBK(iters, edgelength);
  }
  else if (parser.exists("desired-face-area")) {
    // remesh to target face area
    float  desiredFaceArea = parser.retrieve<float>("desired-face-area");
    double avgfacearea = surf->total_area / surf->nfaces;
    double decimationLevel = avgfacearea / desiredFaceArea;
    int nverts = round(surf->nvertices * decimationLevel);
    std::cout << "target face area = " << desiredFaceArea << std::endl;
    std::cout << "average source face area = " << avgfacearea << std::endl;
    std::cout << "decimation level = " << decimationLevel << std::endl;
    std::cout << "target vertices = " << nverts << std::endl;
    remesher.remeshBKV(iters, nverts);
  }
  else {
    fs::fatal() << "must specify target edge length, number of vertices, or face area";
  }

  // convert back to MRIS
  MRIS *remeshed = remesher.toSurface();
  MRIScopyMetadata(surf, remeshed);

  // remove the sea of intersections
  std::cout << "Removing intersections" << std::endl;
  MRISremoveIntersections(remeshed,0);

  // print stats
  double diff = (double)remeshed->nvertices / (double)surf->nvertices;
  printf("Remeshed surface quality stats nv0 = %d  nv = %d  %g\n", surf->nvertices, remeshed->nvertices, diff);
  MRISedges(remeshed);
  MRIScorners(remeshed);
  MRISfaceMetric(remeshed, 0);
  MRISedgeMetric(remeshed, 0);
  MRIScornerMetric(remeshed, 0);
  MRISprettyPrintSurfQualityStats(stdout, remeshed);

  // save
  std::string outputname = parser.retrieve<std::string>("output");
  MRISwrite(remeshed, outputname.c_str());

  MRISfree(&surf);
  MRISfree(&remeshed);
  printf("mris_remesh done\n");

  exit(0);
}
