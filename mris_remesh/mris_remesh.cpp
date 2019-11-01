#include "argparse.h"
#include "mrisurf.h"
#include "remesher.h"
#include "surfgrad.h"


int main(int argc, const char **argv) 
{
  int nv0;

  ArgumentParser parser;
  parser.addArgument("-i", "--input", 1, String, true);
  parser.addArgument("-o", "--output", 1, String, true);
  parser.addArgument("--iters", 1, Int);
  parser.addArgument("--edge-len", 1, Float);
  parser.parse(argc, argv);

  std::string inputname = parser.retrieve<std::string>("input");
  MRIS *surf = MRISread(inputname.c_str());
  if (!surf) fs::fatal() << "could not read input surface " << inputname;
  nv0 = surf->nvertices;

  int iters = parser.exists("iters") ? parser.retrieve<int>("iters") : 5;
  float length = parser.exists("edge-len") ? parser.retrieve<float>("edge-len") : 0.8;

  Remesher remesher = Remesher(surf);
  remesher.remeshBK(iters, length, false);
  remesher.updateSurface(surf);

  printf("Computing metric properties\n"); fflush(stdout);
  MRIScomputeMetricProperties(surf);

  printf("Removing intersectios\n");
  MRISremoveIntersections(surf);

  printf("Remeshed surface quality stats nv0 = %d  nv = %d  %g\n",nv0,surf->nvertices,(double)surf->nvertices/nv0);
  MRISedges(surf);
  MRIScorners(surf);
  MRISfaceMetric(surf,0);
  MRISedgeMetric(surf,0);
  MRIScornerMetric(surf,0);
  MRISprettyPrintSurfQualityStats(stdout, surf);

  std::string outputname = parser.retrieve<std::string>("output");
  MRISwrite(surf, outputname.c_str());

  MRISfree(&surf);
  surf = MRISread(outputname.c_str());

  exit(0);
}
