#include "argparse.h"
#include "mrisurf.h"
#include "remesher.h"
#include "surfgrad.h"


int main(int argc, const char **argv) 
{
  ArgumentParser parser;
  parser.addArgument("-i", "--input", 1, String, true);
  parser.addArgument("-o", "--output", 1, String, true);
  parser.addArgument("--iters", 1, Int);
  parser.addArgument("--edge-len", 1, Float);
  parser.parse(argc, argv);

  std::string inputname = parser.retrieve<std::string>("input");
  MRIS *surf = MRISread(inputname.c_str());
  if (!surf) fs::fatal() << "could not read input surface " << inputname;

  int iters = parser.exists("iters") ? parser.retrieve<int>("iters") : 5;
  float length = parser.exists("edge-len") ? parser.retrieve<float>("edge-len") : 0.8;

  Remesher remesher = Remesher(surf);
  remesher.remeshBK(iters, length, false);

  MRIS *remeshed = remesher.toSurface();
  MRIScopyMetadata(surf, remeshed);

  std::cout << "Removing intersections" << std::endl;
  MRISremoveIntersections(remeshed);

  double diff = (double)remeshed->nvertices / (double)surf->nvertices;
  printf("Remeshed surface quality stats nv0 = %d  nv = %d  %g\n", surf->nvertices, remeshed->nvertices, diff);
  MRISedges(remeshed);
  MRIScorners(remeshed);
  MRISfaceMetric(remeshed, 0);
  MRISedgeMetric(remeshed, 0);
  MRIScornerMetric(remeshed, 0);
  MRISprettyPrintSurfQualityStats(stdout, remeshed);

  std::string outputname = parser.retrieve<std::string>("output");
  MRISwrite(remeshed, outputname.c_str());

  MRISfree(&surf);
  MRISfree(&remeshed);

  exit(0);
}
