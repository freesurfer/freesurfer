#include "argparse.h"
#include "mrisurf.h"
#include "remesher.h"


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
  remesher.updateSurface(surf);

  std::string outputname = parser.retrieve<std::string>("output");
  MRISwrite(surf, outputname.c_str());
  MRISfree(&surf);
}
