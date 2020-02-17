// A very simple tool that computes a variety of brain volume statistics for a given subject.
// The calculated values are cached in a stats file so they can be easily referenced by
// mri_segstats and mris_anatomical_stats later on.

#include <cstdlib>
#include <string>
#include <vector>

#include "argparse.h"
#include "cma.h"


int main(int argc, char **argv) 
{
  // parse arguments
  ArgumentParser parser;
  parser.addArgument("subject");
  parser.addArgument("--sd", 1);
  parser.parse(argc, argv);

  std::string subject = parser.retrieve<std::string>("subject");
  std::string subjdir = parser.retrieve<std::string>("sd");

  // get subjects directory from flag or env var
  if (subjdir.empty()) {
    const char* sd = std::getenv("SUBJECTS_DIR");
    if (!sd) fs::fatal() << "must set SUBJECTS_DIR or specify path with the --sd flag";
    subjdir = std::string(sd);
  }

  // compute and write the vol stats to a cache file at subjects/stats/brainvol.stats
  std::vector<double> stats = ComputeBrainVolumeStats(subject, subjdir);
  CacheBrainVolumeStats(stats, subject, subjdir);

  return 0;
}
