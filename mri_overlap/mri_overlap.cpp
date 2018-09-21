#include <iostream>
#include <string>
#include <iomanip>
#include <map>

#include "mri_overlap.help.xml.h"
#include "argparse.hpp"
#include "json.hpp"
#include "log.hpp"

extern "C" {
#include "mri.h"
}


int main(int argc, const char **argv) 
{
  // parse arguments
  ArgumentParser parser;
  parser.addHelp(mri_overlap_help_xml, mri_overlap_help_xml_len);
  parser.addArgument("seg1");
  parser.addArgument("seg2");
  parser.addArgument("-o", "--out");
  parser.addArgument("-m", "--measures", '+');
  parser.addArgument("-l", "--labels", '+', Int);
  parser.addArgument("-n", "--names", '+');
  parser.addArgument("-q", "--quiet");
  parser.parse(argc, argv);

  // load inputs
  std::string seg1_fname = parser.retrieve<std::string>("seg1");
  std::string seg2_fname = parser.retrieve<std::string>("seg2");
  MRI *seg1 = MRIread(seg1_fname.c_str());
  MRI *seg2 = MRIread(seg2_fname.c_str());
  if (!seg1 || !seg2) fs_fatal(1) << "could not read input volume";

  // check input dimensions
  if ((seg1->width  != seg2->width)  ||
      (seg1->height != seg2->height) ||
      (seg1->depth  != seg2->depth)) {
    fs_fatal(1) << "input volumes must have matching dimensions";
  }

  // check input frames
  if (seg1->nframes != seg2->nframes) {
    fs_fatal(1) << "input volumes must have the same number of frames";
  }

  // --- setup ----

  struct OverlapMetrics {
    int volume1, volume2, intersect, total;
    double dice, jaccard;
  };

  std::map<int, OverlapMetrics> labelmetrics;

  for (int f = 0 ; f < seg1->nframes ; f++) {
    for (int z = 0 ; z < seg1->depth ; z++) {
      for (int y = 0 ; y < seg1->height ; y++) {
        for (int x = 0 ; x < seg1->width ; x++) {
          int v1 = MRIgetVoxVal(seg1, x, y, z, f);
          int v2 = MRIgetVoxVal(seg2, x, y, z, f);

          labelmetrics[v1].volume1++;
          labelmetrics[v2].volume2++;

          if (v1 == v2) {
            labelmetrics[v1].intersect++;
            labelmetrics[v1].total++;
          }
          else {
            labelmetrics[v1].total++;
            labelmetrics[v2].total++;
          }
        }
      }
    }
  }

  // cleanup
  MRIfree(&seg1);
  MRIfree(&seg2);

  // compute the final measures
  for (auto &lm : labelmetrics) {
    OverlapMetrics &overlap = lm.second;
    overlap.dice = ((double)overlap.intersect * 2) / ((double)overlap.volume1 + (double)overlap.volume2);
    overlap.jaccard = (double)overlap.intersect / (double)overlap.total;
  }

  // retrieve any user-specified labels
  std::vector<int> labels;
  if (parser.exists("labels")) {
    labels = parser.retrieve<std::vector<int>>("labels");
    if (parser.exists("names")) {
      labelnames = parser.retrieve<std::vector<std::string>>("names");
      if (labelnames.size() != labels.size()) {
        fs_fatal(1) << "number of label names must match the number the specified labels";
      }
    }
  }

  // retrieve any user-specified measure options
  bool reportDice = false;
  bool reportJaccard = false;
  if (parser.exists("measures")) {
    // report measures specified by the user
    for(auto const &str : parser.retrieve<std::vector<std::string>>("measures")) {
      if (str == "dice") reportDice = true;
      else if (str == "jaccard") reportJaccard = true;
      else fs_fatal(1) << "unknown measure '" << str << "'. Options are: dice, jaccard, voldiff";
    }
  } else {
    // by default, only report dice scores
    reportDice = true;
  }

  // print results to stdout
  if (!parser.exists("quiet")) {
    std::cout << std::fixed << std::left << std::setprecision(4) << fs::term::bold();
    std::cout << std::setw(7) << "label";
    if (reportDice) std::cout << std::setw(9) << "dice";
    if (reportJaccard) std::cout << std::setw(9) << "jaccard";
    std::cout << fs::term::reset() << std::endl;
    for (auto &lm : labelmetrics) {
      std::cout << std::setw(7) << lm.first;
      if (reportDice) std::cout << std::setw(9) << lm.second.dice;
      if (reportJaccard) std::cout << std::setw(9) << lm.second.jaccard;
      std::cout << std::endl;
    }
  }

  // create detailed json report
  if (parser.exists("out")) {
    // nlohmann::json json;
    // json["inputs"] = {seg1_fname, seg2_fname};
    fs_warning << "json output not implemented yet";
    // std::ofstream out(parser.retrieve<std::string>("out"));
    // out << std::setw(4) << json << std::endl;
  }

  exit(0);
}
