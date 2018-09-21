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

  // a few sanity checks on the user input
  if (parser.exists("labelfile") && parser.exists("labels")) {
    fs_fatal(1) << "can't use both the --labels and --labelfile options together";
  } else if (parser.exists("names") && !parser.exists("labels")) {
    fs_fatal(1) << "todo";
  } else if (parser.exists("seg") && (parser.exists("labelfile") || parser.exists("labels"))) {
    fs_fatal(1) << "todo";
  }

  // check if user wants to ignore label names
  bool reportNames = !parser.exists("no-names");

  // determine which labels to report on
  std::vector<int> labels;
  std::vector<std::string> labelnames;
  if (parser.exists("labels")) {
    // user specified labels on the command line
    labels = parser.retrieve<std::vector<int>>("labels");
    if (reportNames && parser.exists("names")) {
      // they provided custom label names as well 
      labelnames = parser.retrieve<std::vector<std::string>>("names");
      if (labelnames.size() != names.size()) {
        fs_fatal(1) << "number of label names must match the number of specified labels";
      }
    }
  } else if (parser.exists("labelfile")) {
    // user specified labels via a file - assume lookup-table format
    // todo: io read color lut
    // foreach label in colorlut;
    //   add to labels
  } else if (parser.exists("seg")) {
    // use the major anatomical segmentation structures
    labels = {
      17, 53,  // Hippocampus
      11, 50,  // Caudate
      12, 51,  // Putamen
      13, 52,  // Pallidum
      18, 54,  // Amygdala
      10, 49,  // Thalamus Proper
      4,  43,  // Lateral Ventricle
      14, 15,  // Third and Fourth Ventricles
      5,  44,  // Inf Lateral Vent
      2,  41,  // Cerebral White Matter
      3,  42,  // Cerebral Cortex
      26, 58   // Accumbens Area
    }
  } else {
    // use all valid labels
  }

  // determine the appropriate label names if not already known
  std::vector<std::string> labelnames;
  if (reportNames && labelnames.empty()) {
    // todo: read in the primary FreeSurferColorLUT.txt and get label names
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

  // free the volumes
  MRIfree(&seg1);
  MRIfree(&seg2);

  // compute the final measures
  for (auto &lm : labelmetrics) {
    OverlapMetrics &overlap = lm.second;
    overlap.dice = ((double)overlap.intersect * 2) / ((double)overlap.volume1 + (double)overlap.volume2);
    overlap.jaccard = (double)overlap.intersect / (double)overlap.total;
  }

  // print results to stdout
  if (!parser.exists("quiet")) {
    std::cout << std::fixed << std::left << std::setprecision(4) << fs::term::bold();
    std::cout << std::setw(7) << "label";
    if (reportNames) std::cout << std::setw(9) << "names";
    if (reportDice) std::cout << std::setw(9) << "dice";
    if (reportJaccard) std::cout << std::setw(9) << "jaccard";
    std::cout << fs::term::reset() << std::endl;
    for (auto &lm : labelmetrics) {
      std::cout << std::setw(7) << lm.first;
      if (reportNames) std::cout << std::setw(9) << "names";
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
