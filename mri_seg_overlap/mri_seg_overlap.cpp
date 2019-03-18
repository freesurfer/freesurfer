#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <map>

#include "mri_seg_overlap.help.xml.h"
#include "argparse.h"
#include "lut.h"
#include "json.h"
#include "log.h"

 
#include "mri.h"
#include "mri2.h"



struct IntermediateMetrics {
  double volume1 = 0, volume2 = 0, intersect = 0, total = 0;
};


struct OverlapMeasure {
  typedef double (*computefunc)(IntermediateMetrics&);
  OverlapMeasure(std::string _name, computefunc _func) : name(_name), compute(_func) {};
  std::string name;
  computefunc compute;
  std::map<int, double> labels;  // measures for each label
  double mean = 0, std = 0;
  double wmean = 0;  // volume-weighted mean
  double wmean_sc = 0;  // exclude wm and cortex
};


static double computeDice(IntermediateMetrics &im) {
  return (im.intersect * 2) / (im.volume1 + im.volume2);
}

static double computeJaccard(IntermediateMetrics &im) {
  return im.intersect / im.total;
}


int main(int argc, const char **argv) 
{
  // ------ parse arguments ------

  ArgumentParser parser;
  parser.addHelp(mri_seg_overlap_help_xml, mri_seg_overlap_help_xml_len);
  parser.addArgument("seg1");
  parser.addArgument("seg2");
  parser.addArgument("-o", "--out", 1);
  parser.addArgument("-m", "--measures", '+');
  parser.addArgument("-l", "--labels", '+', Int);
  parser.addArgument("-f", "--labelfile", 1);
  parser.addArgument("-n", "--names", '+');
  parser.addArgument("-x", "--no-names");
  parser.addArgument("-s", "--seg");
  parser.addArgument("-q", "--quiet");
  parser.parse(argc, argv);

  // ------ load inputs ------

  std::string seg1_fname = parser.retrieve<std::string>("seg1");
  MRI *seg1 = MRIread(seg1_fname.c_str());
  if (!seg1) logFatal(1) << "could not read input volume " << seg1;

  std::string seg2_fname = parser.retrieve<std::string>("seg2");
  MRI *seg2 = MRIread(seg2_fname.c_str());
  if (!seg2) logFatal(1) << "could not read input volume " << seg2;

  // check input dimensions
  if (MRIdimMismatch(seg1, seg2, 0))  logFatal(1) << "input volumes must have matching dimensions";
  if (seg1->nframes != seg2->nframes) logFatal(1) << "input volumes must have the same number of frames";

  // ------ retrieve user options ------

  // todoc
  std::vector<OverlapMeasure> measures;
  
  // measure options
  if (parser.exists("measures")) {
    // report measures specified by the user
    for (auto const &str : parser.retrieve<std::vector<std::string>>("measures")) {
      if (str == "dice") measures.push_back(OverlapMeasure("dice", &computeDice));
      else if (str == "jaccard") measures.push_back(OverlapMeasure("jaccard", &computeJaccard));
      else logFatal(1) << "unknown measure '" << str << "'... options are: dice, jaccard";
    }
  } else {
    // by default, only report dice scores
    measures.push_back(OverlapMeasure("dice", &computeDice));
  }

  // a few sanity checks on the user input
  if (parser.exists("labelfile") && parser.exists("labels")) {
    logFatal(1) << "can't use both the --labels and --labelfile options together";
  } else if (parser.exists("names") && !parser.exists("labels")) {
    logFatal(1) << "the --names option must be used along with the --labels option";
  } else if (parser.exists("seg") && (parser.exists("labelfile") || parser.exists("labels"))) {
    logFatal(1) << "can't specify --seg in addition to a custom label list";
  }

  // check if user wants to ignore label names
  bool reportNames = !parser.exists("no-names");

  // determine which labels to report on
  bool all_labels = false;
  std::vector<int> labels;
  std::map<int, std::string> labelnames;
  if (parser.exists("labels")) {
    // user has specified labels on the command line
    labels = parser.retrieve<std::vector<int>>("labels");
    if (reportNames && parser.exists("names")) {
      // user provided custom label names as well 
      std::vector<std::string> names = parser.retrieve<std::vector<std::string>>("names");
      if (names.size() != labels.size()) {
        logFatal(1) << "number of label names (" << labelnames.size() << ") must match "
                       "the number of specified labels (" << labels.size() << ")";
      }
      for (unsigned int i = 0 ; i < labels.size() ; i++) labelnames[labels[i]] = names[i];
    }
  } else if (parser.exists("labelfile")) {
    // user specified labels via a file - assume lookup-table format
    LookupTable lut(parser.retrieve<std::string>("labelfile"));
    if (lut.empty()) logFatal(1) << "provided label file contains no valid labels";
    labels = lut.labels();
    if (lut.hasNameInfo()) for (int i : labels) labelnames[i] = lut[i].name;
  } else if (parser.exists("seg")) {
    // use the major anatomical segmentation structures
    labels = {
      2,  41,  // Cerebral White Matter
      3,  42,  // Cerebral Cortex
      17, 53,  // Hippocampus
      11, 50,  // Caudate
      12, 51,  // Putamen
      13, 52,  // Pallidum
      18, 54,  // Amygdala
      10, 49,  // Thalamus Proper
      4,  43,  // Lateral Ventricle
      14, 15,  // Third and Fourth Ventricles
      5,  44,  // Inf Lateral Vent
      26, 58   // Accumbens Area
    };
  } else {
    // use all labels found in segs (we'll determine these after already looping through each voxel)
    all_labels = true;
  }

  //------ compute intermediate label metrics ------

  std::map<int, IntermediateMetrics> intermediate_metrics;

  for (int z = 0 ; z < seg1->depth ; z++) {
    for (int y = 0 ; y < seg1->height ; y++) {
      for (int x = 0 ; x < seg1->width ; x++) {
        int v1 = MRIgetVoxVal(seg1, x, y, z, 0);
        int v2 = MRIgetVoxVal(seg2, x, y, z, 0);

        intermediate_metrics[v1].volume1 += 1;
        intermediate_metrics[v2].volume2 += 1;

        if (v1 == v2) {
          intermediate_metrics[v1].intersect += 1;
          intermediate_metrics[v1].total += 1;
        }
        else {
          intermediate_metrics[v1].total += 1;
          intermediate_metrics[v2].total += 1;
        }
      }
    }
  }

  // free the volumes
  MRIfree(&seg1);
  MRIfree(&seg2);

  // ------ finalize labels and label names ------

  if (all_labels) {
    // if we're using all of the unique labels in the volumes, populate the labels
    // vector with the keys of the intermediate_metrics map
    for (auto const &m : intermediate_metrics) if (m.first != 0) labels.push_back(m.first);
  } else {
    // otherwise, remove any labels that were not found in the segmentations
    labels.erase(std::remove_if(labels.begin(), labels.end(), 
      [&intermediate_metrics](int l) { return intermediate_metrics.count(l) < 1; }), labels.end());
  }

  // sanity check to make sure the label list isn't empty
  if (labels.empty()) logFatal(1) << "no matching labels to report on";

  // determine default label names (if not already known) via FreeSurferColorLUT
  if (reportNames && labelnames.empty()) {
    LookupTable lut(std::string(std::getenv("FREESURFER_HOME")) + "/FreeSurferColorLUT.txt");
    if (lut.empty()) {
      logWarning << "can't load default FreeSurferColorLUT - is FREESURFER_HOME set?";
      reportNames = false;
    } else if (lut.hasNameInfo()) {
      for (int l : labels) labelnames[l] = lut[l].name;
    } else {
      reportNames = false;
    }
  }

  // ------ compute measures ------

  double total_volume = 0, subcortical_volume = 0;
  for (int l : labels) {
    IntermediateMetrics &im = intermediate_metrics[l];
    // combine label volumes from both segs (used for weighted mean)
    double combined_volume = im.volume1 + im.volume2;
    // check if the label is wm or cortex so we can created a weighted mean of all other structures
    bool subcortical = false;
    if ((parser.exists("seg")) && (l != 2) && (l != 41) && (l != 3) && (l != 42)) subcortical = true;
    // compute measures
    for (auto &measure : measures) {
      double value = measure.compute(im);
      measure.labels[l] = value;
      measure.mean  += value;
      measure.wmean += value * combined_volume;
      if (subcortical) measure.wmean_sc += value * combined_volume;
    }
    total_volume += combined_volume;
    if (subcortical) subcortical_volume += combined_volume;
  }

  // compute final statistics
  for (auto &measure : measures) {
    // mean and std
    measure.mean /= labels.size();
    double sqr_diff_sum = 0;
    for (auto const &v : measure.labels) sqr_diff_sum += pow(v.second - measure.mean, 2);
    measure.std = sqrt(sqr_diff_sum / labels.size());
    // weighted means
    measure.wmean /= total_volume;
    if (parser.exists("seg")) measure.wmean_sc /= subcortical_volume;
  }

  // ------ report results ------

  // print results to stdout
  if (!parser.exists("quiet")) {
    // determine the longest label name for table formatting
    unsigned int name_width = 7;
    if (reportNames) {
      for (int l : labels) if (labelnames[l].length() + 2 > name_width) name_width = labelnames[l].length() + 2;
    }
    // print the table header
    std::cout << std::fixed << std::left << std::setprecision(4) << term::bold();
    std::cout << std::setw(7) << "label";
    if (reportNames) std::cout << std::setw(name_width) << "names";
    for (auto const &measure : measures) std::cout << std::setw(8) << measure.name;
    std::cout << term::reset() << std::endl;
    // print rows of results
    for (int l : labels) {
      std::cout << std::setw(7) << l;
      if (reportNames) std::cout << std::setw(name_width) << labelnames[l];
      for (auto &measure : measures) std::cout << std::setw(8) << measure.labels[l];
      std::cout << std::endl;
    }
  }

  // create detailed json report
  if (parser.exists("out")) {
    nlohmann::json json;
    // measures
    for (auto &measure : measures) {
      nlohmann::json subjson;
      subjson["mean"] = measure.mean;
      subjson["std"] = measure.std;
      subjson["weighted-mean"] = measure.wmean;
      if (parser.exists("seg")) subjson["weighted-subcortical-mean"] = measure.wmean_sc;
      for (int l : labels) subjson["labels"][std::to_string(l)] = measure.labels[l];
      json["measures"][measure.name] = subjson;
    }
    // label names
    if (reportNames) for (int l : labels) json["names"][std::to_string(l)] = labelnames[l];
    // inputs
    json["inputs"] = {seg1_fname, seg2_fname};
    // write to file
    std::ofstream out(parser.retrieve<std::string>("out"));
    out << std::setw(4) << json << std::endl;
    out.close();
  }

  exit(0);
}
