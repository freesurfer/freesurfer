#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <algorithm>
#include <utility>

#include "argparse.h"
#include "pointset.h"

 
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "volcluster.h"


#include "mri_refine_seg.help.xml.h"


typedef std::vector<int> intvec;


// Utility function for finding an int in a vector of ints.
static bool vectorContains(intvec vec, int item)
{
  return std::find(vec.begin(), vec.end(), item) != vec.end();
}


// Simple structure to represent a label component
struct Component {
  int label;  // main label associated with component
  int max;  // max number of clusters allowed (default: 1)
  int additional;  // additional label (like hypointensities) used to define a component
  Component(int l) : label(l), max(1), additional(0) {};
};


// Finds clusters of a label set. Being able to pass a mask buffer prevents having to reallocate
// a new mask volume everytime this function gets called.
static VOLCLUSTER** getLabelClusters(MRI *seg, MRI *mask_buffer, intvec labels, int *nclusters)
{
  // mask away non-matching labels in the segmentation
  MRIbinarizeMatch(seg, &labels[0], labels.size(), 0, mask_buffer);
  return clustGetClusters(mask_buffer, 0, 1, 1, 1, 0, mask_buffer, nclusters, NULL);
}


// Returns the most common neighbor of a voxel. If the most common neighbor is
// of the same value, then the second-most common value (if exists) is returned.
// Diagonal neighbors are not included.
static int commonNeighbor(MRI *seg, int c, int r, int s)
{
  // get frequencies of each neighboring value
  std::map<int,int> frequencies;
  if (c-1 >= 0)          frequencies[MRIgetVoxVal(seg, c-1, r,   s,   0)]++;
  if (r-1 >= 0)          frequencies[MRIgetVoxVal(seg, c,   r-1, s,   0)]++;
  if (s-1 >= 0)          frequencies[MRIgetVoxVal(seg, c,   r,   s-1, 0)]++;
  if (c+1 < seg->width)  frequencies[MRIgetVoxVal(seg, c+1, r,   s,   0)]++;
  if (r+1 < seg->height) frequencies[MRIgetVoxVal(seg, c,   r+1, s,   0)]++;
  if (s+1 < seg->depth)  frequencies[MRIgetVoxVal(seg, c,   r,   s+1, 0)]++;

  // transfer to vector of pairs and sort because maps are annoying
  std::vector<std::pair<int,int>> neighbors;
  for (auto itr = frequencies.begin(); itr != frequencies.end(); ++itr) neighbors.push_back(*itr);
  sort(neighbors.begin(), neighbors.end(), [=](const std::pair<int, int>& a, const std::pair<int, int>& b){return a.second < b.second;});

  // if the most frequent neighbor is the value of the center voxel, return the second-most frequent
  int most_frequent_neighor = neighbors[0].first;
  if ((most_frequent_neighor == MRIgetVoxVal(seg, c, r, s, 0)) && (neighbors.size() > 1)) {
    most_frequent_neighor = neighbors[1].first;
  }
  return most_frequent_neighor;
}


// Recodes cluster voxels (of a particular label) with the values of their most-common neighbors.
static void recodeCluster(MRI *seg, VOLCLUSTER *clr, int label)
{
  // init list of voxel indices to be recoded
  intvec fill_indices(clr->nmembers);
  for (unsigned int i = 0; i < fill_indices.size(); i++) fill_indices[i] = i;

  // continuously loop through voxels in cluster until every valid voxel has been recoded
  while (fill_indices.size() != 0) {
    for (auto it = fill_indices.begin(); it != fill_indices.end(); ) {
      // get coordinates
      int c = clr->col[*it]; int r = clr->row[*it]; int s = clr->slc[*it];
      // ignore voxel if it does not match the label we're recoding
      if (MRIgetVoxVal(seg, c, r, s, 0) != label) {
        it = fill_indices.erase(it);
        continue;
      }
      // find the common neighboring value. If the voxel is surrounded by
      // like-voxels, keep the index in the fill-list and get back to it later
      int majority_vote = commonNeighbor(seg, c, r, s);
      if (majority_vote != label) {
        MRIsetVoxVal(seg, c, r, s, 0, majority_vote);
        it = fill_indices.erase(it);
      } else {
        ++it;
      }
    }
  }
}


// Fills a cluster with a specific value.
static void fillCluster(MRI *seg, VOLCLUSTER *clr, int label)
{
  for (int i = 0 ; i < clr->nmembers ; i++) {
    int c = clr->col[i]; int r = clr->row[i]; int s = clr->slc[i];
    MRIsetVoxVal(seg, c, r, s, 0, label);
  }
}


// Returns true if cluster contains the specified label.
static bool clusterContains(MRI *seg, VOLCLUSTER *clr, int label)
{
  for (int i = 0 ; i < clr->nmembers ; i++) {
    int c = clr->col[i]; int r = clr->row[i]; int s = clr->slc[i];
    if (MRIgetVoxVal(seg, c, r, s, 0) == label) return true;
  }
  return false;
}


// Returns true if cluster is only surrounded by the specified labels.
static bool clusterSurrounded(MRI *seg, VOLCLUSTER *clr, intvec labels)
{
  for (int i = 0 ; i < clr->nmembers ; i++) {
    intvec surrounding;
    int c = clr->col[i]; int r = clr->row[i]; int s = clr->slc[i];
    if (c-1 >= 0)          surrounding.push_back(MRIgetVoxVal(seg, c-1, r,   s,   0));
    if (r-1 >= 0)          surrounding.push_back(MRIgetVoxVal(seg, c,   r-1, s,   0));
    if (s-1 >= 0)          surrounding.push_back(MRIgetVoxVal(seg, c,   r,   s-1, 0));
    if (c+1 < seg->width)  surrounding.push_back(MRIgetVoxVal(seg, c+1, r,   s,   0));
    if (r+1 < seg->height) surrounding.push_back(MRIgetVoxVal(seg, c,   r+1, s,   0));
    if (s+1 < seg->depth)  surrounding.push_back(MRIgetVoxVal(seg, c,   r,   s+1, 0));
    for (int label : surrounding) {
      if (!vectorContains(labels, label)) return false;
    }
  }
  return true;
}


// Returns true if voxel directly touches specific labels.
static bool voxelTouches(MRI *seg, intvec labels, int c, int r, int s)
{
  intvec surrounding;
  if (c-1 >= 0)          surrounding.push_back(MRIgetVoxVal(seg, c-1, r,   s,   0));
  if (r-1 >= 0)          surrounding.push_back(MRIgetVoxVal(seg, c,   r-1, s,   0));
  if (s-1 >= 0)          surrounding.push_back(MRIgetVoxVal(seg, c,   r,   s-1, 0));
  if (c+1 < seg->width)  surrounding.push_back(MRIgetVoxVal(seg, c+1, r,   s,   0));
  if (r+1 < seg->height) surrounding.push_back(MRIgetVoxVal(seg, c,   r+1, s,   0));
  if (s+1 < seg->depth)  surrounding.push_back(MRIgetVoxVal(seg, c,   r,   s+1, 0));
  for (int label : surrounding) {
    if (vectorContains(labels, label)) return true;
  }
  return false;
}


int main(int argc, char **argv) {

  // ------ setup ------

  // parse arguments
  ArgumentParser parser;
  // required
  parser.addArgument("-i", "--in",  1, String, true);
  parser.addArgument("-o", "--out", 1, String, true);
  // optional
  parser.addArgument("--debug");
  // help text
  parser.addHelp(mri_refine_seg_help_xml, mri_refine_seg_help_xml_len);
  parser.parse(argc, argv);

  // load input
  std::string segname = parser.retrieve<std::string>("in");
  MRI *seg = MRIread(segname.c_str());
  if (!seg) fs::fatal() << "could not read input volume " << segname;
  MRI *orig_seg = MRIcopy(seg, NULL);

  // allocate buffer for label mask (used for efficient cluster-finding)
  MRI *mask_buffer = MRIclone(seg, NULL);
  int nclusters;
  VOLCLUSTER **clusters;

  // labels to exclude from stray-cluster refinement
  intvec exclude_labels = {
    5, 44,  // inferior lateral ventricle
    7, 46,  // cerebellum white matter
    30, 62, // vessel
    31, 63, // choroid plexus
    24,     // CSF
    77,     // wm hypointensities
    80,     // non-wm hypointensities
    85,     // optic chiasm
    165,    // skull
    258,    // soft nonbrain tissue
    259     // fluid in eyes
  };

  // structures that might contain (non-wm) hypointensities
  intvec subcortical_gm_labels = {
    10, 49,  // thalamus
    11, 50,  // caudate
    12, 51,  // putamen
    13, 52,  // pallidum
    17, 53,  // hippocampus
    18, 54   // amygdala
  };

  // create list of valid components to refine by looping through the unique labels in the input seg
  int num_unique;
  int *unique_labels = MRIsegIdListNot0(seg, &num_unique, 0);
  std::vector<Component> components;
  for (int i = 0 ; i < num_unique ; i++) {
    int label = unique_labels[i];
    // create a valid component if the label isn't excluded from stray-cluster refinement
    if (vectorContains(exclude_labels, label)) continue;
    Component comp(label);
    // allow for non-wm hypos to be included in the definition of some gm structures
    if (vectorContains(subcortical_gm_labels, label)) comp.additional = 80;
    // allow for wm hypos to be included in the white matter definition
    if (vectorContains({2, 41}, label)) comp.additional = 77;
    // increase the allowed number of clusters for lateral ventricles
    if (vectorContains({4, 43}, label)) comp.max = 2;
    // add component
    components.push_back(comp);
  }

  // ------ main refinement loop ------

  // loop through each component and recode any extra clusters using majority vote

  std::cout << "recoding stray components" << std::endl;

  bool stray_labels = true;
  while (stray_labels) {
  
    stray_labels = false;
    for (Component &comp : components) {

      // find component clusters
      intvec labels = {comp.label};
      if (comp.additional > 0) labels.push_back(comp.additional);
      clusters = getLabelClusters(seg, mask_buffer, labels, &nclusters);

      // recode the extra clusters
      for (int clustid = comp.max ; clustid < nclusters ; clustid++) {
        // for components that are defined by multiple labels (i.e. structures that also
        // include hypointensities), make sure the cluster actually contains some of the
        // primary structure before recoding
        if (clusterContains(seg, clusters[clustid], comp.label)) {
          recodeCluster(seg, clusters[clustid], comp.label);
          stray_labels = true;
        }
      }

      // free clusters
      clustFreeClusterList(&clusters, nclusters);
    }
  }

  // ------ correct hypointensities ------

  // recode any non-wm hypo clusters that are surrounded by wm and ventricle to wm hypos
  clusters = getLabelClusters(seg, mask_buffer, {80}, &nclusters);
  std::cout << "correcting non-wm hypointensities" << std::endl;

  for (int clustid = 0 ; clustid < nclusters ; clustid++) {
    if (clusterSurrounded(seg, clusters[clustid], {80, 2, 41, 4, 43})) fillCluster(seg, clusters[clustid], 77);
  }
  clustFreeClusterList(&clusters, nclusters);

  // ------ correct CSF ------

  // recode any small (< 10 voxels) CSF clusters surrounded by brain structures
  std::cout << "refining stray CSF voxels in brain matter" << std::endl;

  // first create a list of valid brain structures
  intvec brain_structures;
  for (int i = 0 ; i < num_unique ; i++) {
    int label = unique_labels[i];
    if (vectorContains({
      7, 46,  // cerebellum white matter
      8, 47,  // cerebellum cortex
      85,     // optic chiasm
      165,    // skull
      170,    // brainstem
      258,    // soft nonbrain tissue
      259     // fluid in eyes
    }, label)) {
      brain_structures.push_back(label);
    }
  }

  // then loop through the CSF clusters
  clusters = getLabelClusters(seg, mask_buffer, {24}, &nclusters);

  for (int clustid = 0 ; clustid < nclusters ; clustid++) {
    VOLCLUSTER* clr = clusters[clustid];
    if ((clr->nmembers < 10) && clusterSurrounded(seg, clr, brain_structures)) {
      recodeCluster(seg, clr, 24);
    }
  }
  clustFreeClusterList(&clusters, nclusters);

  // ------ correct putamen/insula boundary ------

  // // a small strip of white matter called the extreme capsule seperates the putamen and
  // // insula, and this can be reflected in segmentations by recoding putamen voxels that neighbor insula
  // std::cout << "correcting extreme capsule" << std::endl;
  // clusters = getLabelClusters(seg, mask_buffer, {12, 51}, &nclusters);

  // for (int clustid = 0 ; clustid < nclusters ; clustid++) {
  //   VOLCLUSTER* clr = clusters[clustid];
  //   for (int i = 0 ; i < clr->nmembers ; i++) {
  //     int c = clr->col[i]; int r = clr->row[i]; int s = clr->slc[i];
  //     // if voxel touches cortex, recode to wm
  //     if (voxelTouches(seg, {3, 42}, c, r, s)) {
  //       if (MRIgetVoxVal(seg, c, r, s, 0) == 12) {
  //         // left wm
  //         MRIsetVoxVal(seg, c, r, s, 0, 2);
  //       } else {
  //         // right wm
  //         MRIsetVoxVal(seg, c, r, s, 0, 41);
  //       }
  //     }      
  //   }
  // }
  // clustFreeClusterList(&clusters, nclusters);

  // ------ review what has been changed ------

  std::cout << "reviewing refinement" << std::endl;

  // check how much was changed from the original
  int did_change = false;
  for (int c = 0; c < seg->width; c++) {
    for (int r = 0; r < seg->height; r++) {
      for (int s = 0; s < seg->depth; s++) {
        int origval = MRIgetVoxVal(orig_seg, c, r, s, 0);
        int newval  = MRIgetVoxVal(seg, c, r, s, 0);
        if (origval != newval) {
          did_change = true;
          MRIsetVoxVal(mask_buffer, c, r, s, 0, 1);
        } else {
          MRIsetVoxVal(mask_buffer, c, r, s, 0, 0);
        }
      }
    }
  }

  if (!did_change) {
    std::cout << "note: no voxels were updated during refinement!" << std::endl;
  }

  //  ------ save outputs -----

  // save refined seg
  std::string outname = parser.retrieve<std::string>("out");
  std::cout << "saving refined seg as " << outname << std::endl;
  MRIwrite(seg, outname.c_str());

  // if debugging is turned on, save the original segmentation, the difference mask, and a
  // json pointset that points to the changed regions
  if (parser.exists("debug")) {

    std::cout << "saving debug pointset, mask, and original segmentation" << std::endl;

    // saving a pointset representing each changed voxel will be quite large (and tedious to
    // use), so instead we'll compute clusters of the change-mask and save a point from there
    clusters = clustGetClusters(mask_buffer, 0, 1, 1, 1, 0, mask_buffer, &nclusters, NULL);
    fsPointSet changed_voxels;
    for (int clustid = 0 ; clustid < nclusters ; clustid++) {
      VOLCLUSTER* clr = clusters[clustid];
      double x, y, z;
      MRIvoxelToWorld(seg, clr->col[0], clr->row[0], clr->slc[0], &x, &y, &z);
      changed_voxels.add(fsPointSet::Point(x, y, z));
    }
    clustFreeClusterList(&clusters, nclusters);

    // get output directory
    std::string dirname = fio_dirname(outname.c_str());

    // save debugging outputs
    changed_voxels.save(dirname + "/refinement.points.json");
    MRIwrite(orig_seg,    std::string(dirname + "/refinement.orig.mgz").c_str());
    MRIwrite(mask_buffer, std::string(dirname + "/refinement.mask.mgz").c_str());
  }

  // ------ cleanup ------

  // free
  MRIfree(&seg);
  MRIfree(&orig_seg);
  MRIfree(&mask_buffer);

  std::cout << "mri_refine_seg done" << std::endl;

  return 0;
}
