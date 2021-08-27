/**
 * @brief two fsl aligned label volume and src volume
 *
 * to get the voxel values of the src corresponding to the label volume
 */
/*
 * Original Author: Doug Greve
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <iostream>
#include <iomanip>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>
#include <fstream>
#include <vector>

 
#include "mri.h"
#include "macros.h"
  const char *Progname="fsl_label2voxel";


using namespace std;

void usage() {
  cout << "Usage: fsl_label2voxel <labelval> <labeled volume> <src volume> <output filename>" << endl;
}

class Row {
public:
  Row(int id, int x, int y, int z, double val):id_(id), x_(x), y_(y), z_(z), val_(val) {}

  void print(ostream &o) {
    o << id_ << "\t" << x_ << "\t" << y_ << "\t" << z_ << "\t" << val_ << endl;
  }

private:
  int id_;
  int x_;
  int y_;
  int z_;
  double val_;
};

int main(int argc, char *argv[]) {
  if (argc < 4) {
    usage();
    return -1;
  }

  double label_val = atof(argv[1]);
  string labelVol = argv[2];
  string srcVol = argv[3];
  string outFile = argv[4];

  MRI *labelMRI = MRIread(const_cast<char *> (labelVol.c_str()));
  if (!labelMRI) {
    cerr << "Could not open " << labelVol.c_str() << endl;
    return -1;
  }
  MRI *srcMRI = MRIread(const_cast<char *> (srcVol.c_str()));
  if (!srcMRI) {
    cerr << "Could not open " << srcVol.c_str() << endl;
    return -1;
  }
  // check the size
  if ((labelMRI->width != srcMRI->width)
      || (labelMRI->height != srcMRI->height)
      || (labelMRI->depth != srcMRI->depth)) {
    cerr << "src and label volume sizes differ" << endl;
    return -1;
  }

  vector<Row> list;
  int counter = 0;
  for (int z = 0; z < srcMRI->depth; z++)
    for (int y=0; y < srcMRI->height; y++)
      for (int x = 0; x < srcMRI->width; x++) {
        float aval = MRIgetVoxVal(labelMRI, x, y, z, 0);
        if (FEQUAL(label_val, aval)) {
          double srcVal = MRIgetVoxVal(srcMRI, x, y, z, 0);
          list.push_back(Row(counter, x, y, z, srcVal));
          counter++;
        }
      }

  // now printout the file
  ofstream ofile(const_cast<char *>(const_cast<char *>(outFile.c_str())));
  ofile << "#!ascii label " << argv[2] << ", from subject " << argv[3] << endl;
  ofile << counter << endl;
  vector<Row>::iterator iter;
  for (iter=list.begin(); iter != list.end(); iter++)
    iter->print(ofile);

  return 0;
}
