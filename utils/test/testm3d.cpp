//
// testm3d.cpp
//

#include <iostream>
#include <iomanip>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>

extern "C" {
#include "mri.h"
#include "gcamorph.h"
#include "transform.h"
#include "error.h"

  char *Progname="testm3d";
}

using namespace std;

int main(int argc, char *argv[])
{
  cout << "m3d file: ";
  string infile;
  cin >> infile;

  cout << "mri file: ";
  string mrifile;
  cin >> mrifile;

  TRANSFORM *transform=0;
  cout << "reading transform file " << infile.c_str() << endl;
  transform = TransformRead(const_cast<char *> (infile.c_str()));
  if (!transform)
    ErrorExit(ERROR_NOFILE, "%s: could not read transform from file %s",
	      Progname, const_cast<char *> (infile.c_str()));
  // change the transform to vox-to-vox
  if (transform->type != MORPH_3D_TYPE)
    ErrorExit(ERROR_BADPARM, "%s: could not read transform from file %s",
	      Progname, const_cast<char *> (infile.c_str()) );

  cout << "reading mri file " << mrifile.c_str() << endl;
  MRI *mri = MRIread(const_cast<char *> (mrifile.c_str()));

  // modify transform to store inverse also
  cout << "TransformInvert processing ..." << endl;
  TransformInvert(transform, mri);

  cout << "Free memory..." << endl;
  MRIfree(&mri);
  TransformFree(&transform);

  cout << "Done" << endl;
}
