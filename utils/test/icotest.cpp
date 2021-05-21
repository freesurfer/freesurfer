/**
 * @brief tests MRISwriteICO.
 *
 * read icosahedron tri file and write test.tri file.
 * read back test.tri and compare values.
 */
/*
 * Original Author: Y. Tosa
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
#include <fstream>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <sstream>
#include <string>

extern "C"
{
#include <mrisurf.h>
const char *Progname = "icotest";
}

using namespace std;

int main()
{
  MRIS *mris;
  const char *mri_dir = "../../distribution";
  for (int index=0; index < 8; ++index)
  {
    ostringstream infile, outfile;
    infile << mri_dir << "/lib/bem/ic" << index << ".tri";
    cout << "Reading " << infile.str() << endl;
    mris = MRISread(const_cast<char *> (infile.str().c_str()));
    outfile << "./test" << index << ".tri";
    cout << "Writing " << outfile.str() << endl;
    MRISwrite(mris, const_cast<char *> (outfile.str().c_str()));
    MRISfree(&mris);
  }
  // compare
  for (int index=0; index < 8; ++index)
  {
    ostringstream infile, outfile;
    infile << mri_dir << "/lib/bem/ic" << index << ".tri";
    outfile << "./test" << index << ".tri";
    cout << "Comparing " << outfile.str() << " with " << infile.str() << endl;

    ifstream origf(infile.str().c_str(), ios::in);
    ifstream outf(outfile.str().c_str(), ios::in);
    char buf[512];
    char buf1[512];
    // first check number of vertices
    origf.getline(buf, sizeof(buf));
    stringstream st(buf);
    int nvertices;
    st >> nvertices;
    outf.getline(buf1, sizeof(buf1));
    stringstream st1(buf1);
    int nvertices1;
    st1 >> nvertices1;
    if (nvertices != nvertices1)
    {
      cerr << "Number of vertices differ: orig " << nvertices << "  new " << nvertices1 << endl;
      return -1;
    }
    // read each vertices
    int idx;
    float x, y, z;
    int idx1;
    float x1, y1, z1;
    for (int i=0; i < nvertices; ++i)
    {
      origf.getline(buf, sizeof(buf));
      stringstream st(buf);
      st >> idx >> x >> y >> z;
      outf.getline(buf1, sizeof(buf1));
      stringstream st1(buf1);
      st1 >> idx1 >> x1 >> y1 >> z1;
      if (idx != idx1 || x != x1 || y != y1 || z != z1)
      {
        cerr << "vertices values do not agree " << endl;
        cerr << "orig: idx " << idx << "  x " << x << " y " << y << " z " << z << endl;
        cerr << "new : idx " << idx1 << "  x " << x1 << " y " << y1 << " z " << z1 << endl;
        return -1;
      }
    }
    // next faces with 3 vertices
    origf.getline(buf, sizeof(buf));
    outf.getline(buf1, sizeof(buf1));
    int faces, faces1;
    stringstream st2(buf);
    st2 >> faces;
    stringstream st3(buf1);
    st3 >> faces1;
    if (faces != faces1)
    {
      cerr << "orig :" << buf << "   new: " << buf1 << endl;
      cerr << "number of faces do not agree: orig " << faces << "   new " << faces1 << endl;
      return -1;
    }
    int v1, v2, v3;
    int v11, v12, v13;
    for (int j=0; j < faces ; ++j)
    {
      origf.getline(buf, sizeof(buf));
      stringstream st4(buf);
      st4 >> idx >> v1 >> v2 >> v3;
      outf.getline(buf1, sizeof(buf1));
      stringstream st5(buf1);
      st5 >> idx1 >> v11 >> v12 >> v13;
      if (idx != idx1 || v1 != v11 || v2 != v12 || v3 != v13)
      {
        cerr << "vertices values do not agree " << endl;
        cerr << "orig: idx " << idx << "  v1 " << v1 << " v2 " << v2 << " v3 " << v3 << endl;
        cerr << "new : idx " << idx1 << "  v11 " << v11 << " v12 " << v12 << " v13 " << v13 << endl;
        return -1;
      }
    }
    // remove outputfile
    stringstream command;
    command << "rm " << outfile.str();
    system(command.str().c_str());
  }
  cout << "No difference found" << endl;
}
