//
// mris_info.cpp
//

#include <iostream>
#include <iomanip>
#include <vector>
extern "C" {
#include "mri.h"
#include "transform.h"
#include "mrisurf.h"

  char *Progname = "mris_info";
}

using namespace std;






int main(int argc, char *argv[])
{
  vector<string> type;
  type.push_back("MRIS_BINARY_QUADRANGLE_FILE");
  type.push_back("MRIS_ASCII_TRIANGLE_FILE");
  type.push_back("MRIS_GEO_TRIANGLE_FILE");
  type.push_back("MRIS_ICO_SURFACE");
  type.push_back("MRIS_ICO_FILE");
  type.push_back("MRIS_VTK_FILE");
  
  if (argc < 1)
  {
    cout << "Usage: mris_info <surface>" << endl;
    return -1;
  }
  MRIS *mris = MRISread(argv[1]);
  if (!mris)
  {
    cerr << "could not open " << argv[1] << endl;
    return -1;
  }
  cout << "SURFACE INFO ================================================== " << endl;
  cout << "type        : " << type[mris->type].c_str() << endl;
  cout << "num vertices: " << mris->nvertices << endl;
  cout << "num faces   : " << mris->nfaces << endl;
  cout << "num stripgs : " << mris->nstrips << endl;
  cout << "ctr         : (" << mris->xctr << ", " << mris->yctr << ", " << mris->zctr << ")" << endl;
  if (mris->lta)
  {
    cout << "talairch.xfm: " << endl;
    MatrixPrint(stdout, mris->lta->xforms[0].m_L);
    cout << "surfaceRAS to talaraiched surfaceRAS: " << endl;
    MatrixPrint(stdout, mris->SRASToTalSRAS_);
    cout << "talairached surfaceRAS to surfaceRAS: " << endl;
    MatrixPrint(stdout, mris->TalSRASToSRAS_);
  }
  vg_print(&mris->vg); 

  MRISfree(&mris);
}
