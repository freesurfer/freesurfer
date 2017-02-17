// 
// unit test for mriBuildVoronoiDiagramFloat - located in utils/mrinorm.c
//

#include <string>
#include <iostream>

#ifdef __cplusplus
extern "C"
{
  #endif

  #include "error.h"
  #include "utils.h"
  #include "macros.h"
  #include "mri.h"
  #include "mrinorm.h"

  #ifdef __cplusplus
}
#endif


int main(int argc, char *argv[])
{
  // check arg count:
  if (argc != 4)
  {
    std::cerr << "ERROR: 3 arguments are required: mri_src mri_ctrl mri_dst\n";
    exit(1);
  }

  // read args:
  std::string Progname = argv[0];
  std::string s_src = argv[1];
  std::string s_ctrl = argv[2];
  std::string s_dst = argv[3];

  std::cout << Progname << std::endl;

  // read MRI input:
  std::cout << "reading " << s_src << std::endl;
  MRI *mri_src = MRIread(s_src.c_str());
  if (!mri_src)
  {
    std::cerr << "ERROR: could not read src '" << s_src << "'!\n";
    exit(ERROR_BADPARM);
  }
  std::cout << "reading " << s_ctrl << std::endl;
  MRI *mri_ctrl = MRIread(s_ctrl.c_str());
  if (!mri_ctrl)
  {
    std::cerr << "ERROR: could not read ctrl '" << s_ctrl << "'!\n";
    exit(ERROR_BADPARM);
  }

  // make sure to test mriBuildVoronoiDiagramFloat:
  if (mri_src->type != MRI_FLOAT)
  {
  	std::cerr << "ERROR: mri_src type is not float!\n";
  	exit(1);
  }

  // run:
  std::cout << "running MRIbuildVoronoiDiagram...\n";
  MRI *mri_dst = MRIbuildVoronoiDiagram(mri_src, mri_ctrl, NULL);
  if (!mri_dst)
  {
    std::cerr << "ERROR: could not run mriBuildVoronoiDiagramFloat!\n";
    exit(1);
  }

  // save:
  std::cout << "writing " << s_dst << std::endl;
  MRIwrite(mri_dst, s_dst.c_str());

  // shut down:
  MRIfree(&mri_src);
  MRIfree(&mri_ctrl);
  MRIfree(&mri_dst);
  
  exit(0);
}