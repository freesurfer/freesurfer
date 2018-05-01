// 
// unit test for MRIScomputeBorderValues - located in utils/mrisurf.c
//

#include <string>
#include <iostream>
#include <iomanip>

#ifdef __cplusplus
extern "C"
{
  #endif

  #include "error.h"
  #include "utils.h"
  #include "macros.h"
  #include "mri.h"
  #include "mrisurf.h"

  #ifdef __cplusplus
}
#endif

int main(int argc, char *argv[])
{
  int err = 0;

  // check arg count:
  if (argc != 5)
  {
    std::cerr << "ERROR: 4 args are required: mris mri_brain mri_smooth mri_aseg\n";
    exit(1);
  }


  // read args:
  std::string Progname = argv[0];
  std::string s_mris = argv[1];
  std::string s_brain = argv[2];
  std::string s_smooth = argv[3];
  std::string s_aseg = argv[4];

  std::cout << Progname << std::endl;


  // read MRIS input:
  std::cout << "reading " << s_mris << std::endl;
  MRIS *mris = MRISread(s_mris.c_str());
  if (!mris)
  {
    std::cerr << "ERROR: could not read mris '" << s_mris << "'!\n";
    exit(ERROR_BADPARM);
  }
  // read MRI_BRAIN input:
  std::cout << "reading " << s_brain << std::endl;
  MRI *mri_brain = MRIread(s_brain.c_str());
  if (!mri_brain)
  {
    std::cerr << "ERROR: could not read mri_brain '" << s_brain << "'!\n";
    exit(ERROR_BADPARM);
  }
  // read MRI_SMOOTH input:
  std::cout << "reading " << s_smooth << std::endl;
  MRI *mri_smooth = MRIread(s_smooth.c_str());
  if (!mri_smooth)
  {
    std::cerr << "ERROR: could not read mri_smooth '" << s_smooth << "'!\n";
    exit(ERROR_BADPARM);
  }
  // read MRI_ASEG input:
  std::cout << "reading " << s_aseg << std::endl;
  MRI *mri_aseg = MRIread(s_aseg.c_str());
  if (!mri_aseg)
  {
    std::cerr << "ERROR: could not read mri_aseg '" << s_aseg << "'!\n";
    exit(ERROR_BADPARM);
  }


  // run:
  std::cout << "running MRIScomputeBorderValues...\n";
  err = MRIScomputeBorderValues(mris, mri_brain, mri_smooth,
                              120.0, 112.28, 64.0, 53.24, 112.28,
                              1.0, 10.0, NULL, 1, NULL, 0.0, 0, mri_aseg);
  if (err)
  {
    std::cerr << "ERROR: could not run MRIScomputeBorderValues!\n";
    exit(err);
  }


  // check averages of val, d, and mean across all vertices with nonzero val:
  int vno, vtot = 0;
  VERTEX *v;
  float val_sum = 0, d_sum = 0, mean_sum = 0;
  float val_avg, d_avg, mean_avg;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno];
    if (!FZERO(v->val))
    {
      val_sum  += v->val;
      d_sum    += v->d;
      mean_sum += v->mean;

      vtot++;
    }
  }

  val_avg = val_sum / vtot;
  d_avg =  d_sum / vtot;
  mean_avg = mean_sum / vtot;

  std::cout << "\ncomputing stats for comparison:\n";
  std::cout << std::fixed;
  std::cout << std::setprecision(8);
  std::cout << "avg val:  " << val_avg << std::endl;
  std::cout << "avg d:    " << d_avg << std::endl;
  std::cout << "avg mean: " << mean_avg << std::endl;

  if (!FEQUAL(val_avg, 80.09034729)) { err = 1; }
  if (!FEQUAL(d_avg, -0.26023445))   { err = 1; }
  if (!FEQUAL(mean_avg, 6.81196499)) { err = 1; }
  if (err == 1)
  {
    std::cout << "surface vertex stats DO NOT match reference data!\n";
  }


  // shut down:
  MRISfree(&mris);
  MRIfree(&mri_brain);
  MRIfree(&mri_smooth);
  MRIfree(&mri_aseg);
  
  exit(err);
}