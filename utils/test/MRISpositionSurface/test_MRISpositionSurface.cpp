// 
// unit test for MRISpositionSurface - located in utils/mrisurf.c
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
  if (argc != 4)
  { 
    std::cerr << "ERROR: 3 args are required: mris mri_brain mri_smooth\n";
    exit(1);
  }


  // read args:
  std::string Progname = argv[0];
  std::string s_mris = argv[1];
  std::string s_brain = argv[2];
  std::string s_smooth = argv[3];

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

  // set up integration parms:
  INTEGRATION_PARMS parms;
  memset(&parms, 0, sizeof(parms));

  parms.fill_interior = 0;
  parms.projection = NO_PROJECTION;
  parms.tol = 1e-4;
  parms.dt = 0.5f;
  parms.base_dt = parms.dt;
  parms.l_curv = 1.0;
  parms.l_intensity = 0.2;
  parms.l_spring = 1.0f;
  parms.l_tspring = 1.0f;
  parms.l_nspring = 0.5f;
  parms.niterations = 20;
  parms.start_t = 0;
  parms.write_iterations = 0;
  parms.integration_type = INTEGRATE_MOMENTUM;
  parms.momentum = 0.0;
  parms.dt_increase = 1.0;
  parms.dt_decrease = 0.50;
  parms.error_ratio = 50.0;
  parms.l_surf_repulse = 0.0;
  parms.l_repulse = 5;
  parms.sigma = 0.2f;
  parms.n_averages = 10;



  // run:
  std::cout << "running MRISpositionSurface...\n";
  err = MRISpositionSurface(mris, mri_brain, mri_smooth, &parms);
  if (err)
  {
    std::cerr << "ERROR: could not run MRISpositionSurface!\n";
    exit(err);
  }


  // compare data:
  VERTEX *v;

  v = &mris->vertices[160];
  if (!FEQUAL(v->dist[0], 1.36991000f)) { err = 1; }
  v = &mris->vertices[1170];
  if (!FEQUAL(v->dist[0], 1.10133839f)) { err = 1; }
  v = &mris->vertices[720];
  if (!FEQUAL(v->dist[0], 1.51904464f)) { err = 1; }
  if (err == 1)
  {
    std::cout << "surface vertex stats DO NOT match reference data!\n";
    std::cout << mris->vertices[160].dist[0] << " "
              << mris->vertices[1170].dist[0] << " "
              << mris->vertices[720].dist[0] << std::endl;
  }


  // shut down:
  MRISwrite(mris, "positionSurfaceNEW");
  MRISfree(&mris);
  MRIfree(&mri_brain);
  MRIfree(&mri_smooth);

  exit(err);
}
