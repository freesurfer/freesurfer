//
// mri_modify.cpp
//
// author: Yasunari Tosa (tosa@nmr.mgh.harvard.edu)
// purpose: modify direction cosine info on the volume
// 
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: tosa $
// Revision Date  : $Date: 2004/08/27 15:23:30 $
// Revision       : $Revision: 1.3 $

#include <iostream>
#include <iomanip>

extern "C" {
#include "macros.h"
#include "mri.h"
#include "transform.h"
#include "version.h"
  char *Progname = "mri_modify";
}

using namespace std;

double gflip_angle=0;
float gte=0;
float gtr=0;
float gti=0;

void print_usage()
{
  cout << "Usage: mri_modify <-xras xr xa xs> <-yras yr ya ys> <-zras zr za zs> <-cras cr ca cs> \\ " << endl;
  cout << "                  <-xsize size> <-ysize size> <-zsize size> \\ " << endl;
  cout << "                  <-tr recoverytime> <-te echotime> <-ti inversiontime> <-fa angledegree> \\ " << endl;
  cout << "                  involume outvolume" << endl;
}

int get_option(int argc, char *argv[], VOL_GEOM &vg)
{
  int  nargs = 0 ;
  char *option ;
  option = argv[1] + 1 ;            /* past '-' */
  if (!strcmp(option, "-help"))
  {
    print_usage();
    exit(0);
  }
  else if (!strcmp(option, "xras"))
  {
    vg.x_r = atof(argv[2]);
    vg.x_a = atof(argv[3]);
    vg.x_s = atof(argv[4]);
    nargs=3;
  }
  else if (!strcmp(option, "yras"))
  {
    vg.y_r = atof(argv[2]);
    vg.y_a = atof(argv[3]);
    vg.y_s = atof(argv[4]);
    nargs=3;
  }
  else if (!strcmp(option, "zras"))
  {
    vg.z_r = atof(argv[2]);
    vg.z_a = atof(argv[3]);
    vg.z_s = atof(argv[4]);
    nargs=3;
  }
  else if (!strcmp(option, "cras"))
  {
    vg.c_r = atof(argv[2]);
    vg.c_a = atof(argv[3]);
    vg.c_s = atof(argv[4]);
    nargs=3;
  }
  else if (!strcmp(option, "xsize"))
  {
    vg.xsize = atof(argv[2]);
    nargs=1;
  }
  else if (!strcmp(option, "ysize"))
  {
    vg.ysize = atof(argv[2]);
    nargs=1;
  }
  else if (!strcmp(option, "zsize"))
  {
    vg.zsize = atof(argv[2]);
    nargs=1;
  }
  else if (!strcmp(option, "tr"))
  {
    gtr=atof(argv[2]);
    nargs=1;
  }
  else if (!strcmp(option, "te"))
  {
    gte=atof(argv[2]);
    nargs=1;
  }
  else if (!strcmp(option, "ti"))
  {
    gti=atof(argv[2]);
    nargs=1;
  }
  else if (!strcmp(option, "fa"))
  {
    // mri stores it as radian
    gflip_angle=RADIANS(atof(argv[2]));
    nargs=1;
  }
  return nargs;
}

int main(int argc, char *argv[])
{
  int nargs;
  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_modify.cpp,v 1.3 2004/08/27 15:23:30 tosa Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  VOL_GEOM vg = {}; // all values are initialized to be zero to detect change

  // argument handling
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv, vg) ;
    argc -= nargs ;
    argv += nargs ;
  }
  if (argc < 3)
  {
    print_usage();
    exit(-1);
  }

  char *invol = argv[argc-2];
  char *outvol = argv[argc-1];

  cerr << "Input volume  is : " << invol << endl;
  cerr << "Output volume is : " << outvol << endl;

  MRI *mri = MRIread(invol);
  if (!mri)
  {
    cerr << "could not open " << invol << endl;
    exit(-1);
  }
  //////////////////////////////////////////////
  VOL_GEOM vgIn;
  getVolGeom(mri, &vgIn);
  VOL_GEOM vgOut;
  copyVolGeom(&vgIn, &vgOut);
  // modify only those which are non-zero in the options
  // x_ras
  if (!FZERO(vg.x_r) || !FZERO(vg.x_a) || !FZERO(vg.x_s))
  {
    // check consistency
    if (!FZERO(vg.x_r*vg.x_r+vg.x_a*vg.x_a+vg.x_s*vg.x_s - 1))
    {
      cerr << "x_(ras) must have the unit length" << endl;
      exit(-1);
    }
    vgOut.x_r = vg.x_r;
    vgOut.x_a = vg.x_a;
    vgOut.x_s = vg.x_s;
  }
  // y_ras
  if (!FZERO(vg.y_r) || !FZERO(vg.y_a) || !FZERO(vg.y_s))
  {
    // check consistency
    if (!FZERO(vg.y_r*vg.y_r+vg.y_a*vg.y_a+vg.y_s*vg.y_s - 1))
    {
      cerr << "y_(ras) must have the unit length" << endl;
      exit(-1);
    }
    vgOut.y_r = vg.y_r;
    vgOut.y_a = vg.y_a;
    vgOut.y_s = vg.y_s;
  }
  // z_ras
  if (!FZERO(vg.z_r) || !FZERO(vg.z_a) || !FZERO(vg.z_s))
  {
    // check consistency
    if (!FZERO(vg.z_r*vg.z_r+vg.z_a*vg.z_a+vg.z_s*vg.z_s - 1))
    {
      cerr << "z_(ras) must have the unit length" << endl;
      exit(-1);
    }
    vgOut.z_r = vg.z_r;
    vgOut.z_a = vg.z_a;
    vgOut.z_s = vg.z_s;
  }
  // c_ras
  if (!FZERO(vg.c_r) || !FZERO(vg.c_a) || !FZERO(vg.c_s))
  {
    vgOut.c_r = vg.c_r;
    vgOut.c_a = vg.c_a;
    vgOut.c_s = vg.c_s;
  }
  // xsize
  if (!FZERO(vg.xsize))
    vgOut.xsize = vg.xsize;
  // ysize
  if (!FZERO(vg.ysize))
    vgOut.ysize = vg.ysize;
  // zsize
  if (!FZERO(vg.zsize))
    vgOut.zsize = vg.zsize;

  useVolGeomToMRI(&vgOut,mri); 

  // now TR, TE, TI, flip_angle
  if (gtr)
    mri->tr = gtr;
  if (gte)
    mri->te = gte;
  if (gti)
    mri->ti = gti;
  if (gflip_angle)
    mri->flip_angle = gflip_angle;

  // write it out
  MRIwrite(mri, outvol);

  return 0;
}
