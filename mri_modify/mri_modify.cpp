/**
 * @brief modify direction cosine info on the volume.
 *
 * also allows changing the 'xform' filename
 */
/*
 * Original Author: Yasunari Tosa
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

 
#include "macros.h"
#include "mri.h"
#include "transform.h"
#include "version.h"

#include "compilerdefs.h"
  const char *Progname = "mri_modify";


using namespace std;

static double gflip_angle=0;
static float gte=0;
static float gtr=0;
static float gti=0;
static char new_transform_fname[STRLEN];

static int tr_specified = 0 ;
static int te_specified = 0 ;
static int ti_specified = 0 ;
static int fa_specified = 0 ;
static int cras_specified = 0 ;

void print_usage() {
  cout << "Usage: mri_modify <-xras xr xa xs> <-yras yr ya ys> <-zras zr za zs> <-cras cr ca cs> \\ " << endl;
  cout << "                  <-xsize size> <-ysize size> <-zsize size> \\ " << endl;
  cout << "                  <-tr recoverytime> <-te echotime> <-ti inversiontime> <-fa angledegree> \\ " << endl;
  cout << "                  <-xform new_file_name> \\ " << endl;
  cout << "                  involume outvolume" << endl;
}

int get_option(int argc, char *argv[], VOL_GEOM &vg) {
  int  nargs = 0 ;
  char *option ;
  option = argv[1] + 1 ;            /* past '-' */
  if (!strcmp(option, "-help")) {
    print_usage();
    exit(0);
  } else if (!stricmp(option, (char*)"xras")) {
    vg.x_r = atof(argv[2]);
    vg.x_a = atof(argv[3]);
    vg.x_s = atof(argv[4]);
    nargs=3;
  } else if (!stricmp(option, (char*)"yras")) {
    vg.y_r = atof(argv[2]);
    vg.y_a = atof(argv[3]);
    vg.y_s = atof(argv[4]);
    nargs=3;
  } else if (!stricmp(option, (char*)"zras")) {
    vg.z_r = atof(argv[2]);
    vg.z_a = atof(argv[3]);
    vg.z_s = atof(argv[4]);
    nargs=3;
  } else if (!stricmp(option, (char*)"cras")) {
    vg.c_r = atof(argv[2]);
    vg.c_a = atof(argv[3]);
    vg.c_s = atof(argv[4]);
    cras_specified = 1 ;
    nargs=3;
  } else if (!stricmp(option, (char*)"xsize")) {
    vg.xsize = atof(argv[2]);
    nargs=1;
  } else if (!stricmp(option, (char*)"ysize")) {
    vg.ysize = atof(argv[2]);
    nargs=1;
  } else if (!stricmp(option, (char*)"zsize")) {
    vg.zsize = atof(argv[2]);
    nargs=1;
  } else if (!stricmp(option, (char*)"tr")) {
    gtr=atof(argv[2]);
    tr_specified = 1 ;
    nargs=1;
  } else if (!stricmp(option, (char*)"te")) {
    gte=atof(argv[2]);
    te_specified = 1 ;
    nargs=1;
  } else if (!stricmp(option, (char*)"ti")) {
    gti=atof(argv[2]);
    ti_specified = 1 ;
    nargs=1;
  } else if (!stricmp(option, (char*)"fa")) {
    // mri stores it as radian
    gflip_angle=RADIANS(atof(argv[2]));
    fa_specified = 1 ;
    nargs=1;
  } else if (!stricmp(option, (char*)"xform")) {
    // get new transform file name
    strcpy(new_transform_fname,argv[2]);
  }
  return nargs;
}

int main(int argc, char *argv[]) {
  int nargs;
  nargs = handleVersionOption(argc, argv, "mri_modify");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmissing-field-initializers"
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wmissing-field-initializers"
#endif
  VOL_GEOM vg = {}; // all values are initialized to be zero to detect change
#if defined(FS_COMP_GNUC)
#pragma GCC diagnostic pop
#elif defined(FS_COMP_CLANG)
#pragma clang diagnostic pop
#endif

  new_transform_fname[0]=0; // null xform filename (assume no change)

  // argument handling
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv, vg) ;
    argc -= nargs ;
    argv += nargs ;
  }
  if (argc < 3) {
    print_usage();
    exit(-1);
  }

  char *invol = argv[argc-2];
  char *outvol = argv[argc-1];

  cerr << "Input volume  is : " << invol << endl;
  cerr << "Output volume is : " << outvol << endl;

  MRI *mri = MRIread(invol);
  if (!mri) {
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
  if (!FZERO(vg.x_r) || !FZERO(vg.x_a) || !FZERO(vg.x_s)) {
    // check consistency
    if (!FZERO(vg.x_r*vg.x_r+vg.x_a*vg.x_a+vg.x_s*vg.x_s - 1)) {
      cerr << "x_(ras) must have the unit length" << endl;
      exit(-1);
    }
    vgOut.x_r = vg.x_r;
    vgOut.x_a = vg.x_a;
    vgOut.x_s = vg.x_s;
  }
  // y_ras
  if (!FZERO(vg.y_r) || !FZERO(vg.y_a) || !FZERO(vg.y_s)) {
    // check consistency
    if (!FZERO(vg.y_r*vg.y_r+vg.y_a*vg.y_a+vg.y_s*vg.y_s - 1)) {
      cerr << "y_(ras) must have the unit length" << endl;
      exit(-1);
    }
    vgOut.y_r = vg.y_r;
    vgOut.y_a = vg.y_a;
    vgOut.y_s = vg.y_s;
  }
  // z_ras
  if (!FZERO(vg.z_r) || !FZERO(vg.z_a) || !FZERO(vg.z_s)) {
    // check consistency
    if (!FZERO(vg.z_r*vg.z_r+vg.z_a*vg.z_a+vg.z_s*vg.z_s - 1)) {
      cerr << "z_(ras) must have the unit length" << endl;
      exit(-1);
    }
    vgOut.z_r = vg.z_r;
    vgOut.z_a = vg.z_a;
    vgOut.z_s = vg.z_s;
  }
  // c_ras
  if (cras_specified){
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
  if (tr_specified)
  {
    printf("setting tr to %2.1f ms\n", gtr) ;
    mri->tr = gtr;
  }
  if (te_specified)
    mri->te = gte;
  if (ti_specified)
    mri->ti = gti;
  if (fa_specified)
    mri->flip_angle = gflip_angle;

  // stuff-in the new transform filename, if one was grabbed from command-line
  if (new_transform_fname[0])
    strcpy(mri->transform_fname,new_transform_fname);

  // write it out
  MRIwrite(mri, outvol);

  return 0;
}
