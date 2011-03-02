/**
 * @file  mri_copy_params.cpp
 * @brief copy volume parameters from template and write out the volume
 *
 */
/*
 * Original Author: Yasunari Tosa
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:23 $
 *    $Revision: 1.5 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

extern "C" {
#include "error.h"
#include "mri.h"
#include "version.h"
#include "macros.h"
  const char *Progname = "mri_copy_params";
}
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_version(void) ;
static char vcid[] =
  "$Id: mri_copy_params.cpp,v 1.5 2011/03/02 00:04:23 nicks Exp $";

static int copy_pulse_params_only = 0 ;

using namespace std;

void print_usage() {
  cout << "Usage: mri_copy_params <in_vol> <template_vol> <out_vol>" << endl;
  cout << "     : where all volume parameters of in_vol are replaced with those of template_vol." << endl;
}

int main(int argc, char *argv[]) {
  bool bVolumeDifferent = false;
  bool bSizeDifferent = false;
  int ac, nargs;
  char **av ;
  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mri_copy_params.cpp,v 1.5 2011/03/02 00:04:23 nicks Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3) {
    print_usage();
    return -1;
  }

  MRI *in = MRIread(argv[1]);
  if (!in) {
    cerr << "could not open " << argv[1] << endl;
    return -1;
  }
  MRI *temp = MRIreadHeader(argv[2], MRI_UNDEFINED);
  if (!temp) {
    cerr << "could not open " << argv[2] << endl;
    return -1;
  }
  MRI *dst = MRIcopy(in, NULL);

  // check few things
  if ((temp->width != in->width)
      || (temp->height != in->height)
      || (temp->depth != in->depth)) {
    cerr << "WARNING: volume sizes are different" << endl;
    cerr << "    in_vol : " << in->width << ", " << in->height << ", " << in->depth << endl;
    cerr << "  temp_vol : " << temp->width << ", " << temp->height << ", " << temp->depth << endl;
    bVolumeDifferent = true;
  }
  if ((temp->xsize != in->xsize)
      || (temp->ysize != in->ysize)
      || (temp->zsize != in->zsize)) {
    cerr << "WARNING: voxel sizes are different" << endl;
    cerr << "    in_vol : " << in->xsize << ", " << in->ysize << ", " << in->zsize << endl;
    cerr << "  temp_vol : " << temp->xsize << ", " << temp->ysize << ", " << temp->zsize << endl;
    bSizeDifferent= true;
  }
  // copy everything in the header from template
  if (copy_pulse_params_only)
    MRIcopyPulseParameters(temp, dst) ;
  else
    MRIcopyHeader(temp, dst);
  // just few things restored
  if (bVolumeDifferent) {
    dst->width = in->width;
    dst->height = in->height;
    dst->depth = in->depth;
  }
  if (bSizeDifferent);
  {
    dst->xsize = in->xsize;
    dst->ysize = in->ysize;
    dst->zsize = in->zsize;
  }
  //
  MRIwrite(dst, argv[3]);

  MRIfree(&in);
  MRIfree(&temp);
  MRIfree(&dst);

  return (NO_ERROR);
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
    usage_exit() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "-pulse") || !stricmp(option, "-mri"))
  {
    printf("only copying pulse parameters\n") ;
    copy_pulse_params_only = 1 ;
  }
  else switch (toupper(*option))
  {
    
  }

  return(nargs) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}


