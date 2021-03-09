/**
 * @brief copy volume parameters from template and write out the volume
 *
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

 
#include "error.h"
#include "mri.h"
#include "version.h"
#include "macros.h"
  const char *Progname = "mri_copy_params";

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_version(void) ;

static int copy_pulse_params_only = 0 ;
static int copy_ras_only = 0 ;
static int copy_voxel_size = 0 ;

using namespace std;

void print_usage() {
  cout << "Usage: mri_copy_params <in_vol> <template_vol> <out_vol>" << endl;
  cout << "     : where all volume parameters of in_vol are replaced with those of template_vol." << endl;
  cout << "use --size to force copying of voxel sizes when resolutions var" << endl;
}

int main(int argc, char *argv[]) {
  bool bVolumeDifferent = false;
  bool bSizeDifferent = false;
  int  nargs;
  nargs = handleVersionOption(argc, argv, "mri_copy_params");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

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
  {
    if (copy_ras_only)
    {
      dst->x_r = temp->x_r;
      dst->x_a = temp->x_a;
      dst->x_s = temp->x_s;
      dst->y_r = temp->y_r;
      dst->y_a = temp->y_a;
      dst->y_s = temp->y_s;
      dst->z_r = temp->z_r;
      dst->z_a = temp->z_a;
      dst->z_s = temp->z_s;
      dst->c_r = temp->c_r;
      dst->c_a = temp->c_a;
      dst->c_s = temp->c_s;
      dst->ras_good_flag = temp->ras_good_flag;
      dst->i_to_r__ = AffineMatrixCopy( temp->i_to_r__,
					dst->i_to_r__ );
      
      dst->r_to_i__ = MatrixCopy(temp->r_to_i__, dst->r_to_i__);
    }
    else
      MRIcopyHeader(temp, dst);
  }
  // just few things restored
  if (bVolumeDifferent) {
    dst->width = in->width;
    dst->height = in->height;
    dst->depth = in->depth;
  }
  if (bSizeDifferent)
  {
    if (copy_voxel_size)
    {
      printf("using template resolution\n");
      dst->xsize = temp->xsize;
      dst->ysize = temp->ysize;
      dst->zsize = temp->zsize;
    }
    else
    {
      printf("retaining input resolution even though voxel sizes vary\n");
      dst->xsize = in->xsize;
      dst->ysize = in->ysize;
      dst->zsize = in->zsize;
    }
    printf("setting destination sizes to (%2.2f, %2.2f %2.2f)\n",
	   dst->xsize, dst->ysize, dst->zsize) ;
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
  else if (!stricmp(option, "-ras"))
  {
    printf("only copying ras2vox matrices\n") ;
    copy_ras_only = 1 ;
  }
  else if (!stricmp(option, "-size"))
  {
    printf("only copying voxel sizes\n") ;
    copy_voxel_size = 1 ;
  }
  else switch (toupper(*option))
  {
    
  }

  return(nargs) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}
static void
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}


