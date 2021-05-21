/*
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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "minc.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mrisurf.h"
#include "label.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void print_usage(void) ;

const char *Progname ;


#define NAME_LEN      100

static char subjects_dir[NAME_LEN] = "" ;

int
main(int argc, char *argv[]) {
  int          ac, nargs ;
  char         *cp, *subject_name, *hemi, *surface, *annot_file, *color_file, *output_file;
  MRI_SURFACE  *mris ;
  MRI          *mri ;
  int          err;
  COLOR_TABLE* ctab;
  char surf_name[NAME_LEN];
  char mri_name[NAME_LEN];
  int vno;
  VERTEX* v;
  int structure;
  float dx, dy, dz, len, d;
  double idxx, idxy, idxz;


  nargs = handleVersionOption(argc, argv, "mris_annot_to_segmentation");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  /* read in command-line options */
  ac = argc ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 7)
    print_usage() ;

  subject_name = argv[1];
  hemi         = argv[2];
  surface      = argv[3];
  annot_file   = argv[4];
  color_file   = argv[5];
  output_file  = argv[6];

  /* Read the surface first. */
  cp = getenv ("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "no subjects directory in environment.\n") ;
  strcpy (subjects_dir, cp) ;
  int req = snprintf(surf_name,NAME_LEN,"%s/%s/surf/%s.%s",subjects_dir,subject_name,hemi,surface);
  if( req >= NAME_LEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf (stderr, "reading %s...\n", surf_name) ;
  mris = MRISread (surf_name) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s\n",
              surf_name) ;

  MRISsaveVertexPositions (mris, TMP_VERTICES);
  err = MRISreadVertexPositions (mris, "pial");
  if (err != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read pial surface\n");
  MRISsaveVertexPositions (mris, CANONICAL_VERTICES);
  MRISrestoreVertexPositions (mris, TMP_VERTICES);

  /* Read the annotation file into the surface. */
  err = MRISreadAnnotation (mris, annot_file);
  if (err != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read annotation %s\n",
              annot_file);

  /* Read the color look up table. */
  ctab = CTABreadASCII ( color_file );
  if (NULL == ctab)
    ErrorExit(ERROR_NOFILE, "%s: could not read color table %s\n",
              color_file);

  /* Read in the T1 for this subject and change its name to the one
     they passed in. Set all values to 0. We'll use this as the
     segmentation volume. */
  req = snprintf (mri_name,NAME_LEN, "%s/%s/mri/T1",subjects_dir,subject_name);
  if( req >= NAME_LEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  mri = MRIread (mri_name);
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read T1 for template volume");
  MRIsetValues (mri, 0);

  /* For every annoted vertex... */
  for (vno = 0; vno < mris->nvertices; vno++) {

    v = &mris->vertices[vno];

    /* Get the color and then the index. */
    if ( 0 != v->annotation ) {
      structure = 0;
      if (mris->ct) {
        CTABfindAnnotation (mris->ct, v->annotation, &structure);
      } else {
        CTABfindAnnotation (ctab, v->annotation, &structure);
      }

      if (0 != structure) {

        /* Set all the voxels between the white (current) vertex
           and the pial vertex coords. Find the direction from the
           pial to the normal vertex positions. Then step from the
           normal outward in that direction in 0.1 increments to a
           distance of 1, convert the coord to ana idx, and set
           the voxel in the seg volume to the structure index we
           found above. */
        dx = v->cx - v->x;
        dy = v->cy - v->y;
        dz = v->cz - v->z;

        len = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (FZERO(len))
          len = 1.0 ;
        dx /= len;
        dy /= len;
        dz /= len;

        for ( d = 0 ; d <= len; d = d+0.1 ) {
          if (mris->useRealRAS) {
            MRIworldToVoxel (mri,
                             v->x + (d * dx),
                             v->y + (d * dy),
                             v->z + (d * dz),
                             &idxx, &idxy, &idxz);
          } else {
            MRIsurfaceRASToVoxel (mri,
                                  v->x + (d * dx),
                                  v->y + (d * dy),
                                  v->z + (d * dz),
                                  &idxx, &idxy, &idxz);
          }
          MRIvox(mri,(int)idxx,(int)idxy,(int)idxz) = (BUFTYPE)structure;
        }

      }
    }

    if ( !(vno % 1000) ) {
      fprintf( stdout, "\rConverting annotation... %.2f%% done",
               ((float)vno / (float)mris->nvertices) * 100.0 );
      fflush( stdout );
    }
  }

  fprintf( stdout, "\rConverting annotation... 100%% done       \n" );


  /* Write the segmentation */
  MRIwrite (mri, output_file);

  exit(0) ;
  return(0) ;
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
  switch (toupper(*option)) {
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
static void
print_usage(void) {
  printf("usage: %s <subject name> <hemi> <surface> <annot file> <color table> <output volume>\n", Progname);
  exit(1) ;
}
