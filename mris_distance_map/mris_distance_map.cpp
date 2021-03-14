/**
 * @brief program for computing a distance map on the surface
 *
 * compute the distance of each point on the surface to a reference point
 * (usually vertex 0).
 */
/*
 * Original Author: Bruce Fischl
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

#include "diag.h"
#include "mrisutils.h"
#include "version.h"

//------------------------------------------------------------------------

/*-------------------------------- CONSTANTS -----------------------------*/

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]);

static int  get_option(int argc, char *argv[]);
static void usage_exit();
static void print_usage();
static void print_help();
static void print_version();
MRI *       MRIScomputeDistanceMap(MRI_SURFACE *mris, MRI *mri_distance,
                                   int ref_vertex_no);

/*-------------------------------- DATA ----------------------------*/

const char *Progname;

static int ref_vertex_no = 0;

/*-------------------------------- FUNCTIONS ----------------------------*/

int main(int argc, char *argv[]) {
  MRI_SURFACE *mris;
  char **      av, *in_fname, *out_fname;
  int          ac, nargs;
  MRI *        mri_distance;

  nargs = handleVersionOption(argc, argv, "mris_distance_map");
  if (nargs && argc - nargs == 1)
    exit(0);
  argc -= nargs;

  Progname = argv[0];
  ErrorInit(NULL, NULL, NULL);
  DiagInit(nullptr, nullptr, nullptr);

  ac = argc;
  av = argv;
  for (; argc > 1 && ISOPTION(*argv[1]); argc--, argv++) {
    nargs = get_option(argc, argv);
    argc -= nargs;
    argv += nargs;
  }

  if (argc < 3)
    usage_exit();

  in_fname  = argv[1];
  out_fname = argv[2];
  mris      = MRISread(in_fname);
  if (mris == nullptr)
    ErrorExit(ERROR_NOFILE, "%s: could not load surface %s", Progname,
              out_fname);

  mri_distance = MRIScomputeDistanceMap(mris, nullptr, ref_vertex_no);

  MRIwrite(mri_distance, out_fname);
  MRISfree(&mris);

  exit(0);

  return (0); /* for ansi */
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int get_option(int argc, char *argv[]) {
  int   nargs = 0;
  char *option;

  option = argv[1] + 1; /* past '-' */
  if (!stricmp(option, "-help"))
    print_help();
  else if (!stricmp(option, "-version"))
    print_version();
  else
    switch (toupper(*option)) {
    case 'V':
      ref_vertex_no = atoi(argv[2]);
      nargs         = 1;
      break;
    case '?':
    case 'U':
    case 'H':
      print_help();
      exit(1);
      break;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]);
      exit(1);
      break;
    }

  return (nargs);
}

static void usage_exit() {
  print_help();
  exit(1);
}

static void print_usage() {
  fprintf(
      stderr,
      "usage: %s [options] <input surface file> <output scalar field (.mgz)>\n",
      Progname);
}

static void print_help() {
  print_usage();
  printf("\nThis program will compute a distance map of each point on the "
         "surface to a\n"
         "reference point (vertex 0 by default)\n");
  printf("\nvalid options are:\n");
  exit(1);
}

static void print_version(void) {
  printf("%s\n", getVersion().c_str());
  exit(1);
}
MRI *MRIScomputeDistanceMap(MRI_SURFACE *mris, MRI *mri_distance,
                            int ref_vertex_no) {
  int     vno;
  VERTEX *v;
  double  circumference, angle, distance;
  VECTOR *v1, *v2;

  if (mri_distance == nullptr)
    mri_distance = MRIalloc(mris->nvertices, 1, 1, MRI_FLOAT);

  v1 = VectorAlloc(3, MATRIX_REAL);
  v2 = VectorAlloc(3, MATRIX_REAL);
  v  = &mris->vertices[ref_vertex_no];
  VECTOR_LOAD(v1, v->x, v->y, v->z); /* radius vector */
  circumference = M_PI * 2.0 * V3_LEN(v1);
  for (vno = 0; vno < mris->nvertices; vno++) {
    v = &mris->vertices[vno];
    if (vno == Gdiag_no)
      DiagBreak();
    VECTOR_LOAD(v2, v->x, v->y, v->z); /* radius vector */
    angle    = fabs(Vector3Angle(v1, v2));
    distance = circumference * angle / (2.0 * M_PI);
    MRIsetVoxVal(mri_distance, vno, 0, 0, 0, distance);
  }

  VectorFree(&v1);
  VectorFree(&v2);
  return (mri_distance);
}
