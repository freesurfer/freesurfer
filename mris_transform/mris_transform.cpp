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

#include "diag.h"
#include "mrisurf.h"
#include "version.h"

int main(int argc, char *argv[]);

static int  get_option(int argc, char *argv[]);
static void print_usage_exit();
static void print_version_exit();
static void write_surface_tidy_up(MRI_SURFACE *mris, const char *fname);

const char *Progname;
static int  inverse_flag = 0;

MRI *mri_src = nullptr;
MRI *mri_dst = nullptr;

int main(int argc, char *argv[]) {
  char **av, *in_fname, *out_fname, *xform_fname;
  int    ac, nargs;

  if (argc == 1)
    print_usage_exit();

  nargs = handleVersionOption(argc, argv, "mris_transform");
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

  if (argc != 4)
    ErrorExit(ERROR_BADPARM, "ERROR: incorrect number of arguments");

  in_fname    = argv[1];
  xform_fname = argv[2];
  out_fname   = argv[3];

  MRI_SURFACE *mris = MRISread(in_fname);
  if (!mris)
    ErrorExit(ERROR_NOFILE, "ERROR: could not read surface file %s", in_fname);

  if (!strcmp(xform_fname, "identity.nofile")) {
    fprintf(stdout, "Transform is identity.nofile. Copying surface to %s.\n",
            out_fname);
    write_surface_tidy_up(mris, out_fname);
    return (EXIT_SUCCESS);
  }

  TRANSFORM *transform = TransformRead(xform_fname);
  if (!transform)
    ErrorExit(ERROR_NOFILE, "ERROR: could not read transform file %s",
              xform_fname);

  if (transform->type != MORPH_3D_TYPE) {
    LTA *tmp = (LTA *)transform->xform;
    LTA *lta = LTAreduce(tmp); // Apply full array, allocation.
    LTAfree(&tmp);
    if (lta->xforms[0].src.valid == 0) {
      if (mri_src == nullptr) {
        fprintf(stderr,
                "The transform does not have the valid src volume info.\n");
        fprintf(stderr, "Either you specify it with --trx-src or\n");
        fprintf(stderr, "make the transform to have the valid src info.\n");
        ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
      }
      getVolGeom(mri_src, &lta->xforms[0].src);
    }
    if (lta->xforms[0].dst.valid == 0) {
      if (mri_dst == nullptr) {
        fprintf(stderr,
                "The transform does not have the valid dst volume info.\n");
        fprintf(stderr, "Either you specify it with --trx-dst or\n");
        fprintf(stderr, "make the transform to have the valid dst info.\n");
        fprintf(stderr, "If the dst was average_305, then you can set\n");
        fprintf(stderr, "environmental variable USE_AVERAGE305 true\n");
        fprintf(stderr,
                "without giving the dst volume for RAS-to-RAS transform.\n");
        ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
      }
      getVolGeom(mri_dst, &lta->xforms[0].dst);
    }
    LTAchangeType(lta, LINEAR_VOX_TO_VOX); // Support more types.
    transform->type  = LINEAR_VOX_TO_VOX;
    transform->xform = (void *)lta;
  }

  // To map source to target image, GCAMs contain a coordinate transform
  // from target to source: need to invert GCAMs as we want to transform
  // coordinates from source to target here.
  int do_invert = transform->type == MORPH_3D_TYPE;
  if (inverse_flag)
    do_invert = !do_invert;
  if (do_invert)
    TransformInvertReplace(transform, mri_dst);

  if (!transform->xform) {
    fprintf(stderr, "ERROR: could not invert transform. Try explicitly");
    fprintf(stderr, " specifying the target geometry of the transform");
    fprintf(stderr, " if the target volume was moved.\n");
    exit(EXIT_FAILURE);
  }

  // MRIStransform() interprets source/target MRIs differently for LTAs and
  // GCAMs. If NULL is passed, it will just figure it out from the transform.
  MRIStransform(mris, nullptr, transform, nullptr);
  write_surface_tidy_up(mris, out_fname);
  return (EXIT_SUCCESS);
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int get_option(int argc, char *argv[]) {
  int   nargs = 0;
  char *option;

  option = argv[1] + 1; /* past '-' */
  if (!stricmp(option, "-help") || !stricmp(option, "h") ||
      !stricmp(option, "u") || !stricmp(option, "?"))
    print_usage_exit();
  else if (!stricmp(option, "-trx-src") || !stricmp(option, "s")) {
    fprintf(stderr, "Reading src volume of transform...\n");
    mri_src = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_src) {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  } else if (!stricmp(option, "-trx-dst") || !stricmp(option, "d")) {
    fprintf(stderr, "Reading dst volume of transform...\n");
    mri_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_dst) {
      ErrorExit(ERROR_BADPARM, "ERROR: could not read file %s\n", argv[2]);
    }
    nargs = 1;
  } else if (!stricmp(option, "-version") || !stricmp(option, "v"))
    print_version_exit();
  else if (!stricmp(option, "-is-inverse") || !stricmp(option, "i"))
    inverse_flag = 1;
  else
    ErrorExit(ERROR_BADPARM, "ERROR: unknown option %s", argv[1]);

  return (nargs);
}

#include "mris_transform.help.xml.h"
static void print_usage_exit() {
  outputHelpXml(mris_transform_help_xml, mris_transform_help_xml_len);
  exit(EXIT_SUCCESS);
}

static void print_version_exit(void) {
  fprintf(stderr, "%s\n", getVersion().c_str());
  exit(EXIT_SUCCESS);
}

static void write_surface_tidy_up(MRI_SURFACE *mris, const char *fname) {
  MRISwrite(mris, fname);
  MRISfree(&mris);
  if (mri_src)
    MRIfree(&mri_src);
  if (mri_dst)
    MRIfree(&mri_dst);
}
