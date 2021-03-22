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
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include "machine.h"
#include "fio.h"
#include "utils.h"
#include "mri.h"
#include "matrix.h"
#include "matfile.h"
#include "version.h"

#define IN_OUT_NAMES  100

const char *Progname;
int verbose_flag = 0;

int read_mat(int argc, char *argv[], int i, MATRIX *in_mat);
int write_mat(int argc, char *argv[], int i, MATRIX *in_mat);

char subject_name[STR_LEN];
char *subjnameuse=NULL;
float ipr, st, brightness;
int register_stuff_defined = 0;
int fsl_flag = 0;
int binarize = 0;

static void usage(int exit_val) {

  fprintf(stderr, "usage: %s [options]\n", Progname);
  fprintf(stderr, "options are:\n");
  fprintf(stderr, "  -v verbose\n");
  fprintf(stderr, "  -fsl : assume input/output are FSL-style matrix files\n");
  fprintf(stderr, "  -bin : 'binarize' output matrix.\n");
  fprintf(stderr, "  -s subject : use subject for subjectname in output reg.dat files\n");
  fprintf(stderr, "\n");
  fprintf(stderr, "'-im file' specifies input matrix files\n");
  fprintf(stderr, "'-iim file' specifies input matrix files to be inverted before multiplication\n");
  fprintf(stderr, "'-om file' specifies output matrix files\n");
  fprintf(stderr, "input and output files may be .dat or .xfm files\n");
  fprintf(stderr, "\n");
  fprintf(stderr, " mri_matrix_multiply -im M1.dat -im M2.dat -iim M3.dat -om M4.dat \n");
  fprintf(stderr, "    will compute M4 = M1*M2*inv(M3)\n");
  fprintf(stderr, "\n");
  exit(exit_val);

} /* end usage() */

/*----------------------------------------------------------------*/
int main(int argc, char *argv[]) {

  int i, r, c;
  int in_names[IN_OUT_NAMES], n_in = 0;
  int out_names[IN_OUT_NAMES], n_out = 0;
  MATRIX *in_mat, *result;
  int nargs;
  double v;

  nargs = handleVersionOption(argc, argv, "mri_matrix_multiply");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  /* ----- get the base executable name ----- */
  Progname = strrchr(argv[0], '/');
  Progname = (Progname == NULL ? argv[0] : Progname + 1);

  if (argc == 1)
    usage(1);

  for (i = 0;i < IN_OUT_NAMES;i++)
    in_names[i] = out_names[i] = 0;

  /* ----- get options ----- */
  for (i = 1;i < argc;i++) {
    if (strcmp(argv[i], "-im") == 0) {
      in_names[n_in] = i+1;
      n_in++;
      i++;
    } else if (strcmp(argv[i], "-iim") == 0) {
      in_names[n_in] = -(i+1);
      n_in++;
      i++;
    } else if (strcmp(argv[i], "-s") == 0) {
      subjnameuse = argv[i+1];
      i++;
    } else if (strcmp(argv[i], "-om") == 0) {
      out_names[n_out] = i+1;
      n_out++;
      i++;
    } else if (strcmp(argv[i], "-v") == 0) {
      verbose_flag = 1;
    } else if (strcmp(argv[i], "-fsl") == 0) {
      fsl_flag = 1;
    } else if (strcmp(argv[i], "-bin") == 0) {
      binarize = 1;
    } else {
      fprintf(stderr, "%s: unknown option %s\n", Progname, argv[i]);
    }
  }

  if (n_in == 0) {
    fprintf(stderr, "%s: no input files specified\n", Progname);
    usage(1);
  }

  if (n_out == 0) {
    fprintf(stderr, "%s: no output files specified\n", Progname);
    usage(1);
  }

  in_mat = MatrixAlloc(4, 4, MATRIX_REAL);
  result = MatrixIdentity(4,NULL);

  /* ----- read input files and keep a running product ----- */
  for (i = 0;i < n_in;i++) {
    if (read_mat(argc, argv, (in_names[i] < 0 ? -in_names[i] : in_names[i]), in_mat) == -1) {
      fprintf(stderr, "%s: exiting...\n", Progname);
      exit(1);
    }
    if (in_names[i] < 0) {

      if (verbose_flag)
        printf("inverting matrix\n");

      if (MatrixInverse(in_mat, in_mat) == NULL) {
        fprintf(stderr, "%s: couldn't invert matrix from file %s\n", Progname, argv[in_names[i]]);
        exit(0);
      }
    }
    MatrixMultiply(result, in_mat, result);
  }

  if (binarize) {
    /* "binarization" sets the rotational elements of the martrix to
       either +1, -1, or 0. +1 if the element is > 0.5. -1 if the element
       is < -0.5. 0 otherwise. */
    for (r=1;r<4;r++) {
      for (c=1;c<4;c++) {
        v = result->rptr[r][c];
        if (v > +0.5) result->rptr[r][c] = +1;
        else if (v < -0.5) result->rptr[r][c] = -1;
        else result->rptr[r][c] = 0.0;
      }
    }
  }

  /* ----- write output files ----- */
  for (i = 0;i < n_out;i++) {
    if (write_mat(argc, argv, out_names[i], result) == -1) {}
  }

  exit(0);

} /* end main() */

/*---------------------------------------------------------------*/
int read_mat(int argc, char *argv[], int i, MATRIX *in_mat) {

  FILE *fin;
  char line[STR_LEN];
  MATRIX *tmpmat;

  line[0] = '\0';

  if (i > argc) {
    fprintf(stderr, "%s: missing input matrix\n", Progname);
    return(-1);
  }

  if( ((strcmp(&argv[i][strlen(argv[i])-4], ".dat") == 0) ||
       (strcmp(&argv[i][strlen(argv[i])-4], ".reg") == 0) ) && !fsl_flag) {
    // tkregister-style

    if ((fin = fopen(argv[i], "r")) == NULL) {
      fprintf(stderr, "%s: error opening file %s\n", Progname, argv[i]);
      return(-1);
    }

    if (verbose_flag)
      printf("reading transform from .dat file %s\n", argv[i]);

    fscanf(fin, "%s", subject_name);
    fscanf(fin, "%f", &ipr);
    fscanf(fin, "%f", &st);
    fscanf(fin, "%f", &brightness);
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 1, 1),
           MATRIX_RELT(in_mat, 1, 2), MATRIX_RELT(in_mat, 1, 3), MATRIX_RELT(in_mat, 1, 4));
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 2, 1),
           MATRIX_RELT(in_mat, 2, 2), MATRIX_RELT(in_mat, 2, 3), MATRIX_RELT(in_mat, 2, 4));
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 3, 1), MATRIX_RELT(in_mat, 3, 2),
           MATRIX_RELT(in_mat, 3, 3), MATRIX_RELT(in_mat, 3, 4));
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 4, 1), MATRIX_RELT(in_mat, 4, 2),
           MATRIX_RELT(in_mat, 4, 3), MATRIX_RELT(in_mat, 4, 4));
    fclose(fin);
    register_stuff_defined = 1;
  } else if ((strcmp(&argv[i][strlen(argv[i])-4], ".xfm") == 0)  && !fsl_flag ) {
    // MINC-style XFM file
    if ((fin = fopen(argv[i], "r")) == NULL) {
      fprintf(stderr, "%s: error opening file %s\n", Progname, argv[i]);
      return(-1);
    }

    if (verbose_flag)
      printf("reading transform from .xfm file %s\n", argv[i]);

    while (strncmp(line, "Linear_Transform", 16) != 0) {
      if (fgets(line, STR_LEN, fin) == NULL) {
        fclose(fin);
        fprintf(stderr, "%s: premature EOF in file %s\n", Progname, argv[i]);
        return(-1);
      }
    }
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 1, 1), MATRIX_RELT(in_mat, 1, 2),
           MATRIX_RELT(in_mat, 1, 3), MATRIX_RELT(in_mat, 1, 4));
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 2, 1), MATRIX_RELT(in_mat, 2, 2),
           MATRIX_RELT(in_mat, 2, 3), MATRIX_RELT(in_mat, 2, 4));
    fscanf(fin, "%f %f %f %f;", MATRIX_RELT(in_mat, 3, 1), MATRIX_RELT(in_mat, 3, 2),
           MATRIX_RELT(in_mat, 3, 3), MATRIX_RELT(in_mat, 3, 4));
    fclose(fin);
  } else if (fsl_flag) {
    if ((fin = fopen(argv[i], "r")) == NULL) {
      fprintf(stderr, "%s: error opening file %s\n", Progname, argv[i]);
      return(-1);
    }

    if (verbose_flag)
      printf("reading matrix from fsl file %s\n", argv[i]);
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 1, 1), MATRIX_RELT(in_mat, 1, 2),
           MATRIX_RELT(in_mat, 1, 3), MATRIX_RELT(in_mat, 1, 4));
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 2, 1), MATRIX_RELT(in_mat, 2, 2),
           MATRIX_RELT(in_mat, 2, 3), MATRIX_RELT(in_mat, 2, 4));
    fscanf(fin, "%f %f %f %f;", MATRIX_RELT(in_mat, 3, 1), MATRIX_RELT(in_mat, 3, 2),
           MATRIX_RELT(in_mat, 3, 3), MATRIX_RELT(in_mat, 3, 4));
    fscanf(fin, "%f %f %f %f", MATRIX_RELT(in_mat, 4, 1), MATRIX_RELT(in_mat, 4, 2),
           MATRIX_RELT(in_mat, 4, 3), MATRIX_RELT(in_mat, 4, 4));
    fclose(fin);
  } else {
    /* try reading as a matlab file */
    tmpmat = MatlabRead(argv[i]);
    if (tmpmat == NULL) {
      printf("%s: unknown input matrix file type for file %s\n",
             Progname, argv[i]);
      return(-1);
    }
    if (verbose_flag) {
      printf("---------- %s ------------\n",argv[i]);
      MatrixPrint(stdout,tmpmat);
    }
    MatrixCopy(tmpmat,in_mat);
    MatrixFree(&tmpmat);
  }

  return(0);

} /* end read_mat() */

/*-----------------------------------------------------------------*/
int write_mat(int argc, char *argv[], int i, MATRIX *out_mat) {

  FILE *fout;

  if (i > argc) {
    fprintf(stderr, "%s: missing output matrix\n", Progname);
    return(-1);
  }

  if ((strcmp(&argv[i][strlen(argv[i])-4], ".dat") == 0) && !fsl_flag) {

    if (verbose_flag)
      printf("writing transform to .dat file %s\n", argv[i]);

    if ((fout = fopen(argv[i], "w")) == NULL) {
      fprintf(stderr, "%s: error opening file %s\n", Progname, argv[i]);
      return(-1);
    }

    if (!register_stuff_defined) {
      fprintf(stderr, "%s: extra parameters for .dat file %s undefined\n", Progname, argv[i]);
      fprintf(stderr, "%s: (not writing this file)\n", Progname);
      fclose(fout);
      return(-1);
    }

    if(subjnameuse == NULL) subjnameuse = subject_name;

    fprintf(fout, "%s\n", subjnameuse);
    fprintf(fout, "%f\n", ipr);
    fprintf(fout, "%f\n", st);
    fprintf(fout, "%f\n", brightness);
    fprintf(fout, "%f %f %f %f\n", *MATRIX_RELT(out_mat, 1, 1), *MATRIX_RELT(out_mat, 1, 2), *MATRIX_RELT(out_mat, 1, 3), *MATRIX_RELT(out_mat, 1, 4));
    fprintf(fout, "%f %f %f %f\n", *MATRIX_RELT(out_mat, 2, 1), *MATRIX_RELT(out_mat, 2, 2), *MATRIX_RELT(out_mat, 2, 3), *MATRIX_RELT(out_mat, 2, 4));
    fprintf(fout, "%f %f %f %f\n", *MATRIX_RELT(out_mat, 3, 1), *MATRIX_RELT(out_mat, 3, 2), *MATRIX_RELT(out_mat, 3, 3), *MATRIX_RELT(out_mat, 3, 4));
    fprintf(fout, "%f %f %f %f\n", *MATRIX_RELT(out_mat, 4, 1), *MATRIX_RELT(out_mat, 4, 2), *MATRIX_RELT(out_mat, 4, 3), *MATRIX_RELT(out_mat, 4, 4));
    fprintf(fout,"round\n");
    fclose(fout);

  } else if ((strcmp(&argv[i][strlen(argv[i])-4], ".xfm") == 0) && !fsl_flag) {

    if ((fout = fopen(argv[i], "w")) == NULL) {
      fprintf(stderr, "%s: error opening file %s\n", Progname, argv[i]);
      return(-1);
    }

    if (verbose_flag)
      printf("writing transform to .xfm file %s\n", argv[i]);

    fprintf(fout, "MNI Transform File\n");
    fprintf(fout, "%% Generated by %s\n", Progname);
    fprintf(fout, "\n");
    fprintf(fout, "Transform_Type = Linear;\n");
    fprintf(fout, "Linear_Transform =\n");
    fprintf(fout, "   %e %e %e %e\n",  *MATRIX_RELT(out_mat, 1, 1), *MATRIX_RELT(out_mat, 1, 2), *MATRIX_RELT(out_mat, 1, 3), *MATRIX_RELT(out_mat, 1, 4));
    fprintf(fout, "   %e %e %e %e\n",  *MATRIX_RELT(out_mat, 2, 1), *MATRIX_RELT(out_mat, 2, 2), *MATRIX_RELT(out_mat, 2, 3), *MATRIX_RELT(out_mat, 2, 4));
    fprintf(fout, "   %e %e %e %e;\n", *MATRIX_RELT(out_mat, 3, 1), *MATRIX_RELT(out_mat, 3, 2), *MATRIX_RELT(out_mat, 3, 3), *MATRIX_RELT(out_mat, 3, 4));

    fclose(fout);

  } else if (fsl_flag) {

    if ((fout = fopen(argv[i], "w")) == NULL) {
      fprintf(stderr, "%s: error opening file %s\n", Progname, argv[i]);
      return(-1);
    }

    if (verbose_flag)
      printf("writing transform to fsl file %s\n", argv[i]);

    fprintf(fout, "%f %f %f %f\n", *MATRIX_RELT(out_mat, 1, 1), *MATRIX_RELT(out_mat, 1, 2),
            *MATRIX_RELT(out_mat, 1, 3), *MATRIX_RELT(out_mat, 1, 4));
    fprintf(fout, "%f %f %f %f\n", *MATRIX_RELT(out_mat, 2, 1), *MATRIX_RELT(out_mat, 2, 2),
            *MATRIX_RELT(out_mat, 2, 3), *MATRIX_RELT(out_mat, 2, 4));
    fprintf(fout, "%f %f %f %f\n", *MATRIX_RELT(out_mat, 3, 1), *MATRIX_RELT(out_mat, 3, 2),
            *MATRIX_RELT(out_mat, 3, 3), *MATRIX_RELT(out_mat, 3, 4));
    fprintf(fout, "%f %f %f %f\n", *MATRIX_RELT(out_mat, 4, 1), *MATRIX_RELT(out_mat, 4, 2),
            *MATRIX_RELT(out_mat, 4, 3), *MATRIX_RELT(out_mat, 4, 4));

    fclose(fout);



  } else {
    fprintf(stderr, "%s: unknown output matrix file type for file %s\n", Progname, argv[i]);
    return(-1);
  }

  return(0);

} /* end read_mat() */

/* EOF */
