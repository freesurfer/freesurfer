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


/* Convert an transformation matrix from flirt output to lta format */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "fio.h"
#include "version.h"
#include "transform.h"

void usage(int exit_val);

static LTA  *ltaFSLread(const char *fname) ;

const char *Progname;

static int get_option(int argc, char *argv[]) ;

static int invert_flag = 0;

int main(int argc, char *argv[]) {
  char **av;
  char *fslfn, *ltafn;
  MRI *mri_src, *mri_tgt;
  int ac, nargs;
  LTA *mylta = 0;
  FILE *fo;
  VOL_GEOM srcG, dstG;
  MATRIX *invTgt, *Dsrc, *V_to_V;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_fslmat_to_lta");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc != 5)
    usage(1);

  printf("Read source volume file %s\n", argv[1]);
  mri_src = MRIread(argv[1]) ;

  if (!mri_src)
    ErrorExit(ERROR_BADPARM, "%s: could not read source volume %s",
              Progname, argv[1]) ;

  printf("Read destination volume file %s\n", argv[2]);
  mri_tgt = MRIread(argv[2]) ;
  if (!mri_tgt)
    ErrorExit(ERROR_BADPARM, "%s: could not read label volume %s",
              Progname, argv[2]) ;

  fslfn = argv[3];
  ltafn = argv[4];

  printf("Read fsl transformation from file %s\n", fslfn);
  mylta = ltaFSLread(fslfn);

  getVolGeom(mri_src, &srcG);
  getVolGeom(mri_tgt, &dstG);

  if (invert_flag) {
    mylta->xforms[0].src = dstG;
    mylta->xforms[0].dst = srcG;
  } else {
    mylta->xforms[0].src = srcG;
    mylta->xforms[0].dst = dstG;
  }

  invTgt = MatrixAlloc(4,4,MATRIX_REAL);
  invTgt->rptr[1][1] = 1.0/dstG.xsize;
  invTgt->rptr[2][2] = 1.0/dstG.ysize;
  invTgt->rptr[3][3] = 1.0/dstG.zsize;
  invTgt->rptr[4][4] = 1.0;

  Dsrc = MatrixAlloc(4,4,MATRIX_REAL);
  Dsrc->rptr[1][1] = srcG.xsize;
  Dsrc->rptr[2][2] = srcG.ysize;
  Dsrc->rptr[3][3] = srcG.zsize;
  Dsrc->rptr[4][4] = 1.0;

  V_to_V = MatrixMultiply(invTgt, mylta->xforms[0].m_L, NULL);
  V_to_V = MatrixMultiply(V_to_V, Dsrc, V_to_V);

  if (invert_flag) {
    mylta->xforms[0].m_L = MatrixInverse(V_to_V, mylta->xforms[0].m_L);
  } else {
    mylta->xforms[0].m_L = MatrixCopy(V_to_V, mylta->xforms[0].m_L);
  }

  if (!mylta) {
    MRIfree(&mri_src);
    MRIfree(&mri_tgt);
    ErrorExit(ERROR_BADFILE, "%s: can't read in transformation",Progname);
  }

  printf("Write transformation to lta file %s\n", ltafn);
  fo = fopen(ltafn,"w");
  if (fo==NULL)
    ErrorExit(ERROR_BADFILE, "%s: can't create file %s",Progname, ltafn);

  LTAprint(fo, mylta);

  fclose(fo);


  MRIfree(&mri_src);
  MRIfree(&mri_tgt);

  LTAfree(&mylta);

  MatrixFree(&invTgt);
  MatrixFree(&Dsrc);
  MatrixFree(&V_to_V);

  return(0);

}  /*  end main()  */

void usage(int exit_val) {

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <src vol> <target vol> <fslmat_file> <lta_file> \n", Progname);
  fprintf(fout, "this program creates the lta-transformation file using \n") ;
  fprintf(fout, "information of src and target volumes and fslmat_file \n") ;

  exit(exit_val);

}  /*  end usage()  */

static LTA  *
ltaFSLread(const char *fname) {
  LTA              *lta ;
  LINEAR_TRANSFORM *lt ;
  char             *cp, line[1000] ;
  FILE             *fp ;
  int              row ;
  MATRIX           *m_L ;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NOFILE, "ltFSLread: could not open file %s",fname));

  lta = LTAalloc(1, NULL) ;
  lt = &lta->xforms[0] ;
  lt->sigma = 1.0f ;
  lt->x0 = lt->y0 = lt->z0 = 0 ;

  m_L = lt->m_L ;
  for (row = 1 ; row <= 3 ; row++) {
    cp = fgetl(line, 900, fp) ;
    if (!cp) {
      LTAfree(&lta) ;
      ErrorReturn(NULL,
                  (ERROR_BADFILE, "ltFSLread: could not read row %d from %s",
                   row, fname)) ;
    }
    sscanf(cp, "%f %f %f %f",
           MATRIX_RELT(m_L,row,1), MATRIX_RELT(m_L,row,2),
           MATRIX_RELT(m_L,row,3), MATRIX_RELT(m_L,row,4)) ;
  }
  fclose(fp) ;
  lta->type = LINEAR_VOX_TO_VOX;
  return(lta) ;
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

  if (!stricmp(option, "inverse")) {
    invert_flag = 1;
    printf("convert to lta and then compute its inverse \n");
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      usage(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

/*  EOF  */
