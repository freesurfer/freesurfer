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


/* Concatenate two (or more?) LTAs into one final LTA */

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

static int get_option(int argc, char *argv[]);

LTA *ltaReadFileEx(const char *fname);
LTA *ltaMNIreadEx(const char *fname);
int  ltaMNIwrite(LTA *lta, char *fname);

static int invert1 = 0;
static int invert2 = 0;

static int out_type = 0;

const char *Progname;
char *tal_src_file = 0;
char *tal_dst_file = 0;
MRI *tal_src = 0;
MRI *tal_dst = 0;

int main(int argc, char *argv[])
{

  char **av, *ltafn1, *ltafn2, *ltafn_total;
  LTA *lta1, *lta2, *lta_total;
  FILE *fo;
  MATRIX *r_to_i_1, *i_to_r_1, *i_to_r_2, *r_to_i_2;
  MATRIX *RAS_1_to_1, *RAS_2_to_2, *m_tmp;
  int nargs, ac;
  int type = 0;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_concatenate_lta");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }


  if (argc != 4)
    usage(1);

  ltafn1 = argv[1];
  ltafn2 = argv[2];
  ltafn_total = argv[3];

  printf("Read individual LTAs\n");
  lta1 = ltaReadFileEx(ltafn1);
  if (!lta1)
    ErrorExit(ERROR_BADFILE, "%s: can't read file %s",Progname, ltafn1);

  if (invert1)
  {
    VOL_GEOM vgtmp;
    LT *lt;
    MATRIX *m_tmp = lta1->xforms[0].m_L ;
    lta1->xforms[0].m_L = MatrixInverse(lta1->xforms[0].m_L, NULL) ;
    MatrixFree(&m_tmp) ;
    lt = &lta1->xforms[0];
    if (lt->dst.valid == 0 || lt->src.valid == 0)
    {
      fprintf(stderr, "WARNING:***************************************************************\n");
      fprintf(stderr, "WARNING:dst volume infor is invalid.  Most likely produce wrong inverse.\n");
      fprintf(stderr, "WARNING:***************************************************************\n");
    }
    copyVolGeom(&lt->dst, &vgtmp);
    copyVolGeom(&lt->src, &lt->dst);
    copyVolGeom(&vgtmp, &lt->src);
  }

  type = TransformFileNameType(ltafn2);
  if (type == MNI_TRANSFORM_TYPE)
  {

    if (invert2 != 0)
      ErrorExit(ERROR_BADFILE, "%s: LTA2 is talairach.xfm, and shouldn't be inverted ",Progname);

    lta2 = ltaMNIreadEx(ltafn2) ;
    //the talairach xform is supposed to be linear_RAS_TO_RAS, right? Yes
    lta2->type =  LINEAR_RAS_TO_RAS;

    if (tal_src_file == 0 && lta2->xforms[0].src.valid == 0)
      ErrorExit(ERROR_BADFILE, "%s: pls use -tal option to give talairach src and template filenames",Progname);
    if (tal_dst_file == 0 && lta2->xforms[0].dst.valid == 0)
      ErrorExit(ERROR_BADFILE, "%s: pls use -tal option to give talairach src and template filenames",Progname);

    if (tal_src_file != 0)
      LTAmodifySrcDstGeom(lta2, tal_src, NULL); // add src and dst information
    if (tal_dst_file != 0)
      LTAmodifySrcDstGeom(lta2, NULL, tal_dst); // add src and dst information
  }
  else
    lta2 = ltaReadFileEx(ltafn2);

  if (!lta2)
    ErrorExit(ERROR_BADFILE, "%s: can't read file %s",Progname, ltafn2);

  if (invert2)
  {
    VOL_GEOM vgtmp;
    LT *lt;
    MATRIX *m_tmp = lta2->xforms[0].m_L ;
    lta2->xforms[0].m_L = MatrixInverse(lta2->xforms[0].m_L, NULL) ;
    MatrixFree(&m_tmp) ;
    lt = &lta2->xforms[0];
    if (lt->dst.valid == 0 || lt->src.valid == 0)
    {
      fprintf(stderr, "WARNING:***************************************************************\n");
      fprintf(stderr, "WARNING:dst volume infor is invalid.  Most likely produce wrong inverse.\n");
      fprintf(stderr, "WARNING:***************************************************************\n");
    }
    copyVolGeom(&lt->dst, &vgtmp);
    copyVolGeom(&lt->src, &lt->dst);
    copyVolGeom(&vgtmp, &lt->src);
  }

  if (vg_isEqual(&lta1->xforms[0].dst, &lta2->xforms[0].src) == 0)
  {
    /*    ErrorExit(ERROR_BADFILE, "%s: dst volume of lta1 doesn't match src volume of lta2",Progname);*/
    printf("Warning: dst volume of lta1 doesn't match src volume of lta2\n");
    printf("Volume geometry for lta1-dst:\n");
    vg_print(&lta1->xforms[0].dst);
    printf("Volume geometry for lta2-src:\n");
    vg_print(&lta2->xforms[0].src);
  }

  printf("Combine the two LTAs to get a RAS-to-RAS from src of LTA1 to dst of LTA2\n");

  if (lta1->type == LINEAR_RAS_TO_RAS)
  {
    RAS_1_to_1 =  MatrixCopy(lta1->xforms[0].m_L, NULL);
  }
  else if (lta1->type == LINEAR_VOX_TO_VOX)
  {
    r_to_i_1 = vg_r_to_i(&lta1->xforms[0].src);
    i_to_r_1 = vg_i_to_r(&lta1->xforms[0].dst);
    if (!r_to_i_1 || !i_to_r_1)
      ErrorExit(ERROR_BADFILE, "%s: failed to convert LTA1 to RAS_to_RAS",Progname);
    m_tmp = MatrixMultiply(lta1->xforms[0].m_L, r_to_i_1, NULL);
    RAS_1_to_1 = MatrixMultiply(i_to_r_1, m_tmp, NULL);
    MatrixFree(&m_tmp);
  }
  else
  {
    ErrorExit(ERROR_BADFILE, "%s: unknown transform type for LTA1",Progname);
  }

  if (lta2->type == LINEAR_RAS_TO_RAS)
  {
    RAS_2_to_2 =  MatrixCopy(lta2->xforms[0].m_L, NULL);
  }
  else if (lta2->type == LINEAR_VOX_TO_VOX)
  {
    r_to_i_2 = vg_r_to_i(&lta2->xforms[0].src);
    i_to_r_2 = vg_i_to_r(&lta2->xforms[0].dst);
    if (!r_to_i_2 || !i_to_r_2)
      ErrorExit(ERROR_BADFILE, "%s: failed to convert LTA1 to RAS_to_RAS",Progname);
    m_tmp = MatrixMultiply(lta2->xforms[0].m_L, r_to_i_2, NULL);
    RAS_2_to_2 = MatrixMultiply(i_to_r_2, m_tmp, NULL);
    MatrixFree(&m_tmp);
  }
  else
  {
    ErrorExit(ERROR_BADFILE, "%s: unknown transform type for LTA1",Progname);
  }

  lta_total = LTAalloc(1, NULL);
  lta_total->type = LINEAR_RAS_TO_RAS;
  MatrixMultiply(RAS_2_to_2, RAS_1_to_1, lta_total->xforms[0].m_L);
  lta_total->xforms[0].src = lta1->xforms[0].src;
  lta_total->xforms[0].dst = lta2->xforms[0].dst;
  lta_total->xforms[0].x0 = 0;
  lta_total->xforms[0].y0 = 0;
  lta_total->xforms[0].z0 = 0;
  lta_total->xforms[0].sigma = 1.0f;

  type = TransformFileNameType(ltafn_total);
  if (type == MNI_TRANSFORM_TYPE)
  {
    ltaMNIwrite(lta_total, ltafn_total);
  }
  else
  {
    //change type to VOXEL_VOXEL
    if (lta_total->type != out_type)
      LTAchangeType(lta_total, out_type);

    printf("Write combined LTA to file %s\n", ltafn_total);
    fo = fopen(ltafn_total,"w");
    if (fo==NULL)
      ErrorExit(ERROR_BADFILE, "%s: can't create file %s",Progname, ltafn_total);

    LTAprint(fo, lta_total);

    fclose(fo);
  }

  LTAfree(&lta1);
  LTAfree(&lta2);
  LTAfree(&lta_total);
  MatrixFree(&RAS_1_to_1);
  MatrixFree(&RAS_2_to_2);

  if (tal_src) MRIfree(&tal_src);
  if (tal_dst) MRIfree(&tal_dst);

  return(0);

}  /*  end main()  */

void usage(int exit_val)
{

  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s  <lta_1> <lta_2> <lta_final> \n", Progname);
  fprintf(fout, "this program concatenate two consecutive LTA transformations \n") ;
  fprintf(fout, "into one overall transformation. \n") ;
  fprintf(fout, "lta_1 maps src1 to dst1, lta_2 maps dst1(src2) to dst2 \n") ;
  fprintf(fout, "The combined LTA maps src1 to dst2 \n") ;
  fprintf(fout, "If lta2 is talairach.xfm, use -tal file1 file2 to specify src (file1) and template (file2) for the talairach xfm \n") ;
  exit(exit_val);

}  /*  end usage()  */

LTA *ltaReadFileEx(const char *fname)
{
  FILE             *fp;
  LINEAR_TRANSFORM *lt ;
  int              i, nxforms, type ;
  char             line[STRLEN], *cp ;
  LTA              *lta ;

  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorReturn(NULL,
                (ERROR_BADFILE, "ltaReadFile(%s): can't open file",fname));
  cp = fgetl(line, 199, fp) ;
  if (cp == NULL)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "ltaReadFile(%s): can't read data",fname));
  }
  sscanf(cp, "type      = %d\n", &type) ;
  cp = fgetl(line, 199, fp) ;
  sscanf(cp, "nxforms   = %d\n", &nxforms) ;
  lta = LTAalloc(nxforms, NULL) ;
  lta->type = type ;
  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    fscanf(fp, "mean      = %f %f %f\n", &lt->x0, &lt->y0, &lt->z0) ;
    fscanf(fp, "sigma     = %f\n", &lt->sigma) ;
    MatrixAsciiReadFrom(fp, lt->m_L) ;
  }
  // oh, well this is the added part
  for (i=0; i < lta->num_xforms; i++)
  {
    if (fgets(line, 199, fp))
    {
      if (strncmp(line, "src volume info", 15)==0)
      {
        char *p;
        readVolGeom(fp, &lta->xforms[i].src);
        p = fgets(line, 199, fp);
        if (strncmp(line, "dst volume info", 15)==0)
          readVolGeom(fp, &lta->xforms[i].dst);
      }
    }
  }
  fclose(fp) ;
  return(lta) ;
}

/*  EOF  */


static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "help"))
    usage(0) ;
  else if (!stricmp(option, "invert1"))
  {
    invert1 = 1;
    fprintf(stderr, "invert the first LTA before applying it \n");
  }
  else if (!stricmp(option, "out_type"))
  {
    out_type = atoi(argv[2]) ;
    nargs = 1;
    fprintf(stderr, "set final LTA type to %d\n", out_type);
  }
  else if (!stricmp(option, "tal"))
  {
    tal_src_file = argv[2];
    tal_dst_file = argv[3];
    nargs = 2;
    fprintf(stderr, "Talairach xfrm src file is %s\n", tal_src_file);
    fprintf(stderr, "Talairach xfrm dst file is %s\n", tal_dst_file);
    tal_src = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!tal_src)
    {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    tal_dst = MRIreadHeader(argv[3], MRI_VOLUME_TYPE_UNKNOWN);
    if (!tal_dst)
    {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[3]);
    }
  }
  else if (!stricmp(option, "invert2"))
  {
    invert2 = 1;
    fprintf(stderr, "invert the second LTA before applying it \n");
  }
  else
  {
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    usage(1) ;
    exit(1) ;
  }

  return(nargs) ;
}


LTA *ltaMNIreadEx(const char *fname)
{
  LTA *lta = 0;
  LINEAR_TRANSFORM *lt ;
  char             *cp, line[1000], infoline[1024], infoline2[1024];
  FILE             *fp ;
  int              row ;
  MATRIX           *m_L ;
  int             no_volinfo = 0;

  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorReturn(NULL,
                (ERROR_NOFILE, "ltMNIreadEx: could not open file %s",fname));

  lta = LTAalloc(1, NULL) ;
  lt = &lta->xforms[0] ;
  lt->sigma = 1.0f ;
  lt->x0 = lt->y0 = lt->z0 = 0 ;

  fgetl(line, 900, fp) ;   /* MNI Transform File */
  if (strncmp("MNI Transform File", line, 18))
    ErrorReturn(NULL,
                (ERROR_NOFILE, "ltMNIreadEx:%s does not start as 'MNI Transform File'",fname));

  fgetl(line, 900, fp) ;   /* fileinfo line */
  if (line[0] == '%')
    strcpy(infoline, line);
  else
  {
    no_volinfo = 1;
    if (!strncmp("Transform_Type", line, 14))
    {
      fgetl(line,900,fp);
      goto get_transform;
    }
  }
  // second line in %
  fgetl(line, 900, fp);
  if (line[0] == '%')
  {
    strcpy(infoline2, line);
    while (line[0] == '%')
      fgetl(line, 900, fp) ; /* variable # of comments */
    fgetl(line, 900, fp) ;
    if (!strncmp("Transform_Type", line, 14))
    {
      fgetl(line,900,fp);
      goto get_transform;
    }
  }
  else
  {
    if (!strncmp("Transform_Type", line, 14))
    {
      fgetl(line, 900, fp);
      goto get_transform;
    }
    while (line[0] == '%')
      fgetl(line, 900, fp) ; /* variable # of comments */
  }

get_transform:
  m_L = lt->m_L ;
  for (row = 1 ; row <= 3 ; row++)
  {
    cp = fgetl(line, 900, fp) ;
    if (!cp)
    {
      LTAfree(&lta) ;
      ErrorReturn(NULL,
                  (ERROR_BADFILE, "ltMNIreadEx: could not read row %d from %s (%s)",
                   row, fname, line)) ;
    }
    sscanf(cp, "%f %f %f %f",
           MATRIX_RELT(m_L,row,1), MATRIX_RELT(m_L,row,2),
           MATRIX_RELT(m_L,row,3), MATRIX_RELT(m_L,row,4)) ;
  }
  if (!lta)
  {
    fclose(fp);
    return NULL;
  }
  fclose(fp);

  // add original src and dst information
  if (no_volinfo == 0)
    mincGetVolInfo(infoline, infoline2, &lta->xforms[0].src, &lta->xforms[0].dst);
  lta->type = LINEAR_RAS_TO_RAS;
  return lta;
}


#include "minc.h"

int  ltaMNIwrite(LTA *lta, char *fname)
{
  FILE             *fp ;
  int              row ;
  MATRIX           *m_L ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE, "ltMNIwrite: could not open file %s",fname));

  fprintf(fp, "MNI Transform File\n") ;
  // now saves src and dst in comment line
  fprintf(fp, "%%Generated by %s src %s dst %s\n",
          Progname, lta->xforms[0].src.fname, lta->xforms[0].dst.fname) ;
  fprintf(fp, "\n") ;
  fprintf(fp, "Transform_Type = Linear;\n") ;
  fprintf(fp, "Linear_Transform =\n") ;

  if (lta->type == LINEAR_RAS_TO_RAS)
  {
    m_L = lta->xforms[0].m_L ;
    for (row = 1 ; row <= 3 ; row++)
    {
      fprintf(fp, "      %f       %f       %f       %f",
              *MATRIX_RELT(m_L,row,1), *MATRIX_RELT(m_L,row,2),
              *MATRIX_RELT(m_L,row,3), *MATRIX_RELT(m_L,row,4)) ;
      if (row == 3)
        fprintf(fp, ";") ;
      fprintf(fp, "\n") ;
    }
  }
  else if (lta->type == LINEAR_VOX_TO_VOX)
  {
    // we use src and dst info to create RAS_TO_RAS xfm
    MATRIX *voxFromRAS = 0;
    MATRIX *rasFromVoxel = 0;
    MATRIX *tmp = 0;
    MATRIX *rasToRAS = 0;
    MRI *src = 0;
    MRI *dst = 0;
    LT  *lt = 0;
    lt = &lta->xforms[0];
    src = MRIallocHeader(lt->src.width, lt->src.height, lt->src.depth, MRI_UCHAR, 0);
    useVolGeomToMRI(&lt->src, src);
    dst = MRIallocHeader(lt->dst.width, lt->dst.height, lt->dst.depth, MRI_UCHAR, 0);
    useVolGeomToMRI(&lt->dst, dst);
    voxFromRAS = extract_r_to_i(src);
    tmp = MatrixMultiply(lta->xforms[0].m_L, voxFromRAS, NULL);
    rasFromVoxel = extract_i_to_r(dst);
    rasToRAS = MatrixMultiply(rasFromVoxel, tmp, NULL);
    for (row = 1 ; row <= 3 ; row++)
    {
      fprintf(fp, "      %f       %f       %f       %f",
              *MATRIX_RELT(rasToRAS,row,1), *MATRIX_RELT(rasToRAS,row,2),
              *MATRIX_RELT(rasToRAS,row,3), *MATRIX_RELT(rasToRAS,row,4)) ;
      if (row == 3)
        fprintf(fp, ";") ;
      fprintf(fp, "\n") ;
    }
    MatrixFree(&voxFromRAS);
    MatrixFree(&rasFromVoxel);
    MatrixFree(&tmp);
    MatrixFree(&rasToRAS);
    MRIfree(&src);
    MRIfree(&dst);
  }
  fclose(fp) ;
  return(NO_ERROR);
}
