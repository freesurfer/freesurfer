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

LTA *ltaReadFileEx(const char *fname);

char *Progname;

int main(int argc, char *argv[])
{

  char *ltafn1, *ltafn2, *ltafn_total;
  int nargs;
  LTA *lta1, *lta2, *lta_total;
  FILE *fo;
  MATRIX *r_to_i_1, *i_to_r_1, *i_to_r_2, *r_to_i_2;
  MATRIX *RAS_1_to_1, *RAS_2_to_2, *m_tmp;

  Progname = argv[0];

  nargs = handle_version_option (argc, argv, "$Id: mri_concatenate_lta.c,v 1.1 2005/02/08 16:26:18 xhan Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs ;
 
  if(argc != 4)
    usage(1);
  
  ltafn1 = argv[1];
  ltafn2 = argv[2];
  ltafn_total = argv[3];

  printf("Read individual LTAs\n");
  lta1 = ltaReadFileEx(ltafn1);
  if(!lta1)
    ErrorExit(ERROR_BADFILE, "%s: can't read file %s",Progname, ltafn1);

  lta2 = ltaReadFileEx(ltafn2);
  if(!lta2)
    ErrorExit(ERROR_BADFILE, "%s: can't read file %s",Progname, ltafn2);

  if(vg_isEqual(&lta1->xforms[0].dst, &lta2->xforms[0].src) == 0){
    /*    ErrorExit(ERROR_BADFILE, "%s: dst volume of lta1 doesn't match src volume of lta2",Progname);*/
    printf("Warning: dst volume of lta1 doesn't match src volume of lta2\n");
  }

  printf("Combine the two LTAs to get a RAS-to-RAS from src of LTA1 to dst of LTA2\n");

  if(lta1->type == LINEAR_RAS_TO_RAS){
    RAS_1_to_1 =  MatrixCopy(lta1->xforms[0].m_L, NULL);
  }else if(lta1->type == LINEAR_VOX_TO_VOX){
    r_to_i_1 = vg_r_to_i(&lta1->xforms[0].src);
    i_to_r_1 = vg_i_to_r(&lta1->xforms[0].dst);
    if(!r_to_i_1 || !i_to_r_1)
      ErrorExit(ERROR_BADFILE, "%s: failed to convert LTA1 to RAS_to_RAS",Progname);
    m_tmp = MatrixMultiply(lta1->xforms[0].m_L, r_to_i_1, NULL);
    RAS_1_to_1 = MatrixMultiply(i_to_r_1, m_tmp, NULL);
    MatrixFree(&m_tmp);
  }else{
    ErrorExit(ERROR_BADFILE, "%s: unknown transform type for LTA1",Progname);
  }

  if(lta2->type == LINEAR_RAS_TO_RAS){
    RAS_2_to_2 =  MatrixCopy(lta2->xforms[0].m_L, NULL);
  }else if(lta2->type == LINEAR_VOX_TO_VOX){
    r_to_i_2 = vg_r_to_i(&lta2->xforms[0].src);
    i_to_r_2 = vg_i_to_r(&lta2->xforms[0].dst);
    if(!r_to_i_2 || !i_to_r_2)
      ErrorExit(ERROR_BADFILE, "%s: failed to convert LTA1 to RAS_to_RAS",Progname);
    m_tmp = MatrixMultiply(lta2->xforms[0].m_L, r_to_i_2, NULL);
    RAS_2_to_2 = MatrixMultiply(i_to_r_2, m_tmp, NULL);
    MatrixFree(&m_tmp);
  }else{
    ErrorExit(ERROR_BADFILE, "%s: unknown transform type for LTA1",Progname);
  }
  
  lta_total = LTAalloc(1, NULL);
  lta_total->type = LINEAR_RAS_TO_RAS;
  MatrixMultiply(RAS_2_to_2, RAS_1_to_1, lta_total->xforms[0].m_L);
  lta_total->xforms[0].src = lta1->xforms[0].src;
  lta_total->xforms[0].dst = lta2->xforms[0].dst;
  lta_total->xforms[0].x0 = 0;   lta_total->xforms[0].y0 = 0;
  lta_total->xforms[0].z0 = 0;   lta_total->xforms[0].sigma = 1.0f; 

  printf("Write combined LTA to file %s\n", ltafn_total);
  fo = fopen(ltafn_total,"w");
  if (fo==NULL) 
    ErrorExit(ERROR_BADFILE, "%s: can't create file %s",Progname, ltafn_total);
  
  LTAprint(fo, lta_total);
  
  fclose(fo);

  LTAfree(&lta1);
  LTAfree(&lta2);
  LTAfree(&lta_total);
  MatrixFree(&RAS_1_to_1);
  MatrixFree(&RAS_2_to_2);

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
  cp = fgetl(line, 199, fp) ; sscanf(cp, "nxforms   = %d\n", &nxforms) ;
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
