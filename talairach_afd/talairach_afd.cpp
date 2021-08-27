/**
 * @brief Automatically detect Talairach alignment failures.
 *
 * Detect bad Talairach alignment by checking subjects Talairach transform
 * against a range of known good values.
 * See: dev/docs/Automatic_Failure_Detection.doc
 */
/*
 * Original Author: Laurence Wastiaux
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
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrimorph.h"
#include "mri_conform.h"
#include "utils.h"
#include "const.h"
#include "timer.h"
#include "version.h"
#include "stats.h"
#include "fio.h"
#include "transform.h"

static int get_option(int argc, char *argv[]) ;
static void usage(int exit_value) ;
static char *subject_name = NULL;
static char *xfm_fname = NULL;
static char *afd_dir = NULL;
#define DEFAULT_THRESHOLD 0.01f
static float threshold = DEFAULT_THRESHOLD;
static int verbose=0;
VECTOR *ReadMeanVect(char *fname);
MATRIX *ReadCovMat(char *fname);
VECTOR *Load_xfm(char *fname);
double *LoadProbas(char *fname, int *nb_samples);
float mvnpdf(VECTOR *t, VECTOR *v_mu, MATRIX *m_sigma);
float ComputeArea(HISTO *h, int nbin);

int main(int argc, char *argv[]) ;

const char *Progname ;

int main(int argc, char *argv[])
{
  char **av;
  char *fsenv=NULL, *sname=NULL;
  char fsafd[1000], tmf[1000], cvf[1000], xfm[1000], probasf[1000];
  int ac, nargs;
  int b;
  int nsamples, bin;
  float p, pval;
  double *ts_probas;
  VECTOR *mu;
  MATRIX *sigma;
  VECTOR *txfm;
  HISTO *h;
  int ret_code=0; // assume 'passed' return code

  mu=NULL;
  sigma=NULL;
  txfm=NULL;

  nargs = handleVersionOption(argc, argv, "talairach_afd");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 1)
  {
    usage(1);
  }

  if ((subject_name==NULL) && (xfm_fname==NULL))
  {
    printf("ERROR: a subject name or .xfm filename must be provided!\n");
    usage(1);
  }

  if (subject_name != NULL)
  {
    sname=fio_basename(subject_name, NULL);
  }
  else if (xfm_fname)
  {
    sname=xfm_fname;
  }

  if (xfm_fname == NULL)
  {
    sprintf(xfm, "%s/mri/transforms/talairach.xfm", subject_name);
  }
  else
  {
    // use command-line specified transform
    sprintf(xfm,"%s",xfm_fname);
  }

  if (afd_dir == NULL)
  {
    fsenv=getenv("FREESURFER_HOME");
    if(!fsenv)
    {
      printf("ERROR: FREESURFER_HOME not defined!\n");
      printf("INFO: optionally use the -afd <afd directory> flag\n");
      exit(1);
    }
    sprintf(fsafd, "%s/fsafd", fsenv);
  }
  else
  {
    // path to .afd files was specified on the command line
    strcpy(fsafd, afd_dir);
  }

  int req = snprintf(tmf, 1000, "%s/TalairachingMean.adf", fsafd);
  if (req >= 1000) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  req = snprintf(cvf, 1000, "%s/TalairachingCovariance.adf", fsafd);
  if (req >= 1000) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  req = snprintf(probasf, 1000, "%s/TalairachingProbas.adf", fsafd);
  if (req >= 1000) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  /*----- load 1x9 mean vector computed from training set -----*/
  mu = ReadMeanVect(tmf);

  /*----- load 9x9 cov matrix computed from training set -----*/
  sigma = ReadCovMat(cvf);

  /*--- Load training set's probas ---*/
  ts_probas = LoadProbas(probasf, &nsamples);

  h = HISTObins(21, 0,1);
  HISTOcount(h, ts_probas, nsamples);

  for (b = 0 ; b < h->nbins ; b++)
  {
    h->counts[b] = (float)h->counts[b]/(float)nsamples ;
  }

  /*---Load Talairach.xfm---*/
  txfm = Load_xfm(xfm);

  /*--- compute mvnpdf ---*/
  p = mvnpdf(txfm, mu, sigma);

  bin = HISTOvalToBin(h,p);

  pval = ComputeArea(h, bin);

  if(pval < threshold)
  {
    ret_code=1; // return a failure code when exiting
    printf("ERROR: talairach_afd: Talairach Transform: %s ***FAILED***"
           " (p=%.4f, pval=%.4f < threshold=%.4f)\n",
           sname, p, pval, threshold);
  }
  else
  {
    printf("talairach_afd: Talairach Transform: %s OK "
           "(p=%.4f, pval=%.4f >= threshold=%.4f)\n",
           sname, p, pval, threshold);
  }

  if (verbose)
  {
    int i,j;
    if (fsenv)
    {
      printf("\nFREESURFER_HOME = %s\n", fsenv);
    }
    printf("\nfsafdDir = %s\n", fsafd);

    printf("\nmu:\n");
    for (i=1; i<=9; i++)
    {
      float m= mu->rptr[1][i];
      if (m<0)
      {
        printf("%.4f  ", m);
      }
      else
      {
        printf(" %.4f  ", m);
      }
    }
    printf("\n\n");
    printf("sigma:\n");
    for(i=1; i<=9; i++)
    {
      for(j=1; j<=9; j++)
      {
        float s=sigma->rptr[i][j];
        if (s<0)
        {
          printf("%.4f  ", s);
        }
        else
        {
          printf(" %.4f  ", s);
        }
      }
      printf("\n");
    }
    printf("\n");
    printf("xfm:\n");
    for(i=1; i<=9; i++)
    {
      float t=txfm->rptr[1][i];
      if (t<0)
      {
        printf("%.4f  ", t);
      }
      else
      {
        printf(" %.4f  ", t);
      }
    }
    printf("\n\n");

    printf("proba = %.4f\n\n", p);
  }

  exit(ret_code) ;
  return(ret_code) ;
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void print_version(void)
{
  fprintf(stdout, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    usage(1) ;
  }
  else if (!stricmp(option, "version"))
  {
    print_version() ;
  }
  if (!stricmp(option, "subj"))
  {
    subject_name = argv[2];
    nargs = 1;
  }
  else if (!stricmp(option, "xfm"))
  {
    xfm_fname = argv[2];
    nargs = 1;
  }
  else if (!stricmp(option, "afd"))
  {
    afd_dir = argv[2];
    nargs = 1;
  }
  else if (!stricmp(option, "T") || !stricmp(option, "threshold"))
  {
    threshold = (float)atof(argv[2]);
    nargs = 1;
  }
  else switch (toupper(*option))
    {
    case 'T':
    case '?':
    case 'H':
    case 'U':
      usage(0) ;
      break ;
    case 'V':
      verbose=1;
      break;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
#include "talairach_afd.help.xml.h"
static void
usage(int exit_value)
{
  outputHelpXml(talairach_afd_help_xml,talairach_afd_help_xml_len);
  exit(exit_value) ;
}


/*----------------------------------------------------------------------
  Parameters:char *fname

  Description:reads the 9 componants of Mean Vector mu
  from "$FREESURFER_HOME/fsafd/TalairachingMean.adf"
  ----------------------------------------------------------------------*/
VECTOR *
ReadMeanVect(char *fname)
{
  FILE *fp;
  char tag[1000], tmpstr[1000];
  int nth, r;
  VECTOR *v;

  v = MatrixAlloc(1,9, MATRIX_REAL) ;

  fp = fopen(fname,"r");
  if(fp == NULL)
  {
    printf("ERROR: talairach_afd::ReadMeanVect(): could not open %s\n",
           fname);
    exit(1);
  }

  nth=1;
  while(1)
  {
    r = fscanf(fp,"%s",tag);
    if(r==EOF)
    {
      break;
    }
    if(!strcmp(tag,"#"))   // line starts with #:
    {
      // header info -> skip the whole line
      fgets(tmpstr,1000,fp);
    }
    else
    {
      v->rptr[1][nth]=(float)atof(tag);
      nth++;
    }
  }
  if (nth!=10)
  {
    printf("ERROR: talairach_afd::ReadMeanVect(): "
           "9 components could not be loaded\n");
    exit(1);
  }

  fclose(fp);
  return(v);
}


/*----------------------------------------------------------------------
  Parameters:char *fname

  Description:reads the 9x9 componants of covariance matrix sigma[9][9]
  from "$FREESURFER_HOME/fsafd/TalairachingCovariance.adf"
  ----------------------------------------------------------------------*/
MATRIX *
ReadCovMat(char *fname)
{
  FILE *fp;
  char tag[1000], tmpstr[1000];
  int nth, r, i, j;
  MATRIX *s;

  s = MatrixAlloc(9,9, MATRIX_REAL) ;

  fp = fopen(fname,"r");
  if(fp == NULL)
  {
    printf("ERROR: talairach_afd::ReadCovMAt(): could not open %s\n",fname);
    exit(1);
  }

  nth=0;
  while(1)
  {
    r = fscanf(fp,"%s",tag);
    if(r==EOF)
    {
      break;
    }
    if(!strcmp(tag,"#"))   // line starts with #:
    {
      // header info -> skip the whole line
      fgets(tmpstr,1000,fp);
    }
    else
    {
      i=nth/9;
      j=nth-(i*9);
      s->rptr[i+1][j+1]=(float)atof(tag);
      nth++;
    }
  }
  if (nth!=81)
  {
    printf("ERROR: talairach_afd::ReadCovMAt(): "
           "9x9 components could not be loaded\n");
    exit(1);
  }

  fclose(fp);
  return(s);
}


/*----------------------------------------------------------------------
  Parameters:char *fname

  Description:Load xfm matrix
  ----------------------------------------------------------------------*/
VECTOR *
Load_xfm(char *fname)
{
  FILE *fp;
  char *sdir;
  char tmp_name[1000];
  VECTOR *v;
  int i,j;
  LTA *lta = LTAreadEx(fname);

  v = MatrixAlloc(1,9, MATRIX_REAL) ;

  fp = fopen(fname,"r");
  if(!fp)
  {
    // try to find subject in SUBJECTS_DIR
    sdir=getenv("SUBJECTS_DIR");
    if(sdir)
    {
      sprintf(tmp_name,"%s/%s", sdir, fname);
      fp = fopen(tmp_name,"r");
      if(!fp)
      {
        printf("ERROR: talairach_afd::Load_xfm(): could not open %s\n",fname);
        exit(1);
      }
      else
      {
        fname = tmp_name;
      }
    }
    else
    {
      printf("ERROR: talairach_afd::Load_xfm(): could not open %s\n",fname);
      if (subject_name)
      {
        printf("SUBJECTS_DIR is not defined!\n");
      }
      exit(1);
    }
  }
  fclose(fp);
  // at this point, 'fname' will successfully open

  // parse .xfm file
  if (lta==NULL)
  {
    printf("ERROR: talairach_afd::Load_xfm(): could not parse %s\n",fname);
    exit(1);
  }
  //MatrixPrint(stdout,lta->xforms->m_L);
  for(i=0; i<3; i++)
  {
    for(j=0; j<3; j++)
    {
      v->rptr[1][i*3+j+1]=lta->xforms->m_L->rptr[i+1][j+1];
    }
  }
  //MatrixPrint(stdout,v);

  return(v);
}


/*----------------------------------------------------------------------
  Parameters:char *fname

  Description:Load probas computed from the training data set
  and contained in probas table
  ----------------------------------------------------------------------*/
double *
LoadProbas(char *fname, int *nb_samples)
{
  FILE *fp;
  char tag[1000], tmpstr[1000], c[1000];
  int  nth, r,  len = 0;
  double *p=NULL;

  fp = fopen(fname,"r");
  if(fp == NULL)
  {
    printf("ERROR: talairach_afd::LoadProbas(): could not open %s\n",fname);
    exit(1);
  }

  nth=0;

  while(1)
  {
    r = fscanf(fp,"%s",tag);
    if(r==EOF)
    {
      break;
    }
    if(!strcmp(tag,"#"))   // line starts with #:
    {
      // header info -> skip the whole line
      fgets(tmpstr,1000,fp);
      if(strcmp(tmpstr, " nrows")>0)
      {
        sscanf(tmpstr, "%s %d", c, &len);
        break;
      }
    }
  }

  p=(double *)calloc(len, sizeof(double));
  if(!p)
  {
    printf("ERROR: talairach_afd::LoadProbas(): "
           "could not allocate memory for probas vector\n");
    exit(1);
  }

  fp=freopen(fname,"r", fp);
  if(fp == NULL)
  {
    printf("ERROR: talairach_afd::LoadProbas(): could not re-open %s\n",
           fname);
    exit(1);
  }

  while(1)
  {
    r = fscanf(fp,"%s",tag);
    if(r==EOF)
    {
      break;
    }
    if(!strcmp(tag,"#"))   // line starts with #:
    {
      // header info -> skip the whole line
      fgets(tmpstr,1000,fp);
    }
    else
    {
      p[nth] = (double)atof(tag);
      nth++;
    }
  }

  if(nth!=len)
  {
    printf("ERROR: talairach_afd::LoadProbas(): "
           "%d components could not be loaded\n", len);
    exit(1);
  }
  *(nb_samples)=len;

  fclose(fp);
  return(p);

}

/*----------------------------------------------------------------------
  Parameters:vector t, mean vector mu, covariance matrix sigma

  Description:computes p = mvnpdf(t, mu, sigma)
  ----------------------------------------------------------------------*/

float
mvnpdf(VECTOR *t, VECTOR *v_mu, MATRIX *m_sigma)
{
  VECTOR *sub, *t_sub;
  MATRIX *m_sinv, *m_tmp;
  float  num, det, den, tmp;
  int d;

  if(!t || !v_mu || !m_sigma)
  {
    fprintf(stderr, "ERROR: talairach_afd::mvnpdf(): cannot compute mvnpdf\n");
    exit(1);
  }

  d = t->cols;

  sub = VectorSubtract(t,v_mu,NULL);
  m_sinv = MatrixInverse(m_sigma, NULL);
  t_sub = VectorTranspose(sub, NULL);

  m_tmp = MatrixMultiply(m_sinv, t_sub, NULL);
  m_tmp = MatrixMultiply(sub, m_tmp, NULL);
  tmp = m_tmp->rptr[1][1];
  num = exp(- (tmp/2));
  det = MatrixDeterminant(m_sigma);
  den = 1/(sqrt(det));

  num = 1/sqrt(pow(2*M_PI,d)) * num * den;

  return(num);
}


/*----------------------------------------------------------------------
  Parameters:

  Description:computes area under the curve for b<=nbin
  ----------------------------------------------------------------------*/
float
ComputeArea(HISTO *h, int nbin)
{
  float a, s;
  int i;

  s=0;
  for(i=1 ; i<=nbin ; i++)
  {
    a = (h->counts[i] + h->counts[i-1]) /2;
    s += a;
  }

  return(s);
}

