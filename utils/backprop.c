/**
* @file  backprop.c
* @brief back-propagation neural net
*
* backprop files can store multiple networks of varying size. In order to
* accomodate this, I mimic the tiff file structure with a header indicating
* the # of networks, and a pointer to the 1st one. Each one then starts with
* a pointer to the next network, or NULL for the last network in the file.
*/
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:42 $
 *    $Revision: 1.13 $
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

/*-----------------------------------------------------------------
              INCLUDE FILES
-----------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <memory.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include "machine.h"
#include "backprop.h"
#include "diag.h"
#include "utils.h"
#include "macros.h"
#include "proto.h"
#include "error.h"

/*-----------------------------------------------------------------
              MACROS AND CONSTANTS
-----------------------------------------------------------------*/

#define InsFree          free
#define InsCalloc        calloc
#define InsHalloc(a,b)   calloc((int)a, (int)b)
#define InsHfree         free
#define far
#define huge

/*
  can use either sigmoid or hyperbolic tangent activation function.
*/
#define SIGMOID 0

#if SIGMOID

#define f(net)      (1.0f / (1.0f + (float)exp((double)-net)))
#define fprime(o)   (o * (1.0f - o))

#define SIG_MAX      0.75f
#define SIG_MIN      0.25f

#else

/* use hyperbolic tangent */

#define T_A  1.0f /*1.716f*/   /* push limits of sigmoid away from +-1 */
#define T_B  0.666f   /* stretch sigmoid out a bit */

#define f(net)      (T_A * tanh(T_B * (net)))
#define fprime(o)   (T_A * T_B * (1.0f - o*o))

#define SIG_MAX      0.5f
#define SIG_MIN      -0.5f

#endif

#define SIG_RANGE    (SIG_MAX - SIG_MIN)
#define MIN_TRATE    0.01f

/*----------------------------------------------------------------------
                STRUCTURES
----------------------------------------------------------------------*/

/*-----------------------------------------------------------------
                PROTOTYPES
-----------------------------------------------------------------*/

static void bpInitLayerWeights(LAYER *layer) ;
static void bpInitLayer(LAYER *backprop, int ninputs, int nunits) ;
static void bpFreeLayer(LAYER *layer) ;
static void bpWriteLayer(FILE *fp, LAYER *layer) ;
static void bpReadLayer(FILE *fp, LAYER *layer) ;
static void bpLayerFeedForward(float *I, LAYER *layer, int nlin) ;
static void bpCalculateOutputDeltas(BACKPROP *backprop, float *targets) ;
static void bpCalculateHiddenDeltas(BACKPROP *backprop) ;
static void bpUpdateLayerWeights(LAYER *layer, float *I,float trate,
                                 float momentum);
static void bpUnnormalizeOutputs(BACKPROP *backprop) ;
static void bpNormalizeTargets(BACKPROP *backprop, float *targets) ;
static void bpUnnormalizeTargets(BACKPROP *backprop, float *targets) ;
static void bpCopyLayer(LAYER *lsrc, LAYER *ldst) ;

/* file maintainance stuff */
static long bpFileNewEnd(FILE *fp, BPFILE_HEADER *hd, int swapped) ;
static long bpFileSeekEndPtr(FILE *fp, BPFILE_HEADER *hd, int swapped) ;
static long bpFileSeekNet(FILE *fp, BPFILE_HEADER *hd, int netno, int swapped);
static int  bpChangeNumberOfNets(FILE *fp, int nnets) ;

/*-----------------------------------------------------------------
                   MACROS
-----------------------------------------------------------------*/

/* from i to j */
#define Wij(l, i, j)  (*((l)->w + (((long)j) * (l)->ninputs) + (long)i))
#define DWij(l, i, j)  (*((l)->dw + (((long)j) * (l)->ninputs) + (long)i))

/*-----------------------------------------------------------------
                STATIC DATA
-----------------------------------------------------------------*/

/*-----------------------------------------------------------------
                GLOBAL DATA
-----------------------------------------------------------------*/

/*-----------------------------------------------------------------
                FUNCTIONS
-----------------------------------------------------------------*/

/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int coef = -1 ;

BACKPROP *
BackpropAlloc(int ninputs, int noutputs, int nhidden, float trate,
              float momentum, float *mean_out, float *std_out)
{
  BACKPROP  *backprop ;

  backprop = (BACKPROP *)InsCalloc(1, sizeof(BACKPROP)) ;

  backprop->ninputs = ninputs ;
  backprop->noutputs = noutputs ;
  backprop->nhidden = nhidden ;
  backprop->old_momentum = backprop->momentum = momentum ;
  backprop->mean_out = (float *)calloc(noutputs, sizeof(float)) ;
  backprop->error_ratio = BP_ERROR_RATIO ;
  backprop->trate_up = BP_TRATE_INCREASE ;
  backprop->trate_down = BP_TRATE_DECREASE ;
  if (!backprop->mean_out)
    ErrorExit(ERROR_BAD_FILE,  "BackpropAlloc: could not output range vector\n") ;
  backprop->std_out = (float *)calloc(noutputs, sizeof(float)) ;
  if (!backprop->std_out)
    ErrorExit(ERROR_BAD_FILE,  "BackpropAlloc: could not output range vector\n") ;
  memmove(backprop->mean_out, mean_out, noutputs*sizeof(float)) ;
  memmove(backprop->std_out, std_out, noutputs*sizeof(float)) ;
  backprop->errors = (float *)calloc(noutputs, sizeof(float)) ;
  if (!backprop->errors)
    ErrorExit(ERROR_BAD_FILE,  "BackpropAlloc: could not allocate error vector\n") ;

  backprop->trate = trate ;
  bpInitLayer(&backprop->hidden, ninputs, nhidden) ;
  bpInitLayer(&backprop->output, nhidden, noutputs) ;

  return(backprop) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
BackpropFileNumberOfNets(char *fname)
{
  BPFILE_HEADER  hd ;
  FILE           *fp ;

  if (strcmp(fname, "-"))
  {
    fp = fopen(fname, "r") ;
    if (!fp)
      return(-1) ;
    fclose(fp) ;
    fp = fopen(fname, "r+") ;
  }
  else
    fp = stdin ;

  if (fread(&hd, sizeof(hd), 1, fp) != 1)
    ErrorExit(ERROR_BAD_FILE,
              "BackpropNumberOfNets(%s): could not read header from file\n",
              fname) ;

  DiagPrintf(DIAG_WRITE,"BackpropNumberOfNets(%s): %d nets, 1st at %ld\n",
             fname, hd.nnets, hd.first) ;
  if (hd.magic != BP_MAGIC)  /* try changing byte order */
  {
    hd.magic = swapLong32(hd.magic) ;
    hd.nnets = swapLong32(hd.nnets) ;
    hd.first = swapLong32(hd.first) ;
  }

  if (hd.magic != BP_MAGIC)
    ErrorReturn(ERROR_BAD_FILE, (ERROR_BAD_FILE, "BackpropFileNumberOfNets(%s): bad magic #\n",
                                 fname)) ;

  return((int)hd.nnets) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
BACKPROP *
BackpropRead(char *fname, int netno)
{
  BACKPROP       *backprop ;
  int            ninputs, noutputs, nhidden, ubytes, swapped = 0 ;
  FILE           *fp ;
  float          momentum, trate, *min_out, *max_out ;
  BPFILE_HEADER  hd ;

  if (strcmp(fname, "-"))
  {
    fp = fopen(fname, "r") ;
    if (!fp)
      return(NULL) ;
    fclose(fp) ;
    fp = fopen(fname, "r+") ;
  }
  else
    fp = stdin ;

  if (fread(&hd, sizeof(hd), 1, fp) != 1)
    ErrorReturn(NULL, (ERROR_BAD_FILE,
                       "BackpropRead(%s): could not read header from file\n", fname)) ;

  if (hd.magic != BP_MAGIC)  /* try changing byte order */
  {
    hd.magic = swapLong32(hd.magic) ;
    hd.nnets = swapLong32(hd.nnets) ;
    hd.first = swapLong32(hd.first) ;
    swapped = 1 ;
  }
  if (hd.magic != BP_MAGIC)
    ErrorReturn(NULL, (ERROR_BAD_FILE, "BackpropRead(%s): bad magic #\n", fname)) ;

  if (netno < 0)        /* read last net in file */
    netno = (int)hd.nnets ;

  DiagPrintf(DIAG_WRITE, "BackpropRead(%s, %d): %d nets, 1st at %ld\n",
             fname, netno, hd.nnets, hd.first) ;
  if (netno >= hd.nnets)
  {
    fclose(fp) ;
    return(NULL) ;
  }

  bpFileSeekNet(fp, &hd, netno, swapped) ;
  DiagPrintf(DIAG_WRITE, "after bpFileSeekNet: fpos @ %ld\n", ftell(fp)) ;

  if (fscanf(fp, "%d %d %d %f %f %d\n", &ninputs, &noutputs,
             &nhidden, &trate, &momentum, &ubytes) != 6)
  {
    if (fp != stdin)
      fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BAD_FILE, "could not scan parameters from bp file %s\n", fname)) ;
  }
  min_out = (float *)calloc(noutputs, sizeof(float)) ;
  if (!min_out)
  {
    if (fp != stdin) fclose(fp) ;
    ErrorExit(ERROR_BAD_FILE,  "BackpropRead(%s): could not allocate min_out\n",
              fname) ;
  }
  max_out = (float *)calloc(noutputs, sizeof(float)) ;
  if (!max_out)
  {
    if (fp != stdin)
      fclose(fp) ;
    ErrorExit(ERROR_BAD_FILE,  "BackpropRead(%s): could not allocate max_out\n",
              fname) ;
  }
  if (fread(min_out, sizeof(float), noutputs, fp) != (size_t)noutputs)
  {
    if (fp != stdin)
      fclose(fp) ;
    ErrorExit(ERROR_BAD_FILE,  "BackpropRead(%s): could not read min_out\n",
              fname) ;
  }
  if (fread(max_out, sizeof(float), noutputs, fp) != (size_t)noutputs)
  {
    if (fp != stdin)
      fclose(fp) ;
    ErrorExit(ERROR_BAD_FILE,  "BackpropRead(%s): could not read max_out\n",
              fname) ;
  }
  if (swapped)
  {
    *min_out = swapFloat(*min_out) ;
    *max_out = swapFloat(*max_out) ;
  }
  backprop =
    BackpropAlloc(ninputs, noutputs, nhidden, momentum, trate,min_out,max_out);

  free(min_out) ;
  free(max_out) ;
  if (ubytes > 0)
  {
    backprop->user = (char *)calloc(ubytes, sizeof(char)) ;
    backprop->user_bytes = ubytes ;
    if (fread(backprop->user, sizeof(char), ubytes, fp) != (size_t)ubytes)
    {
      if (fp != stdin)
        fclose(fp) ;
      ErrorReturn(NULL,(ERROR_BAD_FILE, "BackpropRead: could not read %d user bytes\n",ubytes));
    }
  }
  bpReadLayer(fp, &backprop->hidden) ;
  bpReadLayer(fp, &backprop->output) ;

  if (fp != stdin)
    fclose(fp) ;

  return(backprop) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
BACKPROP *
BackpropCopy(BACKPROP *bp_src, BACKPROP *bp_dst)
{
  if (bp_dst &&
      ((bp_dst->ninputs != bp_src->ninputs) ||
       (bp_dst->noutputs != bp_src->noutputs) ||
       (bp_dst->nhidden != bp_src->nhidden)))
    BackpropFree(&bp_dst) ;

  if (!bp_dst)
    bp_dst = BackpropAlloc(bp_src->ninputs, bp_src->noutputs,
                           bp_src->nhidden, bp_src->trate,
                           bp_src->momentum, bp_src->mean_out,
                           bp_src->std_out) ;

  bp_dst->momentum = bp_src->momentum ;      /* momentum coefficient */
  bp_dst->trate = bp_src->trate ;         /* learning rate */
  bp_dst->trate_up = bp_src->trate_up ;      /* adaptive step size increase */
  bp_dst->trate_down = bp_src->trate_down ; /* adaptive step size decrease */

  /* ratio of new/old sse for step size modification */
  bp_dst->error_ratio = bp_src->error_ratio ;
  bp_dst->nepochs = bp_src->nepochs ;    /* # of epochs of training */
  bp_dst->ntrials = bp_src->ntrials ;
  bp_dst->old_momentum = bp_src->old_momentum ;
  bp_dst->sse = bp_src->sse ;

  memmove(bp_dst->mean_out, bp_src->mean_out, bp_src->noutputs * sizeof(float));
  memmove(bp_dst->std_out, bp_src->std_out, bp_src->noutputs * sizeof(float)) ;
  memmove(bp_dst->errors, bp_src->errors, bp_src->noutputs * sizeof(float)) ;
  if ((bp_dst->user_bytes = bp_src->user_bytes) != 0)
  {
    bp_dst->user = (void *)calloc(bp_src->user_bytes, sizeof(char)) ;
    memmove(bp_dst->user, bp_src->user, bp_src->user_bytes) ;
  }

  bpCopyLayer(&bp_src->hidden, &bp_dst->hidden) ;
  bpCopyLayer(&bp_src->output, &bp_dst->output) ;
  return(bp_dst) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpCopyLayer(LAYER *lsrc, LAYER *ldst)
{
  int nweights ;

  nweights = lsrc->nunits * lsrc->ninputs ;
  memmove(ldst->w, lsrc->w, nweights*sizeof(*(lsrc->w))) ;
  memmove(ldst->dw, lsrc->dw, nweights*sizeof(*(lsrc->dw))) ;
  memmove(ldst->biases, lsrc->biases, lsrc->nunits*sizeof(*(lsrc->biases))) ;
  memmove(ldst->db, lsrc->db, lsrc->nunits*sizeof(*(lsrc->db))) ;
  if (lsrc->deltas)
  {
    if (!ldst->deltas)
      ldst->deltas = (float *)calloc(lsrc->nunits, sizeof(float)) ;
    memmove(ldst->deltas, lsrc->deltas, lsrc->nunits*sizeof(*(lsrc->deltas))) ;
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
BackpropWrite(BACKPROP *backprop, char *fname, int argc, char *argv[],
              char *comments, int mode)
{
  FILE           *fp = NULL ;
  int            i ;
  char           *user, *time_str ;
  time_t         tt ;
  BPFILE_HEADER  hd ;
  long           fpos ;

  DiagPrintf(DIAG_WRITE, "BackpropWrite(%s): mode = %s\n",
             fname, mode == BP_WRITE ? "write" : "append") ;

  switch (mode)
  {
  case BP_WRITE:   /* create a new file */
    fp = fopen(fname, "w") ;
    if (!fp)
      ErrorExit(ERROR_BAD_FILE,  "BackpropWrite(%s): could not create file\n",fname) ;
    hd.nnets = 1 ;
    hd.magic = BP_MAGIC ;
    hd.first = 0L ;   /* will be filled in with ptr to 1st network */
    if (fwrite(&hd, sizeof(hd), 1, fp) != 1)
      ErrorExit(ERROR_BAD_FILE,  "BackpropWrite(%s): could not write header\n", fname) ;
    fclose(fp) ;
    fp = fopen(fname, "r+") ;
    if (!fp)
      ErrorExit(ERROR_BAD_FILE, "BackpropWrite(%s): could not create file for appending\n",
                fname) ;
    break ;
  case BP_APPEND:
    fp = fopen(fname, "r+") ;
    if (!fp)
      ErrorReturn(-3, (ERROR_BAD_FILE,
                       "BackpropWrite(%s): could not create file for appending\n", fname)) ;

    if (fread(&hd, sizeof(hd), 1, fp) != 1)
      ErrorReturn(-3, (ERROR_BAD_FILE,
                       "BackpropWrite(%s): could not read header from file\n", fname)) ;

    if (hd.magic != BP_MAGIC)
      ErrorReturn(-4, (ERROR_BAD_FILE, "BackpropWrite(%s): bad magic #\n", fname)) ;
    break ;
  }

  fpos = bpFileNewEnd(fp, &hd, 0) ;/* seek to the last net in the file */
  DiagPrintf(DIAG_WRITE, "after bpFileNewEnd, fpos @ %ld\n", fpos) ;

  /* put NULL pointer to indicate the end of the ptr chain */
  if (fseek(fp, 0L, SEEK_END))
    ErrorExit(ERROR_BAD_FILE, "BackpropWrite(%s): could not seek to end of file\n",fname);

  fpos = 0L ;
  if (fwrite(&fpos, sizeof(fpos), 1, fp) != 1)
    ErrorExit(ERROR_BAD_FILE,  "BackpropWrite(%s): could not write null\n", fname) ;

  DiagPrintf(DIAG_WRITE, "nulls written in file, fpos @ %ld\n", ftell(fp)) ;
  fprintf(fp, "%d %d %d %f %f %d\n",
          backprop->ninputs, backprop->noutputs, backprop->nhidden,
          backprop->trate, backprop->momentum, backprop->user_bytes) ;

  if (fwrite(backprop->mean_out, sizeof(float), backprop->noutputs, fp) !=
      (size_t)backprop->noutputs)
    ErrorExit(ERROR_BAD_FILE,  "BackpropWrite(%s): could not write min output vector\n",
              fname) ;

  if (fwrite(backprop->std_out, sizeof(float), backprop->noutputs, fp) !=
      (size_t)backprop->noutputs)
    ErrorExit(ERROR_BAD_FILE,  "BackpropWrite(%s): could not write max output vector\n",
              fname) ;

  if (backprop->user_bytes > 0)
  {
    if (fwrite(backprop->user, sizeof(char), backprop->user_bytes, fp) !=
        (size_t)backprop->user_bytes)
      ErrorReturn(ERROR_BAD_FILE,
                  (ERROR_BAD_FILE, "BackpropWrite: could not write %d user bytes\n",
                   backprop->user_bytes)) ;
  }

  bpWriteLayer(fp, &backprop->hidden) ;
  bpWriteLayer(fp, &backprop->output) ;

  user = getenv("USER") ;
  if (!user)
    user = getenv("LOGNAME") ;
  if (!user)
    user = "UNKNOWN" ;

  tt = time(&tt) ;
  time_str = ctime(&tt) ;
  fprintf(fp, "\ncreated by %s on %s\n", user, time_str) ;
  for (i = 0 ; i < argc ; i++)
    fprintf(fp, "%s ", argv[i]) ;
  fprintf(fp, "\n") ;
  if (comments)
    fprintf(fp, "%s\n", comments) ;

  if (mode == BP_APPEND)
    bpChangeNumberOfNets(fp, 1) ;

  if (fp != stdout)
    fclose(fp) ;
  return(0) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
BackpropProcess(BACKPROP *backprop, float *I)
{
  int     i, class = -1 ;
  LAYER   *hidden, *output ;
  float   maxX ;

  if (Gdiag & DIAG_BACKPROP)
  {
    if (!backprop->learn)
      printf("\n") ;
    printf("BackpropProcess(") ;
    for (i = 0 ; i < backprop->ninputs ; i++)
      printf("%2.3f ", I[i]) ;
    printf(")\n") ;
  }

  hidden = &backprop->hidden ;
  output = &backprop->output ;
  bpLayerFeedForward(I, hidden, 1) ;
  bpLayerFeedForward(hidden->x, output, 0) ;

  for (maxX = 0.0f, i = 0 ; i < backprop->noutputs ; i++)
  {
    if (output->x[i] > maxX)
    {
      maxX = output->x[i] ;
      class = i ;
    }
  }
  /* scale outputs to desired range */
  if (!backprop->learn)
    bpUnnormalizeOutputs(backprop) ;

  if (Gdiag & DIAG_BACKPROP)
  {
    printf("BackpropProcess returning: ") ;
    for (i = 0 ; i  < backprop->noutputs ; i++)
      printf("%2.5f ", backprop->output.x[i]) ;
    printf("\n") ;
  }
  return(class) ;
}
#if 1
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
       convert the output range to unnormalized format
----------------------------------------------------------------------*/
#define NO_CONVERT 0
static void
bpUnnormalizeOutputs(BACKPROP *backprop)
{
  float   mean, std ;
  int     i ;
  LAYER   *output ;

  /* scale outputs to desired range */
  output = &backprop->output ;
  for (i = 0 ; i < backprop->noutputs ; i++)
  {
    mean = backprop->mean_out[i] ;
    std = backprop->std_out[i] ;
    output->x[i] = output->x[i] * std + mean ;
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
      scale targets to desired range
----------------------------------------------------------------------*/
static void
bpNormalizeTargets(BACKPROP *backprop, float *targets)
{
  float   mean, std ;
  int     i ;

  /* scale targets to desired range */
  for (i = 0 ; i < backprop->noutputs ; i++)
  {
    mean = backprop->mean_out[i] ;
    std = backprop->std_out[i] ;
    targets[i] = (targets[i] - mean) / std ;
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
       convert targets back to original unnormalized range
----------------------------------------------------------------------*/
static void
bpUnnormalizeTargets(BACKPROP *backprop, float *targets)
{
  float   mean, std ;
  int     i ;

  /* scale targets to desired range */
  for (i = 0 ; i < backprop->noutputs ; i++)
  {
    mean = backprop->mean_out[i] ;
    std = backprop->std_out[i] ;
    targets[i] = (targets[i] * std) + mean ;
  }
}
#else
/*----------------------------------------------------------------------
Parameters:

Description:

Returns:
scale outputs to desired range
----------------------------------------------------------------------*/
#define NO_CONVERT 0
static void
bpUnnormalizeOutputs(BACKPROP *backprop)
{
  float   scale, min_out ;
  int     i ;
  LAYER   *output ;

#if NO_CONVERT
  return ;
#endif

  /* scale outputs to desired range */
  output = &backprop->output ;
  for (i = 0 ; i < backprop->noutputs ; i++)
  {
    min_out = backprop->min_out[i] ;
    scale = (backprop->max_out[i] - min_out) / SIG_RANGE ;
    output->x[i] = (output->x[i] - SIG_MIN) * scale + min_out ;
  }
}
/*----------------------------------------------------------------------
Parameters:

Description:

Returns:
scale outputs to desired range
----------------------------------------------------------------------*/
static void
bpNormalizeTargets(BACKPROP *backprop, float *targets)
{
  float   scale, min_out ;
  int     i ;

#if NO_CONVERT
  return ;
#endif

  /* scale outputs to desired range */
  for (i = 0 ; i < backprop->noutputs ; i++)
  {
    min_out = backprop->min_out[i] ;
    scale = (backprop->max_out[i] - min_out) / SIG_RANGE ;
    targets[i] = (targets[i] - min_out) / scale + SIG_MIN ;
  }
}
/*----------------------------------------------------------------------
Parameters:

Description:

Returns:
scale outputs to desired range
----------------------------------------------------------------------*/
static void
bpUnnormalizeTargets(BACKPROP *backprop, float *targets)
{
  float   scale, min_out ;
  int     i ;

#if NO_CONVERT
  return ;
#endif

  /* scale outputs to desired range */
  for (i = 0 ; i < backprop->noutputs ; i++)
  {
    min_out = backprop->min_out[i] ;
    scale = (backprop->max_out[i] - min_out) / SIG_RANGE ;
    targets[i] = (targets[i] - SIG_MIN) * scale + min_out ;
  }
}
#endif
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpLayerFeedForward(float *I, LAYER *layer, int nlin)
{
  int             i, j, nunits, ninputs ;
  float           *pbias, *px ;
  register float  net, *Ii, *wij ;

  if (Gdiag & DIAG_BACKPROP)
    printf("bpLayerFeedForward()\n") ;

  nunits = layer->nunits ;
  ninputs = layer->ninputs ;
  pbias = &layer->biases[0] ;
  px = &layer->x[0] ;
  for (j = 0 ; j < nunits ; j++)
  {
    net = *pbias++ ;
    if (Gdiag & DIAG_BACKPROP)
      printf("unit %d: %2.5f\n", j, net) ;

    Ii  = &I[0] ;
    wij = &Wij(layer, 0, j) ;
    for (i = 0 ; i < ninputs ; i++)
    {
      net += *wij++ * *Ii++ ;
#if 0
      if (Gdiag & DIAG_BACKPROP)
        printf("net += %2.5f * %2.3f --> %2.5f\n",
               Wij(layer, i, j), I[i], net) ;
#endif
    }

    if (nlin)
      *px++ = (float)f(net) ;
    else
      *px++ = net ;
#if 0
    if (Gdiag & DIAG_BACKPROP)
      printf("unit %d output f(%2.5f) = %2.5f\n", j, net, f(net)) ;
#endif
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:
        train the network for one entire epoch, randomizing the training
        order.

    Returns:
----------------------------------------------------------------------*/
extern FILE *logfp ;
float
BackpropTrainEpoch(BACKPROP *bp, int ntrials,
                   int (*io_func)(float *inputs, float *targets, int index,
                                  int ninputs, int noutputs, void *user),
                   void *user)
{
  char      *input_tested ;
  float     *inputs, *targets, error_ratio ;
  int       tested, index, min_trial, max_trial, i ;
  BACKPROP  *new_bp = NULL ;

  new_bp = BackpropCopy(bp, NULL) ;
  input_tested = (char *)calloc(ntrials, sizeof(char)) ;
  inputs = (float *)calloc(bp->ninputs, sizeof(float)) ;
  targets = (float *)calloc(bp->noutputs, sizeof(float)) ;
  if (!input_tested | !inputs | !targets)
    ErrorExit(ERROR_BAD_FILE,  "BackpropTrainEpoch(%d): could not allocated arrays\n",
              ntrials) ;

  /* reset training parameters */
  for (i = 0 ; i < bp->noutputs ; i++)
    new_bp->errors[i] = 0.0f ;

  new_bp->ntrials = 0 ;
  min_trial = 0 ;
  max_trial = ntrials - 1 ;
  new_bp->sse = 0.0f ;

  for (tested = 0 ; tested < ntrials ; )   /* go through training set once */
  {
    index = nint(randomNumber(0.0, (double)(max_trial))) ;
    if (index == max_trial)
      while (input_tested[--max_trial] && (max_trial > min_trial))
      {}
    if (index == min_trial)
      while (input_tested[++min_trial] && (max_trial > min_trial))
      {}

    if (input_tested[index])
      continue ;
    tested++ ;
    if (!(*io_func)(inputs, targets, index, bp->ninputs, bp->noutputs, user))
      continue ;
    input_tested[index] = 1 ;

    BackpropLearn(new_bp, inputs, targets) ;  /* train candidate network */
  }

#if 1
  new_bp->sse /= (float)ntrials ;
#endif

  if (!bp->nepochs)
    bp->sse = new_bp->sse ;  /* nothing to compare it to */


  /*
    for adaptive learning rate: If sse increases by more than ERROR_RATIO,
    then update the step size according to:

    learning_rate = learning_rate * TRATE_DECREASE

    otherwise, if the sse decresed, the increase the step size by:

    learning_rate = learning_rate * TRATE_INCREASE.
  */

  error_ratio = bp->sse / new_bp->sse ;

#if 0
  /* only use new weights if error didn't increase too much */
  if (new_bp->sse < bp->sse * bp->error_ratio)
    BackpropCopy(new_bp, bp) ;

  if (error_ratio >= 1.0f)   /* error decreased */
  {
    bp->momentum = bp->old_momentum ;  /* restore original momentum */
  }
  else                       /* error increased */
  {
    bp->momentum = 0.0f ;
    error_ratio *= 2.0f ;    /* decrease more than increase */
  }
  bp->trate = bp->trate * error_ratio ;

#else
  if (new_bp->sse > bp->sse * bp->error_ratio)  /* don't accept change */
  {
    bp->trate = bp->trate * bp->trate_down ;
    bp->momentum = 0.0f ;
  }
  else   /* new error is acceptable, use new network */
  {
    if (new_bp->sse < bp->sse)   /* error decreased, increase trate */
    {
      new_bp->trate = new_bp->trate * new_bp->trate_up ;
    }
    else  /* error increased by a little   NEW!!!! */
    {
      new_bp->trate = new_bp->trate * new_bp->trate_down ;
      /*      new_bp->momentum = 0.0f ;*/
    }
    /* only use new weights if error didn't increase too much */
    BackpropCopy(new_bp, bp) ;

    bp->momentum = bp->old_momentum ;  /* restore original momentum */
  }
  if (bp->trate < MIN_TRATE)
    bp->trate = MIN_TRATE ;
#endif

  bp->nepochs++ ;

  BackpropFree(&new_bp) ;
  free(input_tested) ;
  free(inputs) ;
  free(targets) ;
  return(bp->sse) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
float
BackpropLearn(BACKPROP *backprop, float *inputs,  float *targets)
{
  int     bpClass, i ;
  LAYER   *hidden, *output ;
  float   bpError ;

  if (Gdiag & DIAG_BACKPROP)
  {
    printf("\nBackpropLearn ") ;
    for (i = 0 ; i < backprop->ninputs ; i++)
      printf("%2.3f ", inputs[i]) ;
    printf("--> ") ;
    for (i = 0 ; i  < backprop->noutputs ; i++)
      printf("%2.3f ", targets[i]) ;
    printf("\n") ;
  }

  if (!backprop->hidden.deltas)
  {
    backprop->hidden.deltas =
      (float *)InsCalloc(backprop->nhidden, sizeof(float)) ;
    backprop->output.deltas =
      (float *)InsCalloc(backprop->noutputs, sizeof(float)) ;
  }

  backprop->learn = 1 ;

  bpClass = BackpropProcess(backprop, inputs) ;

  hidden = &backprop->hidden ;
  output = &backprop->output ;

  bpNormalizeTargets(backprop, targets) ;  /* convert targets to -1 --> 1 */
  bpCalculateOutputDeltas(backprop, targets) ;
  bpCalculateHiddenDeltas(backprop) ;
  bpUpdateLayerWeights(hidden, inputs, backprop->trate, backprop->momentum) ;
  bpUpdateLayerWeights(output, hidden->x, backprop->trate, backprop->momentum);

  bpError = BackpropError(backprop, targets) ;

  /*
    moved the next 2 lines to after the error call.
  */
  bpUnnormalizeOutputs(backprop) ;
  bpUnnormalizeTargets(backprop, targets) ;

  backprop->learn = 0 ;
  backprop->ntrials++ ;

  return(bpError) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpCalculateOutputDeltas(BACKPROP *backprop, float *targets)
{
  LAYER           *output ;
  float           *ptarget, *px, *pdelta, target ;
  int             j, nunits ;
  register float  out ;

  /* put in deltas */
  output = &backprop->output ;
  ptarget = &targets[0] ;
  px = &output->x[0] ;
  pdelta = &output->deltas[0] ;
  nunits = output->nunits ;
  for (j = 0 ; j < nunits ; j++)
  {
    target = *ptarget++ ;
    out = *px++ ;
    *pdelta++ = (target - out) /* * fprime(out) */ ;

#if 0
    if (Gdiag & DIAG_BACKPROP)
      printf("Output Delta %d: (%2.3f - %2.5f) * %2.5f = %2.5f\n",
             j, target, out, fprime(out), output->deltas[j]) ;
#endif
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpCalculateHiddenDeltas(BACKPROP *backprop)
{
  LAYER          *output, *hidden ;
  int            i, j, hnunits, onunits ;
  register float *pdeli, *pdelj ;

  /* put in deltas */
  hidden = &backprop->hidden ;
  output = &backprop->output ;

  hnunits = hidden->nunits ;
  onunits = output->nunits ;
  pdeli = &hidden->deltas[0] ;
  for (i = 0 ; i < hnunits ; i++)
  {
    *pdeli = 0.0f ;

#if 0
    if (Gdiag & DIAG_BACKPROP)
      printf("HiddenDelta %d: (", j) ;
#endif

    pdelj = &output->deltas[0] ;
    for (j = 0 ; j < onunits ; j++)
    {
      *pdeli += *pdelj++ * Wij(output, i, j) ;
#if 0
      if (Gdiag & DIAG_BACKPROP)
      {
        printf("%2.5f * %2.5f", output->deltas[j], Wij(output,i,j)) ;
        if (j < output->nunits-1)
          printf(" + ") ;
      }
#endif
    }

    *pdeli++ *= fprime(hidden->x[i]) ;
#if 0
    if (Gdiag & DIAG_BACKPROP)
      printf(") * %2.5f = %2.5f\n",
             fprime(hidden->x[i]), hidden->deltas[i]) ;
#endif
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpInitLayer(LAYER *layer, int ninputs, int nunits)
{
  long  nweights ;

  layer->ninputs = ninputs ;
  layer->nunits = nunits ;
  layer->x = (float *)InsCalloc(nunits, sizeof(*(layer->x))) ;

  nweights = (long)ninputs * (long)nunits ;
  layer->w = (float huge *)InsHalloc(nweights, sizeof(*(layer->w))) ;
  layer->dw = (float huge *)InsHalloc(nweights, sizeof(*(layer->dw))) ;
  layer->biases = (float *)InsCalloc(nunits, sizeof(*(layer->biases))) ;
  layer->db = (float *)InsCalloc(nunits, sizeof(*(layer->db))) ;
  bpInitLayerWeights(layer) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpUpdateLayerWeights(LAYER *layer, float *I, float trate, float momentum)
{
  int            i, j, nunits, ninputs ;
  float          delta, one_minus_momentum, db, *Ii ;
  register float dw, *dwij, *wij ;

  one_minus_momentum = trate * (1.0f - momentum) ;
  nunits = layer->nunits ;
  ninputs = layer->ninputs ;

  for (j = 0 ; j < nunits ; j++)
  {
    delta = layer->deltas[j] ;
    db = delta * one_minus_momentum + momentum * layer->db[j] ;
    if (Gdiag & DIAG_BACKPROP)
      printf("update bias %d: %2.5f + (%2.3f * %2.5f) = %2.5f\n",
             j, layer->biases[j], trate, delta, layer->biases[j]+trate*delta) ;

    layer->db[j] = db ;
    layer->biases[j] += db ;

    dwij = &DWij(layer, 0, j) ;
    wij = &Wij(layer, 0, j) ;
    Ii = &I[0] ;
    for (i = 0 ; i < ninputs ; i++, wij++, dwij++)
    {
      dw = one_minus_momentum * delta * *Ii++ + momentum * *dwij ;
#if 0
      if (Gdiag & DIAG_BACKPROP)
        printf("update weight %d-->%d: %2.5f + (%2.3f * %2.5f * %2.3f)"
               "= %2.5f\n",
               i,j, Wij(layer,i,j), trate, delta, I[i],
               Wij(layer,i,j)+dw) ;
#endif
      *wij += dw ;
      *dwij = dw ;
    }
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpInitLayerWeights(LAYER *layer)
{
  int     i, j ;
  double  rlim ;

  /* Haykin rule-of-thumb for good initial weight range */
  rlim = 2.4 / (double)(layer->ninputs) ;

  /* initialize weights and biases to random values between -1 and 1 */
  for (j = 0 ; j < layer->nunits ; j++)
    layer->biases[j] = (float)randomNumber(-rlim, rlim) ;

  for (j = 0 ; j < layer->nunits ; j++)
  {
    for (i = 0 ; i < layer->ninputs ; i++)
      Wij(layer, i, j) = (float)randomNumber(-rlim, rlim) ;
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpReadLayer(FILE *fp, LAYER *layer)
{
  int i, j ;

  fscanf(fp, "%d  %d\n", &layer->ninputs, &layer->nunits) ;
  bpInitLayer(layer, layer->ninputs, layer->nunits) ;

  fscanf(fp, "\n") ;

  /* now read in weights */
  for (j = 0 ; j < layer->nunits ; j++)
  {
    fscanf(fp, "%f\n", &layer->biases[j]) ;
    for (i = 0 ; i < layer->ninputs ; i++)
      fscanf(fp, "%f ", &Wij(layer, i, j)) ;

    fscanf(fp, "\n") ;
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpWriteLayer(FILE *fp, LAYER *layer)
{
  int i, j ;

  fprintf(fp, "%d  %d\n", layer->ninputs, layer->nunits) ;

  for (j = 0 ; j < layer->nunits ; j++)
  {
    fprintf(fp, "%f\n", layer->biases[j]) ;
    for (i = 0 ; i < layer->ninputs ; i++)
      fprintf(fp, "%f ", Wij(layer, i, j)) ;

    fprintf(fp, "\n") ;
  }
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
BackpropFree(BACKPROP **pbackprop)
{
  BACKPROP *backprop ;

  backprop = *pbackprop ;
  bpFreeLayer(&backprop->hidden) ;
  bpFreeLayer(&backprop->output) ;
  if (backprop->user)
    free(backprop->user) ;
  InsFree(backprop) ;
  *pbackprop = NULL ;
  return(0) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static void
bpFreeLayer(LAYER *layer)
{
  InsFree(layer->x) ;
  InsHfree(layer->w) ;
  InsHfree(layer->dw) ;
  InsFree(layer->biases) ;
  InsFree(layer->db) ;
  if (layer->deltas)
    InsFree(layer->deltas) ;
  layer->biases = NULL ;
  layer->deltas = NULL ;
  layer->x = NULL ;
  layer->w = NULL ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/

#ifndef TINY
#define TINY  0.000001f
#endif

float
BackpropError(BACKPROP *bp, float *targets)
{
  float  bpError, total, error ;
  int    i, nunits ;
  LAYER  *output ;

  output = &bp->output ;
  nunits = output->nunits ;
#if 0
  for (i = 0, total = bpError = 0.0f ; i < nunits ; i++)
  {
    error = fabs(targets[i] - output->x[i]) ;

    bpError += error ;
    if (FZERO(targets[i]))
      target = TINY ;
    else
      target = fabs(targets[i]) ;

    bp->errors[i] += error / target ;
    ;
    total += target ;
  }

  bpError /= total ;
#else
  for (i = 0, total = bpError = 0.0f ; i < nunits ; i++)
  {
    error = targets[i] - output->x[i] ;
    error *= error ;

    bpError += error ;
    bp->errors[i] += error ;
  }

#endif

  bp->sse += bpError ;

  return(bpError) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
BackpropErrorReset(BACKPROP *bp)
{
  memset(bp->errors, 0, bp->noutputs*sizeof(float)) ;

  return(0) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
int
BackpropSetParms(BACKPROP *bp, float trate_up,float trate_down,
                 float error_ratio)
{
  bp->trate_up = trate_up ;
  bp->trate_down = trate_down ;
  bp->error_ratio = error_ratio ;
  return(0) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
#if 0
int
BackpropEpochComplete(BACKPROP *bp)
{
  int  ntrials ;

  ntrials = bp->ntrials ;

  bp->ntrials = 0 ;
  if (bp->nepochs++)  /* nothing to compare 0th epoch to */
  {
    /*
      for adaptive learning rate: If sse increases by more than ERROR_RATIO,
      then update the step size according to:

      learning_rate = learning_rate * TRATE_DECREASE

      otherwise, if the sse decresed, the increase the step size by:

      learning_rate = learning_rate * TRATE_INCREASE.
    */

    bp->momentum = bp->old_momentum ;
    if (new_bp->sse > bp->sse * bp->error_ratio)
    {
      bp->trate = bp->trate * bp->trate_down ;
      bp->momentum = 0.0f ;
    }
    else if (new_bp->sse < bp->sse)
      bp->trate = bp->trate * bp->trate_up ;
  }
  bp->old_sse = bp->sse ;
  bp->sse = 0.0f ;

  return(0) ;
}
#endif
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static long
bpFileNewEnd(FILE *fp, BPFILE_HEADER *hd, int swapped)
{
  long  end ;
  int   err ;

  if (fseek(fp, 0L, SEEK_END))
    ErrorExit(ERROR_BAD_FILE, "bpFileNewEnd: could not seek to end\n") ;

  end = ftell(fp) ;

  bpFileSeekEndPtr(fp, hd, swapped) ;
  DiagPrintf(DIAG_WRITE, "bpFileNewEnd: putting end %ld at %ld\n",
             end, ftell(fp)) ;
  fseek(fp, 0L, SEEK_CUR) ;
  if ((err = fwrite(&end, sizeof(end), 1, fp)) != 1)
    ErrorExit(ERROR_BAD_FILE,
              "bpFileNewEnd: %d, could not write new end location %ld @ %ld\n",
              err, end, ftell(fp)) ;

  return(ftell(fp)) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static long
bpFileSeekEndPtr(FILE *fp, BPFILE_HEADER *hd, int swapped)
{
  long next ;

  /*
    now follow pointer chain to find position of null pointer that we need
    to update with the location of the new net.
  */
  next = (char *)&(hd->first) - (char *)hd ;
  DiagPrintf(DIAG_WRITE, "bpFileSeekEnd: starting search at %ld\n",next) ;
  while (next)
  {
    if (fseek(fp, next, SEEK_SET))
      ErrorExit(ERROR_BAD_FILE, "bpFileSeekEndPtr: could not seek to %ld\n", next) ;
    if (fread(&next, sizeof(next), 1, fp) != 1)
      ErrorExit(ERROR_BAD_FILE, "bpFileNewEnd: could not read pointer @%ld\n",ftell(fp));

    if (swapped)
      next = swapLong32(next) ;
    DiagPrintf(DIAG_WRITE, "bpFileSeekEnd: next %ld @ %ld\n", next,
               ftell(fp)-sizeof(long));
  }

  /* file ptr should be at location of end of pointer chain, not beyond it */
  if (fseek(fp, -(long)sizeof(long), SEEK_CUR))
    ErrorExit(ERROR_BAD_FILE, "bpFileSeekEndPtr: could not seek back sizeof(long)\n") ;

  return(ftell(fp)) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static long
bpFileSeekNet(FILE *fp, BPFILE_HEADER *hd, int netno, int swapped)
{
  int  index ;
  long next ;

  /*
    now follow pointer chain to find position of null pointer that we need
    to update with the location of the new net.
  */
  next = (char *)&(hd->first) - (char *)hd ;
  DiagPrintf(DIAG_WRITE, "bpFileSeekNet(%d): starting search at %ld\n",
             netno, next) ;
  for (index = 0 ; index <= netno ; index++)
  {
    if (fseek(fp, next, SEEK_SET))
      ErrorExit(ERROR_BAD_FILE, "bpFileSeekNet: could not seek to %ld\n", next) ;
    if (fread(&next, sizeof(next), 1, fp) != 1)
      ErrorExit(ERROR_BAD_FILE, "bpFileSeekNet: could not seek to next pointer\n") ;

    if (swapped)
      next = swapLong32(next) ;
    DiagPrintf(DIAG_WRITE, "bpFileSeekNet: %d at %ld\n", index, next) ;
  }
  if (fseek(fp, next+sizeof(long), SEEK_SET))
    ErrorExit(ERROR_BAD_FILE, "bpFileSeekNet: could not seek forward sizeof(long)\n") ;

  return(next+sizeof(long)) ;
}
/*----------------------------------------------------------------------
    Parameters:

    Description:

    Returns:
----------------------------------------------------------------------*/
static int
bpChangeNumberOfNets(FILE *fp, int nnets)
{
  BPFILE_HEADER  hd ;

  DiagPrintf(DIAG_WRITE, "bpChangeNumberOfNets(%d)\n", nnets) ;

  if (fseek(fp, 0L, SEEK_SET))
    ErrorExit(ERROR_BAD_FILE,
              "bpChangeNumberOfNets(%d): could not seek to start of file\n",
              nnets) ;

  if (fread(&hd, sizeof(hd), 1, fp) != 1)
    ErrorExit(ERROR_BAD_FILE, "bpChangeNumberOfNets(%d): could not read header\n",
              nnets) ;

  hd.nnets += nnets ;
  if (fseek(fp, 0L, SEEK_SET))
    ErrorExit(ERROR_BAD_FILE,
              "bpChangeNumberOfNets(%d): could not seek to start of file\n",
              nnets) ;

  if (fwrite(&hd, sizeof(hd), 1, fp) != 1)
    ErrorExit(ERROR_BAD_FILE, "bpChangeNumberOfNets(%d): could not write header\n",
              nnets) ;

  return(0) ;
}


