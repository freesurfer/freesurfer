

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "mri.h"
#include "matrix.h"
#include "xvutil.h"
#include "image.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "diag.h"

static void showSlice(MRI *mri, int slice, int which) ;
static void createPanel(XV_FRAME *xvf) ;
static void sliderNewDepth(Panel_item item, int value, Event *event);
static void sliderNewThreshold(Panel_item item, int value, Event *event);
static void buttonXcorr(Panel_item item, Event *event) ;
static void buttonThreshold(Panel_item item, Event *event) ;
static void buttonPCA(Panel_item item, Event *event) ;
static void repaint_handler(XV_FRAME *xvf, DIMAGE *dimage) ;
static MRI  *register_mri(MRI *mri_in, MRI *mri_ref, MRI *mri_reg) ;
static void draw_basis(MATRIX *m_basis, int x0, int y0, int z0, int which) ;
static IMAGE *show_mri(MRI *mri, int which, int slice) ;
static int get_option(int argc, char *argv[]) ;
void main(int argc, char *argv[]) ;


#define MRI_REFERENCE  "/space/inverse/4/users/inverse/subjects/anders/mri/T1"
#define MRI_INPUT      "/space/inverse/4/users/inverse/subjects/anders/mri/T1"
#define WINDOW_SIZE    10

#define IMAGE_ROWS           2
#define IMAGE_COLS           3
#define BUTTON_ROWS          6
#define DISPLAY_SIZE         210

#define MRI_REF_IMAGE        0
#define MRI_INPUT_IMAGE      1
#define MRI_XCORR_IMAGE      2
#define MRI_THRESH_IMAGE     2
#define MRI_THRESHOLD_IMAGE  (MRI_THRESH_IMAGE)
#define MRI_REG_IMAGE        3
#define MRI_DIFF_IMAGE       4

#define ROW_HEIGHT           30
#define FIRST_BUTTON_ROW     30
#define SECOND_BUTTON_ROW    60
#define THIRD_BUTTON_ROW     90
#define FOURTH_BUTTON_ROW    120
#define FIFTH_BUTTON_ROW     150
#define SIXTH_BUTTON_ROW     180
#define SEVENTH_BUTTON_ROW   210
#define EIGHTH_BUTTON_ROW    240
#define LAST_BUTTON_ROW      EIGHTH_BUTTON_ROW

#define FIRST_BUTTON_COL     5
#define SECOND_BUTTON_COL    60
#define THIRD_BUTTON_COL     165
#define FOURTH_BUTTON_COL    315
#define FIFTH_BUTTON_COL     450
#define SIXTH_BUTTON_COL     565
#define SEVENTH_BUTTON_COL   665

#define FIRST_FNAME_COL      5
#define SECOND_FNAME_COL     365

static IMAGE *Idisplay[IMAGE_ROWS*IMAGE_COLS] = { NULL } ;
char         *Progname ;

static XV_FRAME *xvf ;

static int          current_depth, window_size ;
static Panel_item   depth_slider_panel ;

static unsigned char threshold = 40 ;
static Panel_item   threshold_slider_panel ;

static MRI   *mri_ref, *mri_in, *mri_xcorr, *mri_reg, *mri_diff, 
             *mri_threshold = NULL ;

static MATRIX *m_ref_evectors = NULL, *m_in_evectors = NULL ;
static float  ref_evalues[3], in_evalues[3] ;
static int    ref_means[3], in_means[3] ;
static int    processed = 0 ;


/* options */
static int    dx = 0, dy = 0, dz = 0 ;
static float  xa = 0.0f, ya = 0.0f, za = 0.0f ;

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrUpper(option) ;
  if (!strcmp(option, "DX"))
  {
    sscanf(argv[2], "%d", &dx) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "DY"))
  {
    sscanf(argv[2], "%d", &dy) ;
    nargs = 1 ;
  }
  else if (!strcmp(option, "DZ"))
  {
    sscanf(argv[2], "%d", &dz) ;
    nargs = 1 ;
  }
  else switch (*option)
  {
  case '?':
  case 'U':
    printf("usage: %s [image file name]\n", argv[0]) ;
    exit(1) ;
    break ;
  case 'X':
    sscanf(argv[2], "%f", &xa) ;
    nargs = 1 ;
    fprintf(stderr, "rotating about x axis by %2.1f degrees\n", xa) ;
    xa = RADIANS(xa) ;
    break ;
  case 'Y':
    sscanf(argv[2], "%f", &ya) ;
    nargs = 1 ;
    fprintf(stderr, "rotating about y axis by %2.1f degrees\n", ya) ;
    ya = RADIANS(ya) ;
    break ;
  case 'Z':
    sscanf(argv[2], "%f", &za) ;
    nargs = 1 ;
    fprintf(stderr, "rotating about z axis by %2.1f degrees\n", za) ;
    za = RADIANS(za) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}
void
main(int argc, char *argv[])
{
  char  ref_fname[STR_LEN], in_fname[STR_LEN] ;
  int   drows, dcols ;
  float scale ;
  MRI   *mri_tmp ;
  
  char **av ;
  int  ac, nargs ;

  Progname = argv[0] ;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }


  if (argc > 3)
  {
    sscanf(argv[3], "%d", &window_size) ;
    fprintf(stderr, "using window size %d\n", 2*window_size+1) ;
  }
  else
    window_size = WINDOW_SIZE ;

  if (argc < 3)
    strcpy(in_fname, MRI_INPUT) ;
  else
    strcpy(in_fname, argv[2]) ;

  if (argc < 2)
    strcpy(ref_fname, MRI_REFERENCE) ;
  else
    strcpy(ref_fname, argv[1]) ;

  fprintf(stderr, "reading '%s'...", ref_fname) ;

  mri_ref = MRIread(ref_fname) ;

  if (mri_ref->height >= mri_ref->width)
  {
    scale = (float)mri_ref->width / (float)mri_ref->height ;
    drows = DISPLAY_SIZE ;
    dcols = nint(scale * DISPLAY_SIZE) ;
  }
  else
  {
    scale = (float)mri_ref->height / (float)mri_ref->width ;
    dcols = DISPLAY_SIZE ;
    drows = nint(scale * DISPLAY_SIZE) ;
  }
  xv_init(XV_INIT_ARGC_PTR_ARGV, &argc, argv, NULL);
  xvf = XValloc(IMAGE_ROWS, IMAGE_COLS, BUTTON_ROWS, drows, dcols, 
                Progname, NULL) ;

  XVsetRepaintHandler(repaint_handler) ;

  fprintf(stderr, "done.\nreading '%s'...", in_fname) ;
  fflush(stderr) ;

  if (!strcmp(ref_fname, in_fname))
    mri_in = MRIcopy(mri_ref, NULL) ;
  else
    mri_in = MRIread(in_fname) ;
  fprintf(stderr, "done.\n") ;

  if (!FZERO(za))
  {
    mri_tmp = MRIrotateZ(mri_in, NULL, za) ;
    MRIcopy(mri_tmp, mri_in) ;
    MRIfree(&mri_tmp) ;
  }
  if (!FZERO(ya))
  {
    mri_tmp = MRIrotateY(mri_in, NULL, ya) ;
    MRIcopy(mri_tmp, mri_in) ;
    MRIfree(&mri_tmp) ;
  }
  if (!FZERO(xa))
  {
    mri_tmp = MRIrotateX(mri_in, NULL, xa) ;
    MRIcopy(mri_tmp, mri_in) ;
    MRIfree(&mri_tmp) ;
  }

  if (dx || dy || dz)
  {
    fprintf(stderr, "translating by (%d, %d, %d)\n", dx, dy, dz) ;
    mri_tmp = MRItranslate(mri_in, NULL, -25, 5, 0) ;
    MRIcopy(mri_tmp, mri_in) ;
    MRIfree(&mri_tmp) ;
  }
  createPanel(xvf) ;
  showSlice(mri_ref, current_depth, MRI_REF_IMAGE) ;
  showSlice(mri_in, current_depth, MRI_INPUT_IMAGE) ;
  buttonPCA(0, NULL) ;
  xv_main_loop(xvf->frame) ;


  MRIfree(&mri_ref) ;
  MRIfree(&mri_in) ;
  if (mri_xcorr)
    MRIfree(&mri_xcorr) ;

  exit(0) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
createPanel(XV_FRAME *xvf)
{
  current_depth = (mri_ref->imnr1 + mri_ref->imnr0) / 2 ;
  depth_slider_panel = (Panel_item)
    xv_create(xvf->panel, PANEL_SLIDER,
              PANEL_LABEL_STRING,    "depth",
              PANEL_VALUE,           current_depth,
              PANEL_MIN_VALUE,       mri_ref->imnr0,
              XV_X,                  FIRST_FNAME_COL,
              XV_Y,                  THIRD_BUTTON_ROW,
              PANEL_MAX_VALUE,       mri_ref->imnr1,
              PANEL_SLIDER_WIDTH,    100,
              PANEL_TICKS,           mri_ref->depth,
              PANEL_NOTIFY_PROC,     sliderNewDepth,
              NULL) ;

  threshold_slider_panel = (Panel_item)
    xv_create(xvf->panel, PANEL_SLIDER,
              PANEL_LABEL_STRING,    "threshold",
              PANEL_VALUE,           (int)threshold,
              PANEL_MIN_VALUE,       0,
              XV_X,                  FIRST_FNAME_COL,
              XV_Y,                  FOURTH_BUTTON_ROW,
              PANEL_MAX_VALUE,       255,
              PANEL_SLIDER_WIDTH,    100,
              PANEL_TICKS,           10,
              PANEL_NOTIFY_PROC,     sliderNewThreshold,
              NULL) ;

  xv_create(xvf->panel, PANEL_BUTTON,
            XV_X,                  SECOND_BUTTON_COL,
            XV_Y,                  FIRST_BUTTON_ROW,
            PANEL_LABEL_STRING,    "CORRELATE",
            PANEL_NOTIFY_PROC,      buttonXcorr,
            NULL);

  xv_create(xvf->panel, PANEL_BUTTON,
            XV_X,                  THIRD_BUTTON_COL,
            XV_Y,                  FIRST_BUTTON_ROW,
            PANEL_LABEL_STRING,    "THRESHOLD",
            PANEL_NOTIFY_PROC,      buttonThreshold,
            NULL);

  xv_create(xvf->panel, PANEL_BUTTON,
            XV_X,                  FOURTH_BUTTON_COL,
            XV_Y,                  FIRST_BUTTON_ROW,
            PANEL_LABEL_STRING,    "PCA",
            PANEL_NOTIFY_PROC,      buttonPCA,
            NULL);

}

static void
sliderNewDepth(Panel_item item, int value, Event *event)
{
  current_depth = value ;
  showSlice(mri_ref, current_depth, MRI_REF_IMAGE) ;
  showSlice(mri_in, current_depth, MRI_INPUT_IMAGE) ;
  showSlice(mri_reg, current_depth, MRI_REG_IMAGE) ;
  showSlice(mri_diff, current_depth, MRI_DIFF_IMAGE) ;
}

static void
sliderNewThreshold(Panel_item item, int value, Event *event)
{
  threshold = (BUFTYPE)value ;
}

static void
showSlice(MRI *mri, int slice, int which)
{
  IMAGE **pI, *I ;
  char  *title = "nothing" ;

  switch (which)
  {
  case MRI_REF_IMAGE:
    title = "reference" ;
    break ;
  case MRI_INPUT_IMAGE:
    title = "input" ;
    break ;
  case MRI_XCORR_IMAGE:
    title = "cross-correlation" ;
    break ;
  case MRI_DIFF_IMAGE:
    title = "difference" ;
    break ;
  case MRI_REG_IMAGE:
    title = "registered" ;
    break ;
  }

  pI = Idisplay + which ;
  I = *pI ;
  if (I)
    ImageFree(&I) ;

  I = MRItoImage(mri, NULL, current_depth - mri_ref->imnr0) ;
  *pI = I ;
  XVshowImage(xvf, which, I, 0) ;
  XVshowImageTitle(xvf, which, "%s image (%d)", 
                   title, current_depth - mri_ref->imnr0) ;
}

static void
buttonThreshold(Panel_item item, Event *event)
{
  mri_threshold = MRIthreshold(mri_in, mri_threshold, threshold) ;
  show_mri(mri_threshold, MRI_THRESH_IMAGE, current_depth) ;
  XVshowImageTitle(xvf, MRI_THRESH_IMAGE, "thresholded") ;
}

static void
buttonPCA(Panel_item item, Event *event)
{
  int    row, col ; 
  float  dot ;

  processed = 1 ;
  if (!m_ref_evectors)
    m_ref_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;
  if (!m_in_evectors)
    m_in_evectors = MatrixAlloc(3,3,MATRIX_REAL) ;

  MRIprincipleComponents(mri_ref, m_ref_evectors, ref_evalues, ref_means,
                         threshold);
  MRIprincipleComponents(mri_in, m_in_evectors, in_evalues,in_means,threshold);


#if 0
{
  int i ;

  fprintf(stderr, "ref_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_ref_evectors->rptr[i][1],
            m_ref_evectors->rptr[i][2],
            m_ref_evectors->rptr[i][3]) ;

  fprintf(stderr, "\nin_evectors = \n") ;
  for (i = 1 ; i <= 3 ; i++)
    fprintf(stderr, "\t\t%2.2f    %2.2f    %2.2f\n",
            m_in_evectors->rptr[i][1],
            m_in_evectors->rptr[i][2],
            m_in_evectors->rptr[i][3]) ;
}
#endif
  /* draw mean and projection of eigenvectors */
  draw_basis(m_ref_evectors, ref_means[0], ref_means[1], ref_means[2],
             MRI_REF_IMAGE) ;

  /* check to make sure eigenvectors aren't reversed */
  for (col = 1 ; col <= 3 ; col++)
  {
#if 0
    float theta ;
#endif

    for (dot = 0.0f, row = 1 ; row <= 3 ; row++)
      dot += m_in_evectors->rptr[row][col] * m_ref_evectors->rptr[row][col] ;

    if (dot < 0.0f)
    {
      dot *= -1.0f ;
      for (row = 1 ; row <= 3 ; row++)
        m_in_evectors->rptr[row][col] *= -1.0f ;
    }
#if 0
    theta = acos(dot) ;
    fprintf(stderr, "angle[%d] = %2.1f\n", col, DEGREES(theta)) ;
#endif
  }

  draw_basis(m_in_evectors, in_means[0], in_means[1], in_means[2],
             MRI_INPUT_IMAGE) ;

  mri_reg = register_mri(mri_in, mri_ref, mri_reg) ;
}

static void
buttonXcorr(Panel_item item, Event *event)
{
  int   x, y, z, x0, y0, z0 ;
  float xcorr ;
  
  fprintf(stderr, "doing cross-correlation...") ;

  mri_xcorr = MRIxcorrWindow(mri_ref, mri_in, NULL, 2*window_size+1) ;
  x0 = (mri_xcorr->width-1) / 2 ;
  y0 = (mri_xcorr->height-1) / 2 ;
  z0 = (mri_xcorr->depth-1) / 2 ;
  
  fprintf(stderr, "done.\n") ;

  MRIpeak(mri_xcorr, &x, &y, &z) ;
  if (x >= 0 && y >= 0 && z >= 0)
    xcorr = ((float *)(mri_xcorr->slices[z][y]))[x] ;
  else
    xcorr = -1.0f ;

  x -= x0 ;
  y -= y0 ;
  z -= z0 ;
  mri_reg = MRItranslate(mri_in, mri_reg, x, y, z) ;
  mri_diff = MRIabsdiff(mri_reg, mri_ref, mri_diff) ;
  show_mri(mri_diff, MRI_DIFF_IMAGE, current_depth) ;

  fprintf(stderr, "peak of cross-correlation occurs at (%d, %d, %d)=%2.0f\n", 
          x, y, z, xcorr) ;

  show_mri(mri_xcorr, MRI_XCORR_IMAGE, current_depth) ;
  XVshowImageTitle(xvf, MRI_XCORR_IMAGE,"cross correlation") ;

  show_mri(mri_diff, MRI_DIFF_IMAGE, current_depth) ;
  XVshowImageTitle(xvf, MRI_DIFF_IMAGE, "difference image") ;
}
static void
repaint_handler(XV_FRAME *xvf, DIMAGE *dimage)
{
  int        arrow_len ;
  static int entered = 0 ;  /* prevent endless recursion */

  if (processed && !entered)
  {
    entered = 1 ;

    arrow_len = mri_ref->width/4 ;

    draw_basis(m_ref_evectors, ref_means[0], ref_means[1], ref_means[2],
               MRI_REF_IMAGE) ;
    draw_basis(m_in_evectors, in_means[0], in_means[1], in_means[2],
               MRI_INPUT_IMAGE) ;
    draw_basis(m_ref_evectors, ref_means[0], ref_means[1], ref_means[2],
               MRI_REG_IMAGE) ;
    entered = 0 ;
  }
}

static MRI *
register_mri(MRI *mri_in, MRI *mri_ref, MRI *mri_reg)
{
  int    dx, dy, dz ;
  MATRIX *mRot, *m_in_T, *mOrigin ;
  MRI    *mri_tmp ;
  float  x_angle, y_angle, z_angle, r11, r21, r31, r32, r33, cosy ;

  m_in_T = MatrixTranspose(m_in_evectors, NULL) ;
  mRot = MatrixMultiply(m_ref_evectors, m_in_T, NULL) ;

  r11 = mRot->rptr[1][1] ;
  r21 = mRot->rptr[2][1] ;
  r31 = mRot->rptr[3][1] ;
  r32 = mRot->rptr[3][2] ;
  r33 = mRot->rptr[3][3] ;
  y_angle = atan2(-r31, sqrt(r11*r11+r21*r21)) ;
  cosy = cos(y_angle) ;
  z_angle = atan2(r21 / cosy, r11 / cosy) ;
  x_angle = atan2(r32 / cosy, r33 / cosy) ;
  fprintf(stderr, "rotation: (%2.2f, %2.2f, %2.2f)\n",
          DEGREES(x_angle), DEGREES(y_angle), DEGREES(z_angle)) ;
  mOrigin = VectorAlloc(3, MATRIX_REAL) ;
  mOrigin->rptr[1][1] = in_means[0] ;
  mOrigin->rptr[2][1] = in_means[1] ;
  mOrigin->rptr[3][1] = in_means[2] ;

  dx = nint(ref_means[0] - in_means[0]) ;
  dy = nint(ref_means[1] - in_means[1]) ;
  dz = nint(ref_means[2] - in_means[2]) ;

  XVprintf(xvf, 0, "translation parms: (%d, %d, %d)", dx, dy, dz) ;
  mri_tmp = MRIrotate(mri_in, NULL, mRot, mOrigin) ;
  mri_reg = MRItranslate(mri_tmp, mri_reg, dx, dy, dz) ;

  show_mri(mri_tmp, MRI_XCORR_IMAGE, current_depth) ;
  XVshowImageTitle(xvf, MRI_XCORR_IMAGE, "rotated image") ;
  XVdrawBox(xvf, MRI_XCORR_IMAGE, mri_tmp->width/2-1, 
            mri_tmp->height/2-1, 3, 3, XGREEN) ;

  show_mri(mri_reg, MRI_REG_IMAGE, current_depth) ;
  XVshowImageTitle(xvf, MRI_REG_IMAGE, "registered image") ;

  draw_basis(m_ref_evectors, ref_means[0], ref_means[1], ref_means[2],
             MRI_REG_IMAGE) ;

  mri_diff = MRIabsdiff(mri_reg, mri_ref, mri_diff) ;
  show_mri(mri_diff, MRI_DIFF_IMAGE, current_depth) ;
  XVshowImageTitle(xvf, MRI_DIFF_IMAGE, "difference image") ;


  if (mri_tmp)
    MRIfree(&mri_tmp) ;
  MatrixFree(&mRot) ;
  VectorFree(&mOrigin) ;
  return(mri_reg) ;
}
static void
draw_basis(MATRIX *m_basis, int x0, int y0, int z0, int which)
{
  int arrow_len ;

  arrow_len = mri_ref->width/4 ;
  XVdrawBox(xvf, which, x0-2, y0-2, 5, 5, XGREEN) ;

  XVdrawArrow(xvf, which, x0, y0, nint(m_basis->rptr[1][1]*arrow_len), 
              nint(m_basis->rptr[2][1]*arrow_len), XCYAN) ;
  XVdrawArrow(xvf, which, x0, y0, nint(m_basis->rptr[1][2]*arrow_len), 
              nint(m_basis->rptr[2][2]*arrow_len), XRED) ;
}

static IMAGE *
show_mri(MRI *mri, int which, int slice)
{
  if (Idisplay[which])
    ImageFree(&Idisplay[which]) ;

  Idisplay[which] = MRItoImage(mri, NULL, slice - mri->imnr0) ;
  XVshowImage(xvf, which, Idisplay[which], 0) ;
  return(Idisplay[which]) ;
}

