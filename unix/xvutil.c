/*----------------------------------------------------------------------
           File Name:

           Author:

           Description:

           Conventions:

----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
                           INCLUDE FILES
----------------------------------------------------------------------*/

#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "image.h"
#include "utils.h"
#include "macros.h"
#include "xvutil.h"
#include "proto.h"
#include "hipsh.h"
#include "error.h"

/*----------------------------------------------------------------------
                           MACROS AND CONSTANTS
----------------------------------------------------------------------*/

#define MAX_TITLE_CHARS(d)   (MIN(STR_LEN-1, \
                        (nint((float)d->dispImage->cols / (float)CHAR_WIDTH))))

#define MAX_DISP_SCALES   4

#ifdef Linux
#define FRAME_X             150
#else
#define FRAME_X             250
#endif
#define FRAME_Y             10

#define DISPLAY_SIZE        128
#ifdef Linux
#define MIN_FRAME_WIDTH     235
#else
#define MIN_FRAME_WIDTH     220
#endif


#define PANEL_HEIGHT        ((xvf->button_rows+1) * ROW_HEIGHT+4*WINDOW_PAD)

#define CHAR_HEIGHT          20
#define CHAR_PAD             5

#define FIRST_BUTTON_COL     5
#define SECOND_BUTTON_COL    60
#define THIRD_BUTTON_COL     200
#define FOURTH_BUTTON_COL    315
#define FIFTH_BUTTON_COL     450
#define SIXTH_BUTTON_COL     565
#define SEVENTH_BUTTON_COL   665

#define FIRST_FNAME_COL      5
#define SECOND_FNAME_COL     400

#define HIPS_CMD_ROWS       30
#define HIPS_CMD_COLS       400

#define DEFAULT_PRECISION   4
/*----------------------------------------------------------------------
                              PROTOTYPES
----------------------------------------------------------------------*/

static void  xvCreateFrame(XV_FRAME *xvf, char *name) ;
static void  xvInitColors(XV_FRAME *xvf) ;
static void xvInitImages(XV_FRAME *xvf) ;
static void xv_dimage_repaint(Canvas canvas, Xv_Window window, 
                                       Rectlist *repaint_area) ;
static DIMAGE *xvGetDimage(XV_FRAME *xvf, int which, int alloc) ;
static void xv_dimage_event_handler(Xv_Window window, Event *event) ;
static XImage *xvCreateXimage(XV_FRAME *xvf, IMAGE *image) ;
static Panel_setting xvHipsCommand(Panel_item item, Event *event) ;
static void xvHipsCmdFrameInit(XV_FRAME *xvf) ;
static void xvCreateImage(XV_FRAME *xvf, DIMAGE *dimage, int x, int y, 
                          int which) ;
static void xvFreeDimage(DIMAGE *dimage) ;
static Panel_setting xvFileNameCommand(Panel_item item, Event *event) ;
static char *xvGetTitle(int which, char *title, int with_value) ;

/*----------------------------------------------------------------------
                              GLOBAL DATA
----------------------------------------------------------------------*/
static int debug_it = 0 ;

static XV_FRAME *xvf ;
static IMAGE *GtmpFloatImage = NULL, *GtmpByteImage = NULL,
  *GtmpByteImage2 = NULL, *GtmpScaledFloatImage = NULL ;

static Frame           hips_cmd_frame ;
static Display        *hips_cmd_display;
static Panel_item      hips_cmd_panel_item ;
static Panel           hips_cmd_panel ;
static char            hips_cmd_str[301] ;
static int             hips_cmd_source = 0 ;
static void            (*XVevent_handler)(Event *event, DIMAGE *dimage) = NULL;
static void            (*XVkb_handler)(Event *event, DIMAGE *dimage) = NULL;
static void            (*XVquit_func)(void) = NULL;
static void            (*XVrepaint_handler)(XV_FRAME *xvf,DIMAGE *dimage)=NULL;

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
XV_FRAME *
XValloc(int rows, int cols, int button_rows, int display_rows, 
        int display_cols, char *name, Notify_value (*poll)(void))
{
  int               row ;
  static struct itimerval  timer;       

  xvf = (XV_FRAME *)calloc(1, sizeof(XV_FRAME)) ;
  if (!xvf)
    return(NULL) ;

  xvf->noprint = 0 ;
  xvf->ydir = -1 ;
  xvf->precision = DEFAULT_PRECISION ;
  xvf->rows = rows ;
  xvf->cols = cols ;
  xvf->button_rows = button_rows ;
  if (!display_rows)
    display_rows = DISPLAY_SIZE ;
  if (!display_cols)
    display_cols = display_rows ;
  xvf->display_rows = display_rows ;
  xvf->display_cols = display_cols ;

  xvf->dimages = (DIMAGE **)calloc(rows, sizeof(DIMAGE *)) ;
  if (!xvf->dimages)
    return(NULL) ;
  for (row = 0 ; row < rows ; row++)
  {
    xvf->dimages[row] = (DIMAGE *)calloc(cols, sizeof(DIMAGE)) ;
    if (!xvf->dimages[row])
      return(NULL) ;
  }

  xvCreateFrame(xvf, name) ;

  window_fit(xvf->frame);

  if (poll)
  {
    timer.it_value.tv_usec = 1;
    timer.it_interval.tv_usec = 1;
    notify_set_itimer_func(xvf->frame, poll, ITIMER_REAL, &timer, NULL);
    if (notify_errno)
    {
      fprintf(stderr, "notifier error %d\n", notify_errno) ;
      notify_perror("notify error installing poll routine:") ;
    }
  }

  xvHipsCmdFrameInit(xvf) ;

  xvInitImages(xvf) ;

  return(xvf) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
xvCreateFrame(XV_FRAME *xvf, char *name)
{
  int width, height ;

  width = xvf->cols * xvf->display_cols + (xvf->cols-1)*WINDOW_PAD ;
  height = PANEL_HEIGHT + 
    xvf->rows * xvf->display_rows + (xvf->rows-1) * CHAR_HEIGHT ;
  if (width < MIN_FRAME_WIDTH)
    width = MIN_FRAME_WIDTH ;

  xvf->frame = (Frame)xv_create((Xv_opaque)NULL, FRAME,
                           FRAME_LABEL, name,
                           XV_X,        FRAME_X,
                           XV_Y,        FRAME_Y,
                           XV_WIDTH,    width,
                           XV_HEIGHT,   height,
                           NULL);
  if (!xvf->frame)
    ErrorExit(ERROR_BADPARM, 
              "xvCreateFrame(%s, (%d, %d), (%d, %d)): could not create frame",
              name, xvf->rows, xvf->cols,xvf->display_rows,xvf->display_cols);

  xvf->display = (Display *)xv_get(xvf->frame,XV_DISPLAY);
  xvf->screen = DefaultScreen(xvf->display);
  xvf->black_pixel = BlackPixel(xvf->display, xvf->screen);
  xvf->white_pixel = WhitePixel(xvf->display, xvf->screen);

  xvInitColors(xvf) ;

  xvf->panel =(Panel)xv_create((Xv_opaque)xvf->frame,PANEL,XV_X,0,XV_Y,0,NULL);
  if (!xvf->panel)
    ErrorExit(ERROR_BADPARM, 
              "xvCreateFrame(%s, (%d, %d), %d): could not create panel",
              name, xvf->rows, xvf->cols, xvf->display_rows) ;

  xvf->gc = DefaultGC(xvf->display,DefaultScreen(xvf->display));
  XSetForeground(xvf->display, xvf->gc, xvf->black_pixel);
  XSetBackground(xvf->display, xvf->gc, xvf->white_pixel);

  xv_create((Xv_opaque)xvf->panel, PANEL_BUTTON,
            XV_X,                  3,
            XV_Y,                  FIRST_BUTTON_ROW,
            PANEL_LABEL_STRING,    "Quit",
            PANEL_NOTIFY_PROC,     buttonQuit,
            NULL);

  /* create a region for messages */
  xvf->msg_item[0] = (Panel_item) 
    xv_create((Xv_opaque)xvf->panel,        PANEL_MESSAGE,
             PANEL_LABEL_BOLD,   TRUE,
             XV_X,               5,
             XV_Y,               3,
             PANEL_LABEL_STRING, 
             xvf->msg_str[0],
             NULL);
  xvf->msg_item[1] = (Panel_item) 
    xv_create((Xv_opaque)xvf->panel,        PANEL_MESSAGE,
             PANEL_LABEL_BOLD,   TRUE,
             PANEL_LABEL_STRING, xvf->msg_str[1],
             XV_X,               SEVENTH_BUTTON_COL,
             XV_Y,               5,
             NULL);

}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVprintf(XV_FRAME *xvf, int which, ...)
{
  va_list  args ;
  char     *fmt ;
  int      len ;

  if (which < 0 || which > 1)
    ErrorReturn(0, (ERROR_BADPARM, "XVprintf(%d): out of range", which)) ;

  va_start(args, which) ;
  fmt = va_arg(args, char *) ;
  vsprintf(xvf->msg_str[which], fmt, args) ;
  len = strlen(xvf->msg_str[which]) ;
  xv_set(xvf->msg_item[which], PANEL_LABEL_STRING, xvf->msg_str[which], NULL);

  va_end(args) ;
  XSync(xvf->display, 0) ;
  return(len) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/

#define MAX_DISP_VAL     MAX_GRAY
#define RESERVED_COLORS  6
#define MAX_GRAY         (MAX_COLORS-RESERVED_COLORS-1)


#define PURPLE_INDEX     (MAX_COLORS-6)
#define YELLOW_INDEX     (MAX_COLORS-5)
#define CYAN_INDEX       (MAX_COLORS-4)
#define GREEN_INDEX      (MAX_COLORS-3)
#define BLUE_INDEX       (MAX_COLORS-2)
#define RED_INDEX        (MAX_COLORS-1)

static XColor        xcolors[MAX_COLORS] ;

static void
xvInitColors(XV_FRAME *xvf)
{
  int           i, shift ;
  unsigned long *pix_vals ;

  shift = MAX_COLORS == 32 ? 11 : MAX_COLORS == 64 ? 10 : 
    MAX_COLORS == 128 ? 9 : 8 ;

  for (i = 0 ; i < MAX_COLORS ; i++)
  {
    xcolors[i].pixel = i ;
    switch (i)
    {
#if 1
    case PURPLE_INDEX:
      xcolors[i].red = xcolors[i].blue = 255<<shift ;
      xcolors[i].green = 200<<shift ;
      break ;
    case YELLOW_INDEX:
      xcolors[i].green = xcolors[i].red = 200<<shift ;
      xcolors[i].blue = 0 ;
      break ;
    case CYAN_INDEX:
      xcolors[i].green = xcolors[i].blue = 255<<shift ;
      xcolors[i].red = 128<<shift ;
      break ;
    case RED_INDEX:
      xcolors[i].green = xcolors[i].blue = 200<<shift ;
      xcolors[i].red = 255<<shift ;
      break ;
    case BLUE_INDEX:
      xcolors[i].green = xcolors[i].red = 200<<shift ;
      xcolors[i].blue = 255<<shift ;
      break ;
    case GREEN_INDEX:
      xcolors[i].red = xcolors[i].blue = 100<<shift ;
      xcolors[i].red = xcolors[i].blue = 75<<shift ;
      xcolors[i].green = 255<<shift ;
      break ;
#endif
    default:
      xcolors[i].red = i<<shift ;
      xcolors[i].green = i<<shift ;
      xcolors[i].blue = i<<shift ;
      break ;
    }
  }

  xvf->cms = (Cms)xv_create(xvf->screen, CMS,
                       CMS_NAME,     "xvutil",
                       CMS_SIZE,      MAX_COLORS,
                       CMS_X_COLORS,  xcolors,
                       CMS_TYPE,      MAX_COLORS == 64 ?
                                      XV_STATIC_CMS : XV_DYNAMIC_CMS,
                       NULL) ;
  
  pix_vals = (unsigned long *)xv_get(xvf->cms, CMS_INDEX_TABLE) ;
  xvf->red_pixel = pix_vals[RED_INDEX] ;
  xvf->blue_pixel = pix_vals[BLUE_INDEX] ;
  xvf->green_pixel = pix_vals[GREEN_INDEX] ;
  xvf->cyan_pixel = pix_vals[CYAN_INDEX] ;
  xvf->purple_pixel = pix_vals[PURPLE_INDEX] ;
  xvf->yellow_pixel = pix_vals[YELLOW_INDEX] ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVbrighten(XV_FRAME *xvf, int which, int offset)
{
  DIMAGE        *dimage ;

  /* global is a hack, but too much trouble otherwise */
  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return(ERROR_BAD_PARM) ;

  dimage->bshift += offset ;
  XVshowImage(xvf, which, dimage->sourceImage, dimage->frame) ;

  return(NO_ERROR) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
xvInitImages(XV_FRAME *xvf)
{
  DIMAGE     *dimage ;
  int        row, col, x, y ;

  y = PANEL_HEIGHT ;
  for (row = 0 ; row < xvf->rows ; row++)
  {
    x = 0 ;
    for (col = 0 ; col < xvf->cols ; col++)
    {
      dimage = &xvf->dimages[row][col] ;
      xvCreateImage(xvf, dimage, x, y, row * xvf->cols + col) ;
      x += WINDOW_PAD + xvf->display_cols ;
    }
    y += CHAR_HEIGHT + xvf->display_rows ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void 
XVclearImageTitle(XV_FRAME *xvf, int which) 
{
  XVshowImageTitle(xvf, which, "                                       ") ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void 
XVclearImage(XV_FRAME *xvf, int which, int dotitle)
{
  DIMAGE    *dimage ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return ;
  XClearArea(xvf->display, dimage->window, 0, 0, 0, 0, False) ;
  if (dotitle)
    XVclearImageTitle(xvf, which) ;
  dimage->used = 0 ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVrepaintImage(XV_FRAME *xvf, int which)
{
  DIMAGE    *dimage ;
  IMAGE *image ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage || dimage->entered)
    return ;
  dimage->entered = 1 ;
  image = dimage->dispImage ;
  if (!dimage->entered)
    ImageWrite(image, "repaint.hipl") ;
  XPutImage(xvf->display, (Drawable)dimage->window, xvf->gc, dimage->ximage, 
            0, 0, 0, 0, image->cols, image->rows);
  if (XVrepaint_handler)
    (*XVrepaint_handler)(xvf, dimage) ;
  dimage->entered = 0 ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVshowImage(XV_FRAME *xvf, int which, IMAGE *image, int frame)
{
  float         xscale, yscale, scale, fmin, fmax ;
  DIMAGE        *dimage ;
  unsigned long *substtable ;
  int           i, rows, cols, srows, scols, c ;
  byte          lut[MAX_COLORS] ;

  if (!image)
    return ;

  dimage = xvGetDimage(xvf, which, 1) ;
  if (!dimage)
    return ;
  dimage->used = 1 ;
  dimage->frame = frame ;
  if ((image != dimage->sourceImage) || (dimage->dx == 0))
  {
    if (dimage->zoomImage != dimage->sourceImage)
      ImageFree(&dimage->zoomImage) ;
    dimage->sourceImage = dimage->zoomImage = image ;
  }
  else   /* image from here on in refers to the zoomed image */
  {
    dimage->zoomImage = 
      ImageExtract(image, NULL, dimage->x0, dimage->y0, dimage->dx,dimage->dy);
    dimage->sourceImage = image ;
    image = dimage->zoomImage ;
  }

  scols = dimage->dispImage->cols ;
  srows = dimage->dispImage->rows ;
  xscale = (float)scols / (float)image->cols ;
  yscale = (float)srows / (float)image->rows ;
  dimage->xscale = xscale ;
  dimage->yscale = yscale ;

  if (!ImageCheckSize(image, GtmpFloatImage, image->rows, image->cols, 1))
  {
    if (GtmpFloatImage)
      ImageFree(&GtmpFloatImage) ;
    GtmpFloatImage = ImageAlloc(image->rows,image->cols, PFFLOAT,1);
  }
  else
    ImageSetSize(GtmpFloatImage, image->rows, image->cols) ;

  if (!ImageCheckSize(dimage->dispImage, GtmpByteImage, srows, scols, 1))
  {
    if (GtmpByteImage)
      ImageFree(&GtmpByteImage) ;
    if (GtmpByteImage2)
      ImageFree(&GtmpScaledFloatImage) ;
    if (GtmpScaledFloatImage)
      ImageFree(&GtmpByteImage2) ;
    GtmpScaledFloatImage = ImageAlloc(srows, scols, PFFLOAT, 1) ;
    GtmpByteImage = ImageAlloc(srows, scols, PFBYTE, 1) ;
    GtmpByteImage2 = ImageAlloc(srows, scols, PFBYTE, 1) ;
  }
  else
  {
    ImageSetSize(GtmpScaledFloatImage, srows, scols) ;
    ImageSetSize(GtmpByteImage, srows, scols) ;
    ImageSetSize(GtmpByteImage2, srows, scols) ;
  }

  if (image->pixel_format == PFDOUBLE)
  {
    IMAGE *Itmp ;

    Itmp = ImageAlloc(image->rows, image->cols, PFFLOAT, image->num_frame) ;
    ImageCopy(image, Itmp) ;
    ImageCopyFrames(Itmp, GtmpFloatImage, frame, 1, 0) ;
    ImageFree(&Itmp) ;
  }
  else
    ImageCopyFrames(image, GtmpFloatImage, frame, 1, 0) ;

  /* scale range of values to be in byte display range (about 0-60) */
  if (!xvf->rescale)
  {
    if (dimage->rescale_range || image->num_frame == 1)
      ImageScale(GtmpFloatImage, GtmpFloatImage, 0, MAX_DISP_VAL) ;
    else   /* use entire sequence to compute display range */
    {
      ImageValRange(image, &fmin, &fmax) ;
      ImageScaleRange(GtmpFloatImage, fmin, fmax, 0, MAX_DISP_VAL) ;
    }
  }
  else   /* find min and max of all synced images of the same pixel format */
  {
    DIMAGE *dimage2 ;
    float  sync_fmin, sync_fmax ;

    ImageValRange(image, &sync_fmin, &sync_fmax) ;
    for (i = 0 ; i < xvf->rows*xvf->cols ; i++)
    {
      dimage2 = xvGetDimage(xvf, i, 0) ;
      if (dimage2&&(image->pixel_format == dimage2->sourceImage->pixel_format))
      {
        ImageValRange(dimage2->sourceImage, &fmin, &fmax) ;
        if (fmin < sync_fmin)
          sync_fmin = fmin ;
        if (fmax > sync_fmax)
          sync_fmax = fmax ;
      }
    }
    ImageScaleRange(GtmpFloatImage, sync_fmin, sync_fmax, 0, MAX_DISP_VAL) ;
  }

  /* resize image to display area */
  scale = (xscale + yscale)/2.0f ;
  rows = nint((float)GtmpFloatImage->rows*scale) ;
  cols = nint((float)GtmpFloatImage->cols*scale) ;

#if 0
  if ((rows != dimage->dispImage->rows) || (cols != dimage->dispImage->cols))
    ImageResize(GtmpFloatImage, GtmpScaledFloatImage, srows, scols) ;
  else
    ImageRescale(GtmpFloatImage, GtmpScaledFloatImage, MIN(xscale,yscale)) ;
#else
  ImageDifferentialScale(GtmpFloatImage, GtmpScaledFloatImage, srows, scols) ;
#endif

  ImageCopy(GtmpScaledFloatImage, GtmpByteImage) ; /* convert to bytes */
  h_invert(GtmpByteImage, dimage->dispImage) ;

  /* use current colormap */
  substtable = (unsigned long *) xv_get(xvf->cms,CMS_INDEX_TABLE);
  for (i=0;i<MAX_COLORS;i++)
  {
    c = (dimage->gamma[i] + dimage->bshift)  ;
    if (c < 0)
      c = 0 ;
    if (c > MAX_GRAY)
      c = MAX_GRAY ;
    lut[i] = (byte)substtable[c];
  }
#if 0
  if (MAX_COLORS == 64)
    h_shift_b(dimage->dispImage, dimage->dispImage, -2);
#endif
  h_applylut(dimage->dispImage, dimage->dispImage, MAX_COLORS, lut);

  XVrepaintImage(xvf, which) ;
  XFlush(xvf->display); 
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
xv_dimage_repaint(Canvas canvas, Xv_Window window, Rectlist *repaint_area)
{
  int    row, col, which = -1 ;
  DIMAGE *dimage ;

  dimage = NULL ;
  for (row = 0 ; row < xvf->rows ; row++)
  {
    for (col = 0 ; col < xvf->cols ; col++)
    {
      if (xvf->dimages[row][col].canvas == canvas)
      {
        dimage = &xvf->dimages[row][col] ;
        which = row * xvf->cols + col ;
        break ;
      }
    }
    if (which >= 0)
      break ;
  }
  if (dimage && dimage->used)
    XVrepaintImage(xvf, which) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static DIMAGE *
xvGetDimage(XV_FRAME *xvf, int which, int alloc)
{
  DIMAGE *dimage ;
  int    row, col ;

  row = which / xvf->cols ;
  col = which % xvf->cols ;
  if (row >= xvf->rows || col >= xvf->cols)
    return(NULL) ;

  dimage = &xvf->dimages[row][col] ;

  if (!alloc && !dimage->used)
    return(NULL) ;

  return(dimage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
/*
  the xview coordinate system is (0,0) in the upper left, and increases
  down and to the right.
*/

static void
xv_dimage_event_handler(Xv_Window xv_window, Event *event)
{
  int        x, y, yprint, i ;
  double     val = 0.0 ;
  Window     window ;
  int        row, col, which = -1, rows, cols ;
  DIMAGE     *dimage, *dimage2 = NULL ;
  char       fmt[150], buf[100], title[40] ;
  float      aspect = 1.0f, scale = 1.0f ;

  window = event_window(event) ;

  dimage = NULL ;
  for (row = 0 ; row < xvf->rows ; row++)
  {
    for (col = 0 ; col < xvf->cols ; col++)
    {
      if (canvas_paint_window(xvf->dimages[row][col].canvas) == xv_window)
      {
        dimage = &xvf->dimages[row][col] ;
        which = row * xvf->cols + col ;
      }
    }
  }
  x = event_x(event) ;
  y = event_y(event) ;  
  
  if (!dimage)
  {
    fprintf(stderr, "(%d, %d) in win %ld\n", x, y, event_window(event)) ;
    fprintf(stderr, "could not find appropriate window in event!\n") ;
    return ;
  }
  if (!dimage->used)
    return ;

  x = (int)((float)x / dimage->xscale) ;
  y = (int)((float)y / dimage->yscale) ;

  
  y = (dimage->zoomImage->rows-1) - y ;/* convert y to hips coordinate sysem */

  x += dimage->x0 ;
  y += dimage->y0 ;

  /* do boundary checking */
  if (y < 0) y = 0 ;
  if (x < 0) x = 0 ;
  if (y >= dimage->sourceImage->rows) y = dimage->sourceImage->rows - 1 ;
  if (x >= dimage->sourceImage->cols) x = dimage->sourceImage->cols - 1 ;

  if (xvf->ydir >= 0)
    yprint = (dimage->sourceImage->rows-1)-y ;
  else
    yprint = y ;

  if (event_ctrl_is_down(event)) 
  {
    switch (event_id(event)) 
    {
    case MS_RIGHT:
      if (event_is_up(event))
      {
        dimage->dx1 = dimage->dy1 = dimage->dx = dimage->dy = 
          dimage->x0 = dimage->y0 = dimage->x1 = dimage->y1 = 0 ;
        rows = xvf->display_rows ;
        cols = xvf->display_cols ;
        aspect = 
          (float)dimage->sourceImage->rows / (float)dimage->sourceImage->cols ;
        if (rows > cols)
          XVsetImageSize(xvf, which, rows, nint((float)rows / aspect)) ;
        else
          XVsetImageSize(xvf, which, nint((float)cols * aspect), cols) ;
        if (dimage->sync)
          XVdoSync(xvf, which) ;
        XVresize(xvf) ;
        XVshowAllSyncedImages(xvf, which) ;
      }
      break ;
    case MS_LEFT:
      if (event_is_down(event))
      {
        dimage->x1 = x ;
        dimage->y1 = y ;
      }
      else if (event_is_up(event))  /* erase box and show zoomed image */
      {
        XVdrawBox(xvf, which, dimage->x1, dimage->y1, dimage->dx1,
                  dimage->dy1, XXOR) ;
        if (dimage->dx1 < 0)
        {
          dimage->x0 = dimage->x1 + dimage->dx1 - 1 ;
          dimage->dx1 *= -1 ;
        }
        else
          dimage->x0 = dimage->x1 ;
        if (dimage->dy1 < 0)
        {
          dimage->y0 = dimage->y1 + dimage->dy1 - 1 ;
          dimage->dy1 *= -1 ;
        }
        else
          dimage->y0 = dimage->y1 ;

        if (dimage->x0 + dimage->dx1 >= dimage->sourceImage->cols)
          dimage->dx1 = dimage->sourceImage->cols - dimage->x0 - 1 ;
        if (dimage->y0 + dimage->dy1 >= dimage->sourceImage->rows)
          dimage->dy1 = dimage->sourceImage->rows - dimage->y0 - 1 ;
        if ((dimage->dx1 >= 2) && (dimage->dy1 >= 2))
        {
          float rscale, cscale ;

          rows = xvf->display_rows ;
          cols = xvf->display_cols ;
          rscale = (float)rows / (float)dimage->dy1 ;
          cscale = (float)cols / (float)dimage->dx1 ;
          aspect = (float)dimage->dy1 / (float)dimage->dx1 ;
          dimage->dx = dimage->dx1 ;
          dimage->dy = dimage->dy1 ;
          if (dimage->sync)
            XVdoSync(xvf, which) ;
          if (rscale < cscale)
            XVsetImageSize(xvf, which, rows, nint((float)rows / aspect)) ;
          else
            XVsetImageSize(xvf, which, nint((float)cols * aspect),cols) ;
          XVshowAllSyncedImages(xvf, which) ;
        }
          
        dimage->x1 = dimage->y1 = 0 ;  /* no longer drawing */
        dimage->dx1 = dimage->dy1 = 0 ;
      }
      break ;
    case LOC_DRAG:
      if (event_left_is_down(event))  /* draw rubber band box */
      {
        if (dimage->x1)   /* erase old box */
          XVdrawBox(xvf, which, dimage->x1, dimage->y1, dimage->dx1,
                    dimage->dy1, XXOR) ;
        dimage->dx1 = x - dimage->x1 ;
        dimage->dy1 = y - dimage->y1 ;

        /* enforce appropriate aspect ratio */
#if 0
        aspect = 
          (float)dimage->sourceImage->cols / (float)dimage->sourceImage->rows;
        dimage->dx1 = ISIGN(dimage->dx1) * abs(nint(aspect * dimage->dy1)) ;
#endif
        XVdrawBox(xvf, which, dimage->x1, dimage->y1, dimage->dx1,
                  dimage->dy1, XXOR) ;
      }
      break ;
    default:
      break ;
    }
  }
  else switch (event_id(event)) 
  {
#ifdef Linux
  case MS_MIDDLE:
    xv_set(hips_cmd_frame, FRAME_CMD_PUSHPIN_IN, TRUE, XV_SHOW, TRUE, NULL) ;
    hips_cmd_source = which ;
    break ;
  case MS_RIGHT:
    if (event_is_down(event))
    {
    }
    break ;
#else
  case MS_RIGHT:
    xv_set(hips_cmd_frame, FRAME_CMD_PUSHPIN_IN, TRUE, XV_SHOW, TRUE, NULL) ;
    hips_cmd_source = which ;
    break ;
  case MS_MIDDLE:
    if (event_is_down(event))
    {
    }
    break ;
#endif
  case LOC_DRAG:
    if (!event_left_is_down(event))
      break ;
  case MS_LEFT:
    switch (dimage->sourceImage->pixel_format)
    {
    case PFDOUBLE:
      val = *IMAGEDseq_pix(dimage->sourceImage, x, y, dimage->frame) ;
      break ;
    case PFFLOAT:
      val = (double)*IMAGEFseq_pix(dimage->sourceImage, x, y, dimage->frame);
      break ;
    case PFBYTE:
      val = (double)*IMAGEseq_pix(dimage->sourceImage, x, y, dimage->frame) ;
      break ;
    case PFINT:
      val = (double)*IMAGEIseq_pix(dimage->sourceImage, x, y, dimage->frame) ;
      break ;
    }

    xvGetTitle(which, title, 0) ;
    sprintf(fmt, "%%10.10s: (%%3d, %%3d) --> %%2.%dlf\n", xvf->precision) ;
    if (!xvf->noprint)
      XVprintf(xvf, 0, fmt, title, x, yprint, val) ;
    if (dimage->sync)
    {
      for (i = 0 ; i < xvf->rows*xvf->cols ; i++)
      {
        dimage2 = xvGetDimage(xvf, i, 0) ;
        if (dimage2 && (dimage2->sync == dimage->sync))
        {
          xvGetTitle(i, title, 0) ;
          switch (dimage2->sourceImage->pixel_format)
          {
          case PFDOUBLE:
            sprintf(buf, "%2.2f", 
                    *IMAGEDseq_pix(dimage2->sourceImage, x,y,dimage2->frame));
            break ;
          case PFFLOAT:
            sprintf(buf, "%2.2f",  
                    *IMAGEFseq_pix(dimage2->sourceImage, x,y,dimage2->frame));
            break ;
          case PFBYTE:
            sprintf(buf, "%d", (int)
                    *IMAGEseq_pix(dimage2->sourceImage, x, y,dimage2->frame));
            break ;
          case PFINT:
            sprintf(buf, "%d",
                    *IMAGEIseq_pix(dimage2->sourceImage, x, y,dimage2->frame));
            break ;
          }
          XVshowImageTitle(xvf, i, "%s (%s)", title, buf) ;
        }
      }
    }

    /* extract this location as center of new template */
    if (event_id(event) == LOC_DRAG || event_is_down(event))
    {
    }
    else if (event_is_up(event))
    {
    }
  default:
    if (event_is_up(event))
    {
      switch ((char)event->ie_code)
      {
      case '8':
        if (xvf->get_next_image)
        {
          if (dimage->sync)
          {
            int    which2 ;

            for (which2 = 0 ; which2 < xvf->rows*xvf->cols ; which2++)
            {
              dimage2 = xvGetDimage(xvf, which2, 0) ;
              if (dimage2 && (dimage2->sync == dimage->sync))
              {
                dimage2->sourceImage = 
                  (*xvf->get_next_image)(dimage2->sourceImage, which2, 1) ;
                XVshowImage(xvf, which2, dimage2->sourceImage, dimage2->frame);
              }
            }
          }
          else
          {
            dimage->sourceImage = 
              (*xvf->get_next_image)(dimage->sourceImage, which, 1) ;
            XVshowImage(xvf, which, dimage->sourceImage, dimage->frame) ;
          }
        }
        break ;
      case '2':
        if (xvf->get_next_image)
        {
          if (dimage->sync)
          {
            int    which2 ;

            for (which2 = 0 ; which2 < xvf->rows*xvf->cols ; which2++)
            {
              dimage2 = xvGetDimage(xvf, which2, 0) ;
              if (dimage2 && (dimage2->sync == dimage->sync))
              {
                dimage2->sourceImage = 
                  (*xvf->get_next_image)(dimage2->sourceImage, which2, -1) ;
                XVshowImage(xvf, which2, dimage2->sourceImage, dimage2->frame);
              }
            }
          }
          else
          {
            dimage->sourceImage = 
              (*xvf->get_next_image)(dimage->sourceImage, which, -1) ;
            XVshowImage(xvf, which, dimage->sourceImage, dimage->frame) ;
          }
        }
        break ;
      case 'r':
        xvf->rescale = 0 ;
        XVshowAll(xvf) ;
        break ;
      case 'R':
        xvf->rescale = 1 ;
        XVshowAll(xvf) ;
        break ;
      case 'S':
        XVsyncAll(xvf, which) ;
        break ;
      case 's':
        XVunsyncAll(xvf, which) ;
        break ;
      case 'm':
        scale = 0.9f ;
      case 'M':
        if (event->ie_code == 'M')
          scale = 1.1f ;
        rows = nint((float)dimage->dispImage->rows*scale) ;
        cols = nint((float)dimage->dispImage->cols*scale) ;
        XVsetImageSize(xvf, which, rows, cols) ;
#if 1
        if (rows > xvf->display_rows || cols > xvf->display_cols)
#endif
          XVresize(xvf) ;
        XVshowAllSyncedImages(xvf, which) ;
        break;
      case '*':
        XVgamma(xvf, which, 0.02f) ;
        break ;
      case '/':
        XVgamma(xvf, which, -0.02f) ;
        break ;
      case '+':
        XVbrighten(xvf, which, 1) ;
        break ;
      case '-':
        XVbrighten(xvf, which, -1) ;
        break ;
      case '\n':
      case '\r':
      {
        int i ;
        for (i = 0 ; i < MAX_COLORS ; i++)
          dimage->gamma[i] = (float)i ;
        
        XVbrighten(xvf, which, dimage->bshift = 0) ;
        break ;
      }
      }
      break ;
    case 'D':
    case 'd':
      debug_it = 1 ;
      XVshowImage(xvf, which, dimage->sourceImage, dimage->frame) ;
      break ;
    }
  }
  if (XVevent_handler && !event_ctrl_is_down(event))
  {
    event_x(event) = x ;
    event_y(event) = yprint ;
    (*XVevent_handler)(event, dimage) ;
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVsetParms(void (*event_handler)(Event *event, DIMAGE *dimage))
{
  XVevent_handler = event_handler ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVsetRepaintHandler(void (*repaint_handler)(XV_FRAME *xvf, DIMAGE *dimage))
{
  XVrepaint_handler = repaint_handler ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVsetQuitFunc(void (*quit_func)(void))
{
  XVquit_func = quit_func ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVshowImageTitle(XV_FRAME *xvf, int which, ...)
{
  Panel_item item ;
  va_list args ;
  char    *fmt, *str ;
  int     len, row, col, i, spaces ;
  DIMAGE  *dimage ;

  va_start(args, which) ;
  fmt = va_arg(args, char *) ;

  row = which / xvf->cols ;
  col = which % xvf->cols ;
  if (row >= xvf->rows)
    return(0) ;

  dimage = &xvf->dimages[row][col] ;
  item = dimage->title_item ;
  str = dimage->title_string ;
  vsprintf(str, fmt, args) ;
  spaces = strlen(str) ;
  spaces = MAX_TITLE_CHARS(dimage) - strlen(str) ;
  spaces = nint((float)spaces/1.0f) ;
  for (i = 0 ; i < spaces ; i++)
    str[i] = ' ' ;
  vsprintf(str+i, fmt, args) ;
  StrUpper(str) ;

  len = strlen(str) ;

  if (len > (MAX_TITLE_CHARS(dimage)+spaces/3))
    str[len = MAX_TITLE_CHARS(dimage)] = 0 ;
  va_end(args) ;
  xv_set(item, PANEL_LABEL_STRING, str, NULL);
  XSync(xvf->display, 0) ;
  return(len) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static XImage *
xvCreateXimage(XV_FRAME *xvf, IMAGE *image)
{
  int    rows, cols ;
  XImage *ximage ;

  if (image->pixel_format != PFBYTE)
  {
    fprintf(stderr,
            "xvCreateXimage: unsupported pixel format %d, must be BYTE\n",
            image->pixel_format) ;
    exit(4) ;
  }
  rows = image->rows ;
  cols = image->cols ;
  ximage = 
    XCreateImage(xvf->display, 
                 DefaultVisual(xvf->display, DefaultScreen(xvf->display)),
                 8,                   /* 8 bits per pixel */
                 ZPixmap,             /* not a bitmap */
                 0, image->image,
                 cols,
                 rows,
                 8,                   /* 8 bits per pixel */
                 cols);
  return(ximage) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static Panel_setting
xvHipsCommand(Panel_item item, Event *event)
{
  DIMAGE  *dimage ;

  strcpy(hips_cmd_str, (char *)xv_get(hips_cmd_panel_item, PANEL_VALUE)) ;
  xv_set(hips_cmd_frame, FRAME_CMD_PUSHPIN_IN, FALSE, XV_SHOW, FALSE, NULL) ;

  dimage = xvGetDimage(xvf, hips_cmd_source, 0) ;
  if (!dimage)
    return(0) ;
#if 1
  ImageWrite(dimage->sourceImage, "out.hipl") ;
#else
  if (ImageWriteFrames(dimage->sourceImage, "out.hipl", dimage->frame, 1) < 0)
  {
    XVprintf(xvf, 0, "write failed\n") ;
    return(0) ;
  }
#endif

  if (strlen(hips_cmd_str) < 4)
    return(0) ;

  fprintf(stderr, "executing hips command '%s'\n", hips_cmd_str) ;

  system(hips_cmd_str) ;
  if (strstr(hips_cmd_str, "in.hipl"))
  {
    ImageReadInto("in.hipl", dimage->sourceImage, 0) ;
    XVshowImage(xvf, hips_cmd_source, dimage->sourceImage, dimage->frame) ;
  }

  XFlush(xvf->display); 
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
xvHipsCmdFrameInit(XV_FRAME *xvf) 
{
  xvf->fname_frame = (Frame)
    xv_create(xvf->frame, FRAME_CMD,
              XV_HEIGHT, HIPS_CMD_ROWS,
              XV_WIDTH,  HIPS_CMD_COLS,
              XV_X,      600,
              XV_Y,      310,
              FRAME_LABEL, "FILE NAME",
              XV_SHOW,   FALSE,
              NULL);
  xvf->fname_panel = 
    (Panel)xv_create((Xv_opaque)xvf->fname_frame, PANEL, XV_X, 0, XV_Y,0,NULL);
  xvf->fname_panel_item = (Panel_item)
    xv_create((Xv_opaque)xvf->fname_panel, PANEL_TEXT,
              PANEL_VALUE_STORED_LENGTH,   300,
              PANEL_VALUE_DISPLAY_LENGTH,  40,
              XV_SHOW,                     TRUE,
              PANEL_NOTIFY_PROC,           xvFileNameCommand,
              PANEL_LABEL_STRING,          "File Name: ",
              PANEL_VALUE,                 xvf->file_name,
              NULL) ;

  hips_cmd_frame = (Frame)
    xv_create(xvf->frame, FRAME_CMD,
              XV_HEIGHT, HIPS_CMD_ROWS,
              XV_WIDTH,  HIPS_CMD_COLS,
              XV_X,      600,
              XV_Y,      310,
              FRAME_LABEL, "HIPS COMMAND",
              XV_SHOW,   FALSE,
              NULL);

  hips_cmd_panel = 
    (Panel)xv_create((Xv_opaque)hips_cmd_frame, PANEL, XV_X, 0, XV_Y,0,NULL);
  hips_cmd_display = (Display *)xv_get(hips_cmd_frame, XV_DISPLAY);
  hips_cmd_panel_item = (Panel_item)
    xv_create((Xv_opaque)hips_cmd_panel, PANEL_TEXT,
              PANEL_VALUE_STORED_LENGTH,   300,
              PANEL_VALUE_DISPLAY_LENGTH,  40,
              XV_SHOW,                     TRUE,
              PANEL_NOTIFY_PROC,           xvHipsCommand,
              PANEL_LABEL_STRING,          "Hips Command: ",
              PANEL_VALUE,                 hips_cmd_str,
              NULL) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void 
buttonQuit(Panel_item item, Event *event)
{
  xv_destroy_safe(xvf->frame);
  xv_destroy_safe(hips_cmd_frame) ;
  if (XVquit_func)
    (*XVquit_func)() ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVdrawPoint(XV_FRAME *xvf, int which, int x, int y, int color)
{
  GC      gc ;
  Display *display ;
  Window  window ;
  int     x0, y0, x1, y1 ;
  float   xscale, yscale ;
  DIMAGE  *dimage ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return ;
  
  display = xvf->display ;
  window = dimage->window ;
  switch (color)
  {
  case XWHITE:
    gc = dimage->whiteGC ;
    break ;
  case XBLUE:
    gc = dimage->blueGC ;
    break ;
  case XCYAN:
    gc = dimage->cyanGC ;
    break ;
  case XGREEN:
    gc = dimage->greenGC ;
    break ;
  case XXOR:
    gc = dimage->xorGC ;
    break ;
  case XRED:
  default:
    gc = dimage->redGC ;
    break ;
  }

  xscale = dimage->xscale ;
  yscale = dimage->yscale ;
  x = nint(((float)(x-dimage->x0)+0.5f) * xscale) ;
  y = nint((float)(((dimage->sourceImage->rows-1) - 
                    (y-dimage->y0)) + 0.5f) * yscale) ;
  XSetLineAttributes(display, gc, 0, LineSolid, CapRound, JoinBevel) ;
  
  x0 = x - 4 ;
  y0 = y - 4 ;
  x1 = x + 4 ;
  y1 = y + 4 ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;

  x0 = x + 4 ;
  y0 = y - 4 ;
  x1 = x - 4 ;
  y1 = y + 4 ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;

/*  XFlush(display);*/
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVdrawBox(XV_FRAME *xvf, int which, int x, int y, int dx, int dy, int color)
{
  GC      gc ;
  Display *display ;
  Window  window ;
  int     x0, y0, x1, y1 ;
  DIMAGE  *dimage ;
  float   xscale, yscale ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return ;
  
  display = xvf->display ;
  window = dimage->window ;
  switch (color)
  {
  case XCYAN:
    gc = dimage->cyanGC ;
    break ;
  case XWHITE:
    gc = dimage->whiteGC ;
    break ;
  case XBLUE:
    gc = dimage->blueGC ;
    break ;
  case XGREEN:
    gc = dimage->greenGC ;
    break ;
  case XYELLOW:
    gc = dimage->yellowGC ;
    break ;
  case XPURPLE:
    gc = dimage->purpleGC ;
    break ;
  case XXOR:
    gc = dimage->xorGC ;
    break ;
  case XRED:
  default:
    gc = dimage->redGC ;
    break ;
  }

  /* convert to window coordinate system */
  xscale = dimage->xscale ;
  yscale = dimage->yscale ;

  x = nint((float)(x-dimage->x0) * xscale) ;
  y = nint((float)(((dimage->zoomImage->rows) - (y-dimage->y0))) * yscale) ;
  dx = nint((float)dx * xscale) ;
  dy = nint((float)-dy * yscale) ;

  XSetLineAttributes(display, gc, 0, LineSolid, CapRound, JoinBevel) ;

  /* top line */
  x0 = x ;
  y0 = y ;
  x1 = x + dx ;
  y1 = y0 ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;

  /* bottom line */
  x0 = x ;
  y0 = y + dy ;
  x1 = x + dx ;
  y1 = y + dy ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;

  /* draw left line */
  x0 = x ;
  y0 = y ;
  x1 = x ;
  y1 = y0 + dy ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;

  /* draw right line */
  x0 = x + dx ;
  y0 = y ;
  x1 = x + dx ;
  y1 = y0 + dy ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;

/*  XFlush(display);*/
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVdrawLine(XV_FRAME *xvf, int which, int x, int y, int dx, int dy, int color)
{
  GC      gc ;
  Display *display ;
  Window  window ;
  int     x0, y0, x1, y1 ;
  DIMAGE  *dimage ;
  float   xscale, yscale ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return ;

  display = xvf->display ;
  window = dimage->window ;

  /* convert to window coordinate system */
  xscale = dimage->xscale ;
  yscale = dimage->yscale ;
  x = nint(((float)(x-dimage->x0) +0.5f)* xscale) ;
  y = nint((float)(((dimage->zoomImage->rows-1) - 
                    (y-dimage->y0)) +0.5f)* yscale) ;
  dx = nint((float)dx * xscale) ;
  dy = nint((float)-dy * yscale) ;

  switch (color)
  {
  case XWHITE:
    gc = dimage->whiteGC ;
    break ;
  case XCYAN:
    gc = dimage->cyanGC ;
    break ;
  case XBLUE:
    gc = dimage->blueGC ;
    break ;
  case XGREEN:
    gc = dimage->greenGC ;
    break ;
  case XXOR:
    gc = dimage->xorGC ;
    break ;
  case XRED:
  default:
    gc = dimage->redGC ;
    break ;
  }

  XSetLineAttributes(display, gc, 0, LineSolid, CapRound, JoinBevel) ;

  /* top line */
  x0 = x ;
  y0 = y ;
  x1 = x + dx ;
  y1 = y + dy ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;


/*  XFlush(display);*/
}
#define ARROW_HEAD_SIZE    8.0f
#define ARROW_HEAD_ANGLE   RADIANS(30.0f)

void
XVdrawArrow(XV_FRAME *xvf, int which, int x, int y,float dx,float dy,int color)
{
  GC      gc ;
  Display *display ;
  Window  window ;
  int     x0, y0, x1, y1 ;
  DIMAGE  *dimage ;
  float   xscale, yscale, theta, theta0 ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return ;

  display = xvf->display ;
  window = dimage->window ;

  /* convert to window coordinate system */
  xscale = dimage->xscale ;
  yscale = dimage->yscale ;
  x = nint(((float)(x-dimage->x0)+0.5f) * xscale) ;
  y = nint(((float)((dimage->zoomImage->rows-1) - 
                    (y-dimage->y0)) + 0.5f) * yscale) ;
  dx = nint(dx * xscale) ;
  dy = nint(-dy * yscale) ;

  switch (color)
  {
  case XWHITE:
    gc = dimage->whiteGC ;
    break ;
  case XCYAN:
    gc = dimage->cyanGC ;
    break ;
  case XBLUE:
    gc = dimage->blueGC ;
    break ;
  case XXOR:
    gc = dimage->xorGC ;
    break ;
  case XGREEN:
    gc = dimage->greenGC ;
    break ;
  case XRED:
  default:
    gc = dimage->redGC ;
    break ;
  }

  XSetLineAttributes(display, gc, 0, LineSolid, CapRound, JoinBevel) ;

  /* top line */
  x0 = x ;
  y0 = y ;
  x1 = nint((float)x + dx) ;
  y1 = nint((float)y + dy) ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;

  /* now draw arrow head */
  if (iszero(dx) && iszero(dy))
    theta0 = 0 ;
  else
    theta0 = atan2(dy, dx) ;

  x0 = x1 ;
  y0 = y1 ;

  /* top portion */
  theta = theta0 + PI-ARROW_HEAD_ANGLE ;
  x1 = x0 + nint(ARROW_HEAD_SIZE * cos(theta)) ;
  y1 = y0 + nint(ARROW_HEAD_SIZE * sin(theta)) ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;

  /* bottom portion */
  theta = theta0 + ARROW_HEAD_ANGLE - PI ;
  x1 = x0 + nint(ARROW_HEAD_SIZE * cos(theta)) ;
  y1 = y0 + nint(ARROW_HEAD_SIZE * sin(theta)) ;
  XDrawLine(display, window, gc, x0, y0, x1, y1) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
              xo,yo is image coordinate of center of area, x1,y1 and
              x2,y2 are left and right side of area in image coords.
----------------------------------------------------------------------*/
void
XVshowVectorImage(XV_FRAME *xvf, int which, int xo, int yo, 
                  int width, int height, int color, IMAGE *image)
{
  int     x, y, rows, cols, x0, y0, x1, y1, half_width, half_height ;
  float   dx, dy, len ;

  rows = image->rows ;
  cols = image->cols ;
  half_width = (width - 1) / 2 ;
  half_height = (height - 1) / 2 ;
  x1 = xo + half_width ;
  y1 = yo + half_height ;
  if (x1 >= cols)
    x1 = cols - 1 ;
  if (y1 >= rows)
    y1 = rows - 1 ;
  x0 = xo - half_width ;
  y0 = yo - half_height ;
  if (x0 < 0) 
    x0 = 0 ;
  if (y0 < 0)
    y0 = 0 ;
  
  for (y = y0 ; y <= y1 ; y++)
  {
    for (x = x0 ; x <= x1 ; x++)
    {
      dx = *IMAGEFseq_pix(image, x, y, 0) ;
      dy = *IMAGEFseq_pix(image, x, y, 1) ;
      len = 2*sqrt (dx * dx + dy * dy) ;
      if (iszero(dx) && iszero(dy))
        XVdrawPoint(xvf, which, x-x0, y-y0, color) ;
      else
        XVdrawArrow(xvf, which, x-x0, y-y0, dx/len, dy/len, color) ;
    }
  }
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
void
XVsetKBhandler(void (*kb_handler)(Event *event, DIMAGE *dimage))
{
  XVkb_handler = kb_handler ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVsetMessagePosition(XV_FRAME *xvf, int which, int col, int row)
{
  if (which < 0 || which > 1)
    ErrorReturn(-1, (ERROR_BADPARM, "XVsetMessagePosition(%d): out of range", 
                    which)) ;

  xv_set((Xv_opaque)xvf->msg_item[which], XV_X, col, XV_Y, row, NULL);
  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVsetImageSize(XV_FRAME *xvf, int which, int rows, int cols)
{
  DIMAGE    *dimage, *dimage2 = NULL ;
  char      title[100] ;
  int       i, done = 0 ;

  dimage = xvGetDimage(xvf, which, 1) ;
  if (!dimage || rows <= 0 || cols <= 0)
    return(-1) ;

  i = which ;
  do
  {
    xv_set(dimage->canvas, CANVAS_WIDTH, cols, CANVAS_HEIGHT, rows, NULL) ;
    xv_set(dimage->canvas, XV_WIDTH, cols, XV_HEIGHT, rows, NULL) ;
    
    xv_set(canvas_paint_window(dimage->canvas),
           WIN_EVENT_PROC, xv_dimage_event_handler, WIN_CONSUME_EVENTS,  
           MS_LEFT, LOC_DRAG, MS_RIGHT, MS_MIDDLE, WIN_ASCII_EVENTS, 
           WIN_META_EVENTS, WIN_LEFT_KEYS, WIN_RIGHT_KEYS,NULL, NULL);
    dimage->window = 
      (Window)xv_get(canvas_paint_window(dimage->canvas),XV_XID);
    
    ImageFree(&dimage->dispImage) ;
    dimage->dispImage = ImageAlloc(rows, cols, PFBYTE, 1) ;
    
    /* should free the ximage and the canvas, but don't know how yet */
    dimage->ximage = xvCreateXimage(xvf, dimage->dispImage) ;
    if (xvGetTitle(i, title, 1))
      XVshowImageTitle(xvf, i, title) ;

/* 
   must show image here to refill dispImage, otherwise repaint will display
   a blank image.
*/
    XVshowImage(xvf, i, dimage->sourceImage, dimage->frame) ;

    if (dimage->sync)
    {
      if (i == which)    /* just did one passed in by caller */
        i = -1 ;

      for (i++ ; i < xvf->rows*xvf->cols ; i++)
      {
        if (i == which)
          continue ;
        dimage2 = xvGetDimage(xvf, i, 0) ;
        if (dimage2 &&  (dimage2->sync == dimage->sync))
          break ;
      }
      if (i < xvf->rows*xvf->cols)
      {
        done = 0 ;
        dimage = dimage2 ;
      }
      else
        done = 1 ;
    }
    else
      done = 1 ;
  } while (!done) ;

  return(0) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
int
XVresize(XV_FRAME *xvf)
{
  DIMAGE     *dimage ;
  int        row, col, x, y, max_height, max_width ;

  y = PANEL_HEIGHT ;
  max_width = 0 ;
  dimage = &xvf->dimages[0][0] ;
  xvf->display_rows = dimage->dispImage->rows ;
  xvf->display_cols = dimage->dispImage->cols ;
  for (row = 0 ; row < xvf->rows ; row++)
  {
    x = 0 ;
    max_height = 0 ;
    for (col = 0 ; col < xvf->cols ; col++)
    {
      dimage = &xvf->dimages[row][col] ;
      dimage->which = row * xvf->cols + col ;
      dimage->x = x ;
      dimage->y = y ;

      xv_set(dimage->canvas,
             XV_X,                  x,
             XV_Y,                  y,
             NULL) ;

      xv_set(dimage->title_item,
             XV_X, x,
             XV_Y, y-CHAR_HEIGHT+CHAR_PAD,
             NULL) ;


      x += WINDOW_PAD + dimage->dispImage->cols ;
      if (dimage->dispImage->rows >= max_height)
        max_height = dimage->dispImage->rows ;
    }
    if (x >= max_width)
      max_width = x ;
    y += CHAR_HEIGHT + max_height ;
  }
  y -= CHAR_HEIGHT ;  /* don't need room for characters at bottom */
  if (xvf->min_panel_width && xvf->min_panel_width > max_width)
    max_width = xvf->min_panel_width ;
  xv_set(xvf->frame, XV_WIDTH, max_width, XV_HEIGHT, y, NULL) ;
  return(0) ;
}

void
XVsetMinPanelWidth(XV_FRAME *xvf, int min_panel_width)
{
  xvf->min_panel_width = min_panel_width ;
  xv_set(xvf->frame, XV_WIDTH, min_panel_width, NULL) ;
}

int
XVaddImageCol(XV_FRAME *xvf)
{
  int     row, col, which, x, y ;
  DIMAGE  *dimages ;

  col = xvf->cols++ ;

  x = col * (WINDOW_PAD + xvf->display_cols) ;
  y = PANEL_HEIGHT ;
  for (row = 0 ; row < xvf->rows ; row++)
  {
    dimages = (DIMAGE *)calloc(xvf->cols, sizeof(DIMAGE)) ;
    for (col = 0 ; col < xvf->cols ; col++)
    {
      if (col < xvf->cols-1)
        memcpy(&dimages[col], &xvf->dimages[row][col], sizeof(DIMAGE)) ;
      else
      {
        which = row * xvf->cols + col ;
        xvCreateImage(xvf, &dimages[col], x, y, which) ;
      }
    }
    free(xvf->dimages[row]) ;
    xvf->dimages[row] = dimages ;
    y += CHAR_HEIGHT + xvf->display_rows ;
  }

  return(0) ;
}

int
XVdeleteImageCol(XV_FRAME *xvf)
{
  int     row, col ;
  DIMAGE  *dimages ;

  xvf->cols-- ;
  
  for (row = 0 ; row < xvf->rows ; row++)
  {
    dimages = (DIMAGE *)calloc(xvf->cols, sizeof(DIMAGE)) ;
    for (col = 0 ; col < xvf->cols ; col++)
      memcpy(&dimages[col], &xvf->dimages[row][col], sizeof(DIMAGE)) ;

    xvFreeDimage(&xvf->dimages[row][xvf->cols]) ;
    free(xvf->dimages[row]) ;
    xvf->dimages[row] = dimages ;
  }

  XVresize(xvf) ;

  return(0) ;
}

static void
xvCreateImage(XV_FRAME *xvf, DIMAGE *dimage, int x, int y, int which)
{
  XGCValues   GCvalues ;
  int         i ;

  for (i = 0 ; i < MAX_COLORS ; i++)
    dimage->gamma[i] = (float)i ;

  strcpy(dimage->title_string, "") ;
  dimage->zoom = 1.0f ;
  dimage->which = which ;
  dimage->x = x ;
  dimage->y = y ;
  dimage->bshift = 0 ;
  
  dimage->canvas =
    (Canvas)xv_create((Xv_opaque)xvf->frame, CANVAS,
                      XV_X,                  x,
                      XV_Y,                  y,
                      XV_HEIGHT,             xvf->display_rows,
                      XV_WIDTH,              xvf->display_cols,
                      CANVAS_X_PAINT_WINDOW, TRUE,
                      CANVAS_REPAINT_PROC,   xv_dimage_repaint,
                      CANVAS_RETAINED,       FALSE,
                      WIN_CMS,               xvf->cms,
                      NULL);
  
  dimage->title_item = (Panel_item)
    xv_create((Xv_opaque)xvf->panel, PANEL_MESSAGE,
              PANEL_LABEL_BOLD, TRUE,
              XV_X, x,
              XV_Y, y-CHAR_HEIGHT+CHAR_PAD,
              PANEL_LABEL_STRING, 
              dimage->title_string,
              NULL);
  
  dimage->dispImage = ImageAlloc(xvf->display_rows, xvf->display_cols, 
                                 PFBYTE, 1) ;
  dimage->ximage = xvCreateXimage(xvf, dimage->dispImage) ;
  xv_set(canvas_paint_window(dimage->canvas),
         WIN_EVENT_PROC, xv_dimage_event_handler,
         WIN_CONSUME_EVENTS,  MS_LEFT, LOC_DRAG, MS_RIGHT, MS_MIDDLE, 
         WIN_ASCII_EVENTS, NULL,
         NULL);
  dimage->window = 
    (Window)xv_get(canvas_paint_window(dimage->canvas),XV_XID);
  
  dimage->clearGC = 
    XCreateGC(xvf->display, dimage->window, (unsigned long )0, &GCvalues);
  XSetFunction(xvf->display, dimage->clearGC, GXclear);
  dimage->xorGC = 
    XCreateGC(xvf->display, dimage->window, (unsigned long)0, &GCvalues);
  XSetFunction(xvf->display, dimage->xorGC, GXinvert);
  
  dimage->greenGC = 
    XCreateGC(xvf->display, dimage->window,(unsigned long )0, &GCvalues);
  XSetFunction(xvf->display, dimage->greenGC, GXcopy) ;
  XSetBackground(xvf->display, dimage->greenGC, xvf->green_pixel) ;
  XSetForeground(xvf->display, dimage->greenGC, xvf->green_pixel) ;
  
  dimage->blueGC = 
    XCreateGC(xvf->display, dimage->window,(unsigned long )0, &GCvalues);
  XSetFunction(xvf->display, dimage->blueGC, GXcopy) ;
  XSetBackground(xvf->display, dimage->blueGC, xvf->blue_pixel) ;
  XSetForeground(xvf->display, dimage->blueGC, xvf->blue_pixel) ;
  
  dimage->purpleGC = 
    XCreateGC(xvf->display, dimage->window,(unsigned long )0, &GCvalues);
  XSetFunction(xvf->display, dimage->purpleGC, GXcopy) ;
  XSetBackground(xvf->display, dimage->purpleGC, xvf->purple_pixel) ;
  XSetForeground(xvf->display, dimage->purpleGC, xvf->purple_pixel) ;
  
  dimage->yellowGC = 
    XCreateGC(xvf->display, dimage->window,(unsigned long )0, &GCvalues);
  XSetFunction(xvf->display, dimage->yellowGC, GXcopy) ;
  XSetBackground(xvf->display, dimage->yellowGC, xvf->yellow_pixel) ;
  XSetForeground(xvf->display, dimage->yellowGC, xvf->yellow_pixel) ;
  
  dimage->cyanGC = 
    XCreateGC(xvf->display, dimage->window,(unsigned long )0, &GCvalues);
  XSetFunction(xvf->display, dimage->cyanGC, GXcopy) ;
  XSetBackground(xvf->display, dimage->cyanGC, xvf->cyan_pixel) ;
  XSetForeground(xvf->display, dimage->cyanGC, xvf->cyan_pixel) ;
  
  dimage->redGC = 
    XCreateGC(xvf->display, dimage->window,(unsigned long )0, &GCvalues);
  XSetFunction(xvf->display, dimage->redGC, GXcopy) ;
  XSetBackground(xvf->display, dimage->redGC, xvf->red_pixel) ;
  XSetForeground(xvf->display, dimage->redGC, xvf->red_pixel) ;
  
  dimage->whiteGC = 
    XCreateGC(xvf->display, dimage->window,(unsigned long )0, &GCvalues);
  XSetForeground(xvf->display, dimage->whiteGC, xvf->white_pixel);
  XSetBackground(xvf->display, dimage->whiteGC, xvf->white_pixel);
}

static void
xvFreeDimage(DIMAGE *dimage)
{
  ImageFree(&dimage->dispImage) ;
  xv_destroy_safe(dimage->canvas) ;
  xv_destroy_safe(dimage->title_item) ;

/*  (*dimage->ximage->destroy_image)(dimage->ximage)*/
  /* lots of other stuff needed here */
}

int
XVsetPrecision(XV_FRAME *xvf, int precision)
{
  xvf->precision = precision ;
  return(0) ;
}

int
XVgamma(XV_FRAME *xvf, int which, float beta)
{
  float         gamma ;
  int           i ;
  DIMAGE        *dimage ;
  double        dval ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return(NO_ERROR) ;

  if (beta < -0.9999f)
    beta = -0.9999f ;
  else if (beta > 1.0f)
    beta = 1.0f ;


  if (beta > 0)
    gamma = 1 - beta ;
  else
    gamma = 1 / (1+beta) ;

  for (i=0;i<MAX_COLORS;i++)
  {
    dval = dimage->gamma[i] / (float)MAX_GRAY ;
    dval = pow(dval, gamma) ;
    dimage->gamma[i] = dval * (float)MAX_GRAY ;
  }

  XVshowImage(xvf, which, dimage->sourceImage, dimage->frame) ;

  return(NO_ERROR) ;
}

int
XVzoom(XV_FRAME *xvf, int which, float zoom)
{
  float   max_dim ;
  DIMAGE  *dimage ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return(NO_ERROR) ;

  max_dim = (float)
    MIN(dimage->sourceImage->cols, dimage->sourceImage->rows) ;
  dimage->zoom *= zoom ;
  if (dimage->zoom > max_dim)
    dimage->zoom = max_dim ;
  if (dimage->zoom < 1.0f)
    dimage->zoom = 1.0f ;

/*  fprintf(stderr, "zoom = %2.3f\n", dimage->zoom) ;*/

  XVshowImage(xvf, which, dimage->sourceImage, dimage->frame) ;
  return(NO_ERROR) ;
}

int
XVsync(XV_FRAME *xvf, int which, int sync)
{
  DIMAGE  *dimage ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return(ERROR_BAD_PARM) ;

  dimage->sync = sync ;
  return(NO_ERROR) ;
}
int
XVdoSync(XV_FRAME *xvf, int which)
{
  DIMAGE  *dimage, *dimage2 ;
  int     i, rows, cols ;
  float   aspect ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return(ERROR_BAD_PARM) ;

  for (i = 0 ; i < xvf->rows*xvf->cols ; i++)
  {
    if (i == which)
      continue ;
    dimage2 = xvGetDimage(xvf, i, 0) ;
    if (dimage2 && (dimage2->sync == dimage->sync))
    {
      dimage2->dx = dimage->dx ;
      dimage2->dy = dimage->dy ;
      dimage2->x0 = dimage->x0 ;
      dimage2->y0 = dimage->y0 ;
      dimage2->x1 = dimage->x1 ;
      dimage2->y1 = dimage->y1 ;
      dimage2->dx1 = dimage->dx1 ;
      dimage2->dy1 = dimage->dy1 ;
      rows = xvf->display_rows ;
      cols = xvf->display_cols ;
      if (dimage->dx)
        aspect = (float)dimage->dy / (float)dimage->dx ;
      else
        aspect = 
          (float)dimage->sourceImage->rows / (float)dimage->sourceImage->cols ;

#if 0
      if (aspect > 1.0f)
        XVsetImageSize(xvf, i, rows, nint((float)rows / aspect)) ;
      else
        XVsetImageSize(xvf, i, nint((float)cols * aspect),cols) ;
#endif
/*      XVshowImage(xvf, i, dimage2->sourceImage, dimage->frame) ;*/
    }
  }
  return(NO_ERROR) ;
}

int
XVsyncAll(XV_FRAME *xvf, int which)
{
  int    rows, cols, which2, sync ;
  DIMAGE *dimage, *dimage2 ;

#if 0
  for (sync = which = 0 ; which < xvf->rows*xvf->cols ; which++)
#else
    sync = 0 ;
#endif
  {
    dimage = xvGetDimage(xvf, which, 0) ;
    if (dimage && (dimage->sync == 0))
    {
      dimage->sync = ++sync ;
      rows = dimage->sourceImage->rows ;
      cols = dimage->sourceImage->cols ;
      for (which2 = 0 ; which2 < xvf->rows*xvf->cols ; which2++)
      {
        dimage2 = xvGetDimage(xvf, which2, 0) ;
        if (dimage2 && (dimage2->sourceImage->rows == rows) && 
            (dimage2->sourceImage->cols == cols))
          dimage2->sync = sync ;
      }
      XVdoSync(xvf, which) ;
    }
  }
  XVshowAllSyncedImages(xvf, which) ;
  return(NO_ERROR) ;
}
int
XVshowAll(XV_FRAME *xvf)
{
  DIMAGE *dimage ;
  int    which ;

  for (which = 0 ; which < xvf->rows*xvf->cols ; which++)
  {
    dimage = xvGetDimage(xvf, which, 0) ;
    if (dimage)
      XVshowImage(xvf, which, dimage->sourceImage, dimage->frame) ;
  }
  
  return(NO_ERROR) ;
}
int
XVunsyncAll(XV_FRAME *xvf, int which)
{
  DIMAGE *dimage ;

  for (which = 0 ; which < xvf->rows*xvf->cols ; which++)
  {
    dimage = xvGetDimage(xvf, which, 0) ;
    if (dimage && dimage->sync)
    {
      XVdoSync(xvf, which) ;
      dimage->sync = 0 ;
    }
  }
  return(NO_ERROR) ;
}
int
XVshowAllSyncedImages(XV_FRAME *xvf, int which)
{
  int     which2 ;
  DIMAGE  *dimage, *dimage2 ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return(ERROR_BADPARM) ;
  XVshowImage(xvf, which, dimage->sourceImage, dimage->frame) ;

  for (which2 = 0 ; which2 < xvf->rows*xvf->cols ; which2++)
  {
    dimage2 = xvGetDimage(xvf, which2, 0) ;
    if (dimage2 && (dimage2->sync == dimage->sync))
      XVshowImage(xvf, which2, dimage2->sourceImage, dimage2->frame) ;
  }
  return(NO_ERROR) ;
}

int
XVsetPrintStatus(XV_FRAME *xvf, int status)
{
  if (!status)
    xvf->noprint = 1 ;
  else
    xvf->noprint = 0 ;
  return(NO_ERROR) ;
}
int
XVsetYDir(XV_FRAME *xvf, int ydir)
{
  xvf->ydir = ydir ;
  return(NO_ERROR) ;
}

int 
XVgetFileName(XV_FRAME *xvf, char *default_fname,
              int (*fname_func)(char *fname), ...)
{
  va_list  args ;
  char     *fmt ;
  int      len ;

  va_start(args, fname_func) ;
  fmt = va_arg(args, char *) ;
  vsprintf(xvf->fname_prompt, fmt, args) ;
  len = strlen(xvf->fname_prompt) ;

  va_end(args) ;

  xvf->fname_func = fname_func ;
  strcpy(xvf->file_name, default_fname) ;
  xv_set(xvf->fname_frame, FRAME_LABEL, xvf->fname_prompt, PANEL_VALUE, 
         xvf->file_name, FRAME_CMD_PUSHPIN_IN, TRUE, XV_SHOW, TRUE, NULL) ;
  return(len) ;
}

static Panel_setting
xvFileNameCommand(Panel_item item, Event *event)
{
  strcpy(xvf->file_name, (char *)xv_get(xvf->fname_panel_item, PANEL_VALUE)) ;
  xv_set(xvf->fname_frame, FRAME_CMD_PUSHPIN_IN, FALSE, XV_SHOW, FALSE, NULL) ;
  XFlush(xvf->display); 

  if (strlen(xvf->file_name) < 4)
    return(0) ;

  if (xvf->fname_func)
    (*xvf->fname_func)(xvf->file_name) ;
  return(0) ;
}

int
XVsetDepthFunc(XV_FRAME *xvf,
               IMAGE *(*get_image)(IMAGE *Iold,int which, int dir))
{
  xvf->get_next_image = get_image ;
  return(NO_ERROR) ;
}

static char *
xvGetTitle(int which, char *title, int with_value)
{
  char   *cp ;
  DIMAGE *dimage ;

  dimage = xvGetDimage(xvf, which, 0) ;
  if (!dimage)
    return(NULL) ;
  strcpy(title, StrRemoveSpaces(dimage->title_string)) ;
  if (!with_value)  /* remove embedded value in title string: <title> (val) */
  {
    cp = strrchr(title, '(') ;
    if (cp)
      *cp-- = 0 ;
  }
  else
    cp = title+strlen(title) ;

  if (cp)   /* strip trailing spaces */
    while (isspace(*cp))
      *cp-- = 0 ;

  return(title) ;
}

