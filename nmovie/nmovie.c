#include <stdio.h>
#include <image.h>
#include <stdlib.h>
#include "error.h"

#include <X11/Xlib.h>
#include <X11/Intrinsic.h>
#include <X11/StringDefs.h>
#include <X11/Shell.h>

#include <X11/Xaw/Form.h>
#include <X11/Xaw/Command.h>
#include <X11/Xaw/Simple.h>

#include <sys/ipc.h>
#include <sys/shm.h>
#include <X11/extensions/XShm.h>

#include <X11/extensions/xf86dga.h>

/* -------- Prototypes -------- */
void useage(void);
void quit_proc(Widget w, XEvent *event, String *pars, Cardinal *npars);
void repaint_proc(Widget w, XEvent *event, String *pars, Cardinal *npars);
static void XInit(int *argc, char ***argv);
static int highbit(unsigned long ul);
static void XSetupDisplay(int nframes);
/* -------- End Prototypes -------- */

/* -------- Typedefs -------- */
typedef struct {
  XtAppContext context;
  Display *disp;
  GC theGC;
  XVisualInfo vi;
  Visual *vis;
  Window root;
  Window canvas;
  Colormap colormap;
  int screenno;
  int depth;
  unsigned long rmask, gmask, bmask;
  int rshift, gshift, bshift;
} XCRUFT;
/* -------- End Typedefs -------- */
  
/* -------- Global Variables ------- */
Widget frame,toplevel,quit_bt,canvas,buttons,stop_bt;
int shmext,dgaext,nframes;
XImage *ximg;
byte *imgdata;
XCRUFT xi;
int pfmt=0,rows=0,cols=0;
/* -------- End Global Variables ------- */

/* -------- Static Variables -------- */
static XtActionsRec actions[] = {
  {"Refresh", repaint_proc},
  {"Quit", quit_proc}
};

static char *fallback_resources[] = {
  "*Shell.visual:          TrueColor",
  "*Background:            black",
  "*BorderColor:           white",
  "*Foreground:            wheat",
  "*Frame*sensitive:       True",
  "*Font:                  -*-helvetica-bold-r-*-*-*-120-*-*-*-*-*-*",
  "*Buttons.borderWidth:   0",
  "*Buttons.hSpace:        0",
  "*Canvas.fromVert:       Buttons",
  "*Quit.Translations:     #replace <Btn1Down>,<Btn1Up>:Quit() notify()",
  "*Quit.sensitive:        True",
  "*Quit.fromHoriz:        Stop",
  "*Canvas.Translations:   #override <Expose>:Refresh()",
  NULL
};
/* -------- End Static Variables -------- */

void useage(void)
{
  fprintf(stderr,"nmovie <image file> <image file> ...\n");
  exit(1);
}

void repaint_proc(Widget w, XEvent *event, String *pars, Cardinal *npars)
{
  if (w == canvas)
    {
      XPutImage(xi.disp, xi.canvas, xi.theGC, ximg, 0, 0, 0, 0, cols, rows);
      fprintf(stderr,"Expose event for canvas\n");
    }
}

void quit_proc(Widget w, XEvent *event, String *pars, Cardinal *npars)
{
  XtDestroyWidget(toplevel);
  exit(0);
}

static int highbit(unsigned long ul)
{
  /* returns position of highest set bit in 'ul' as an integer (0-31),
   or -1 if none */

  int i;  unsigned long hb;

  hb = 0x80;  hb = hb << 24;   /* hb = 0x80000000UL */
  for (i=31; ((ul & hb) == 0) && i>=0;  i--, ul<<=1);
  return i;
}

static void XInit(int *argc, char ***argv)
{
  int ebase,errbase,flags;

  XtToolkitInitialize();
  xi.context = XtCreateApplicationContext();
  XtAppSetFallbackResources(xi.context, fallback_resources);
  xi.disp = XtOpenDisplay(xi.context, NULL, "NMovie", "NMovie", NULL, 0, 
			  argc, *argv);

  shmext = XShmQueryExtension(xi.disp);
  xi.screenno = DefaultScreen(xi.disp);

  if (XF86DGAQueryExtension(xi.disp,&ebase,&errbase))
    {
      XF86DGAQueryDirectVideo(xi.disp, xi.screenno, &flags);
      dgaext = flags & XF86DGADirectPresent;
    }
  else
    dgaext = 0;

}

static void XSetupDisplay(int nframes)
{
  XGCValues xgcv;

  xi.depth = DefaultDepthOfScreen(DefaultScreenOfDisplay(xi.disp));

  if (!XMatchVisualInfo(xi.disp, xi.screenno, xi.depth, TrueColor, &(xi.vi)))
    ErrorExit(ERROR_BADPARM, "Could not find a TrueColor visual");

  xi.vis = xi.vi.visual;
  xi.root = RootWindow(xi.disp, xi.screenno);

  xi.colormap = XCreateColormap(xi.disp, xi.root, xi.vis, AllocNone);

  toplevel = XtVaAppCreateShell("NMovie", "NMovie", 
				applicationShellWidgetClass,
				xi.disp, 
				XtNvisual, xi.vis, 
				XtNcolormap, xi.colormap, 
				NULL);

  XtAppAddActions(xi.context,actions,XtNumber(actions));

  frame = XtVaCreateManagedWidget("Frame", formWidgetClass, toplevel, NULL);
  buttons = XtVaCreateManagedWidget("Buttons", formWidgetClass, frame, NULL );
  stop_bt = XtVaCreateManagedWidget("Stop", commandWidgetClass, 
				    buttons, NULL); 
  quit_bt = XtVaCreateManagedWidget("Quit", commandWidgetClass, 
				    buttons, NULL);
  canvas = XtVaCreateManagedWidget("Canvas", simpleWidgetClass, frame, 
				   XtNwidth, cols,
				   XtNheight, rows, NULL);
  XtInstallAllAccelerators(frame,toplevel);
  XtRealizeWidget(toplevel);
  xi.canvas = XtWindow(canvas);

  xi.theGC = XCreateGC(xi.disp, xi.canvas, 0L, &xgcv);

  xi.rmask = xi.vis->red_mask;
  xi.gmask = xi.vis->green_mask;
  xi.bmask = xi.vis->blue_mask;

  xi.rshift = 7 - highbit(xi.rmask);
  xi.gshift = 7 - highbit(xi.gmask);
  xi.bshift = 7 - highbit(xi.bmask);

  ximg = XCreateImage(xi.disp,xi.vis,xi.depth,ZPixmap, 0, NULL, 
		      cols, rows, 32, 0);

  if ((imgdata = (byte *)calloc((size_t)(rows*ximg->bytes_per_line*nframes),
				sizeof(byte))) 
      ==NULL)
    ErrorExit(ERROR_NO_MEMORY,"Failed to allocate image buffer");

  ximg->data = (char *) imgdata;
}

void rgb2xcol(IMAGE *I, byte *ximgdata, int fnum)
{
  int i,j;
  unsigned long r,g,b,xcol;
  byte *bptr, *xptr, *ip;

  bptr = I->image; 
  xptr = ximgdata+fnum*rows*ximg->bytes_per_line;
  for(i=0;i<rows;i++,xptr += ximg->bytes_per_line)
    for(j=0, ip = xptr; j<cols; j++)
      {
	r = *bptr++; g = *bptr++; b = *bptr++;
	
	if (xi.rshift<0) 
	  r = r << (-xi.rshift);
	else 
	  r = r >> xi.rshift;
	if (xi.gshift<0) 
	  g = g << (-xi.gshift);
	else 
	  g = g >> xi.gshift;
	if (xi.bshift<0) 
	  b = b << (-xi.bshift);
	else 
	  b = b >> xi.bshift;

	r = r & xi.rmask;
	g = g & xi.gmask;
	b = b & xi.bmask; 

	xcol = r | g | b;
	*ip++ = xcol & 0xff;
      }
}

void ConvertImages(int nframes, char **argv)
{
  IMAGE *I;
  int i;
  
  for(i=0;i<nframes;i++)
    {
      I = ImageRead(argv[i+1]);
      rgb2xcol(I,imgdata,i);
    }
}

int main(int argc, char **argv)
{
  IMAGE *I;
  int i;

  XInit(&argc,&argv);

  if (argc<2)
    useage();

  for(i=1;i<argc;i++)
    {
      I = ImageReadHeader(argv[i]);
      pfmt = MAX(pfmt,I->pixel_format);
      rows = MAX(rows,I->orows);
      cols = MAX(cols,I->ocols);
    }

  ImageFree(&I);

  nframes = argc-1;

  XSetupDisplay(nframes);

  ConvertImages(nframes,argv);

  XtMapWidget(toplevel);

  XtAppMainLoop(xi.context);

  return 0; /* Make -Wall happy! */
}



