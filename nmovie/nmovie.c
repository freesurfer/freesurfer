#include <stdio.h>
#include <image.h>
#include <stdlib.h>

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
static int highbit(unsigned long ul);
/* -------- End Prototypes -------- */

/* -------- Global Variables ------- */
XtAppContext context;
Display *disp;
Widget frame,toplevel,quit_bt,canvas,buttons,stop_bt;
int shmext,dgaext;
GC theGC;
IMAGE *I;
XImage *timg;
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
  exit(1);
}

void repaint_proc(Widget w, XEvent *event, String *pars, Cardinal *npars)
{
  if (w == canvas)
    {
      XPutImage(disp, XtWindow(canvas), theGC, timg, 0, 0, 0, 0, I->ocols, 
		I->orows);
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

int main(int argc, char **argv)
{
  int ebase,errbase,flags,pfmt=0,rows=0,cols=0,i,cnt,screenno;
  IMAGE *I2;
  XVisualInfo vi;
  Colormap colormap;
  XGCValues xgcv;
  unsigned long r,g,b,rmask,gmask,bmask;
  int rshift, gshift, bshift, bperpix, bperline, cshift, maplen, xcol;
  byte *imgdata, *bptr, *xptr, *ip;
  int j;

  XtToolkitInitialize();
  context = XtCreateApplicationContext();
  XtAppSetFallbackResources(context, fallback_resources);
  disp = XtOpenDisplay(context, NULL, "NMovie", "NMovie", NULL, 0, 
		       &argc, argv);

  shmext = XShmQueryExtension(disp);

  if (XF86DGAQueryExtension(disp,&ebase,&errbase))
    {
      XF86DGAQueryDirectVideo(disp, DefaultScreen(disp), &flags);
      dgaext = flags & XF86DGADirectPresent;
    }
  else
    dgaext = 0;

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

  I = ImageAlloc(rows, cols, pfmt, argc-1);

  for (i=1;i<argc;i++)
    {
      I2 = ImageRead(argv[i]);
      ImageCopyFrames(I2,I,0,1,i-1);
    }

  ImageFree(&I2);

  screenno = DefaultScreen(disp);
  if (XMatchVisualInfo(disp, screenno, 
		       DefaultDepthOfScreen(DefaultScreenOfDisplay(disp)), 
		       TrueColor, &vi))
    fprintf(stderr,"Found a TrueColor visual for the current screen\n");

  colormap = XCreateColormap(disp, RootWindow(disp, screenno), vi.visual,
			     AllocNone);

  toplevel = XtVaAppCreateShell("NMovie", "NMovie", 
				applicationShellWidgetClass,
				disp, 
				XtNvisual, vi.visual, 
				XtNcolormap,
				colormap, NULL);

  XtAppAddActions(context,actions,XtNumber(actions));

  frame = XtVaCreateManagedWidget("Frame", formWidgetClass, toplevel, NULL);
  buttons = XtVaCreateManagedWidget("Buttons", formWidgetClass, frame, NULL );

  stop_bt = XtVaCreateManagedWidget("Stop", commandWidgetClass, 
				    buttons, NULL); 

  quit_bt = XtVaCreateManagedWidget("Quit", commandWidgetClass, 
				    buttons, NULL);

  canvas = XtVaCreateManagedWidget("Canvas", simpleWidgetClass, frame, 
				   XtNwidth, I->ocols,
				   XtNheight, I->orows, NULL);


  XtInstallAllAccelerators(frame,toplevel);
  XtRealizeWidget(toplevel);
  theGC = XCreateGC(disp, XtWindow(canvas), 0L, &xgcv);

  rmask = vi.visual->red_mask;
  gmask = vi.visual->green_mask;
  bmask = vi.visual->blue_mask;

  rshift = 7 - highbit(rmask);
  gshift = 7 - highbit(gmask);
  bshift = 7 - highbit(bmask);

  maplen = vi.visual->map_entries;
  if (maplen>256) maplen = 256;
  cshift = 7 - highbit((unsigned long)(maplen-1));

  timg = XCreateImage(disp,vi.visual,
		      DefaultDepthOfScreen(DefaultScreenOfDisplay(disp)),
		      ZPixmap, 0, NULL, I->ocols, I->orows, 32, 0);

  imgdata = (byte *)malloc((size_t)(I->orows*timg->bytes_per_line));

  timg->data = (char *) imgdata;

  bptr = I->image; xptr = imgdata;
  for(i=0;i<I->orows;i++,xptr += timg->bytes_per_line)
    for(j=0, ip = xptr; j<I->ocols; j++)
      {
	r = *bptr++; g = *bptr++; b = *bptr++;
	
	if (rshift<0) 
	  r = r << (-rshift);
	else 
	  r = r >> rshift;
	if (gshift<0) 
	  g = g << (-gshift);
	else 
	  g = g >> gshift;
	if (bshift<0) 
	  b = b << (-bshift);
	else 
	  b = b >> bshift;

	r = r & rmask;
	g = g & gmask;
	b = b & bmask; 

	xcol = r | g | b;
	*ip++ = xcol & 0xff;
      }

  XtMapWidget(toplevel);

  XtAppMainLoop(context);

  return 0; /* Make -Wall happy! */
}


