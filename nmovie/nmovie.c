#include <stdio.h>
#include <image.h>

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
/* -------- End Prototypes -------- */

/* -------- Global Variables ------- */
XtAppContext context;
Display *disp;
Widget frame,toplevel,quit_bt,canvas,buttons,stop_bt;
int shmext,dgaext;
GC theGC;
/* -------- End Global Variables ------- */

/* -------- Static Variables -------- */
static XtActionsRec actions[] = {
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
  "*Canvas.width:          256",
  "*Canvas.height:         256",
  NULL
};
/* -------- End Static Variables -------- */

void useage(void)
{
  exit(1);
}

void quit_proc(Widget w, XEvent *event, String *pars, Cardinal *npars)
{
  XtDestroyWidget(toplevel);
  exit(0);
}

int main(int argc, char **argv)
{
  int ebase,errbase,flags,pfmt=0,rows=0,cols=0,i,cnt,screenno;
  IMAGE *I,*I2;
  XVisualInfo vi;
  Colormap colormap;
  XGCValues xgcv;
  XImage *timg;

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

  canvas = XtVaCreateManagedWidget("Canvas", simpleWidgetClass, frame, NULL);


  XtInstallAllAccelerators(frame,toplevel);
  XtRealizeWidget(toplevel);
  theGC = XCreateGC(disp, XtWindow(canvas), 0L, &xgcv);
  timg = XCreateImage(disp,vi.visual,
		      24,
		      ZPixmap, 0, I->image, I->ocols, I->orows, 8,
		      3*I->ocols);

  XPutImage(disp, XtWindow(canvas), theGC, timg, 0, 0, 0, 0, I->ocols, 
	    I->orows);

  XtMapWidget(toplevel);

  XtAppMainLoop(context);

  return 0; /* Make -Wall happy! */
}


