#ifndef XWINDOW_H
#define XWINDOW_H

/*###########################################################################*/
/* tk.h: marty: tk=>tko, TK_=>TKO_; add bstring.h (ffs) */
#include <GL/gl.h>
#include <GL/glu.h>

/*
** Nano Window Toolkit.
*/

/*
** Window Types
*/

#define TKO_RGB   0
#define TKO_INDEX 1
#define TKO_SINGLE  0
#define TKO_DOUBLE  2
#define TKO_DIRECT  0
#define TKO_INDIRECT  4
#define TKO_ACCUM 8
#define TKO_ALPHA 16
#define TKO_DEPTH 32
#define TKO_OVERLAY 64
#define TKO_UNDERLAY  128
#define TKO_STENCIL 512

/*
** Display Mode Selection Criteria
*/

enum {
    TKO_USE_ID = 1,
    TKO_EXACT_MATCH,
    TKO_MINIMUM_CRITERIA
};

/* 
** Window Masks
*/

#define TKO_IS_RGB(x)   (((x) & TKO_INDEX) == 0)
#define TKO_IS_INDEX(x)   (((x) & TKO_INDEX) != 0)
#define TKO_IS_SINGLE(x)  (((x) & TKO_DOUBLE) == 0)
#define TKO_IS_DOUBLE(x)  (((x) & TKO_DOUBLE) != 0)
#define TKO_IS_DIRECT(x)  (((x) & TKO_INDIRECT) == 0)
#define TKO_IS_INDIRECT(x)  (((x) & TKO_INDIRECT) != 0)
#define TKO_HAS_ACCUM(x)  (((x) & TKO_ACCUM) != 0)
#define TKO_HAS_ALPHA(x)  (((x) & TKO_ALPHA) != 0)
#define TKO_HAS_DEPTH(x)  (((x) & TKO_DEPTH) != 0)
#define TKO_HAS_OVERLAY(x)  (((x) & TKO_OVERLAY) != 0)
#define TKO_HAS_UNDERLAY(x) (((x) & TKO_UNDERLAY) != 0)
#define TKO_HAS_STENCIL(x)  (((x) & TKO_STENCIL) != 0)

/*
** Windowing System Specific Gets
*/

enum {
    TKO_X_DISPLAY = 1,
    TKO_X_WINDOW,
    TKO_X_SCREEN
};

/*
** Event Status
*/

#define TKO_LEFTBUTTON    1
#define TKO_RIGHTBUTTON   2
#define TKO_MIDDLEBUTTON    4
#define TKO_SHIFT   1
#define TKO_CONTROL   2

/* 
** Key Codes
*/

#define TKO_RETURN    0x0D
#define TKO_ESCAPE    0x1B
#define TKO_SPACE   0x20
#define TKO_LEFT      0x25
#define TKO_UP      0x26
#define TKO_RIGHT   0x27
#define TKO_DOWN      0x28
#define TKO_A     'A'
#define TKO_B     'B'
#define TKO_C     'C'
#define TKO_D     'D'
#define TKO_E     'E'
#define TKO_F     'F'
#define TKO_G     'G'
#define TKO_H     'H'
#define TKO_I     'I'
#define TKO_J     'J'
#define TKO_K     'K'
#define TKO_L     'L'
#define TKO_M     'M'
#define TKO_N     'N'
#define TKO_O     'O'
#define TKO_P     'P'
#define TKO_Q     'Q'
#define TKO_R     'R'
#define TKO_S     'S'
#define TKO_T     'T'
#define TKO_U     'U'
#define TKO_V     'V'
#define TKO_W     'W'
#define TKO_X     'X'
#define TKO_Y     'Y'
#define TKO_Z     'Z'
#define TKO_a     'a'
#define TKO_b     'b'
#define TKO_c     'c'
#define TKO_d     'd'
#define TKO_e     'e'
#define TKO_f     'f'
#define TKO_g     'g'
#define TKO_h     'h'
#define TKO_i     'i'
#define TKO_j     'j'
#define TKO_k     'k'
#define TKO_l     'l'
#define TKO_m     'm'
#define TKO_n     'n'
#define TKO_o     'o'
#define TKO_p     'p'
#define TKO_q     'q'
#define TKO_r     'r'
#define TKO_s     's'
#define TKO_t     't'
#define TKO_u     'u'
#define TKO_v     'v'
#define TKO_w     'w'
#define TKO_x     'x'
#define TKO_y     'y'
#define TKO_z     'z'
#define TKO_0     '0'
#define TKO_1     '1'
#define TKO_2     '2'
#define TKO_3     '3'
#define TKO_4     '4'
#define TKO_5     '5'
#define TKO_6     '6'
#define TKO_7     '7'
#define TKO_8     '8'
#define TKO_9     '9'

/*
** Color Macros
*/

enum {
    TKO_BLACK = 0,
    TKO_RED,
    TKO_GREEN,
    TKO_YELLOW,
    TKO_BLUE,
    TKO_MAGENTA,
    TKO_CYAN,
    TKO_WHITE
};

extern float tkoRGBMap[8][3];

#define TKO_SETCOLOR(x, y) (TKO_IS_RGB((x)) ? \
               glColor3fv(tkoRGBMap[(y)]) : glIndexf((y)))

/*
** RGB Image Structure
*/

typedef struct _TKO_RGBImageRec {
    GLint sizeX, sizeY;
    unsigned char *data;
} TKO_RGBImageRec;

/*
** Prototypes
*/

extern GLenum tkoInitPixmap(long x, long y);
extern GLenum tkoInitDisplay(void);
extern void tkoInitDisplayMode(GLenum);
extern void tkoInitDisplayModePolicy(GLenum type);
extern GLenum tkoInitDisplayModeID(GLint id);
extern void tkoInitPosition(int, int, int, int);
extern GLenum tkoInitWindow(char *);
extern void tkoCloseWindow(void);
extern void tkoQuit(void);

extern GLenum tkoSetWindowLevel(GLenum);
extern void tkoSwapBuffers(void);

extern void tkoExec(void);
extern void tkoExposeFunc(void (*)(int, int));
extern void tkoReshapeFunc(void (*)(int, int));
extern void tkoDisplayFunc(void (*)(void));
extern void tkoKeyDownFunc(GLenum (*)(int, GLenum));
extern void tkoMouseDownFunc(GLenum (*)(int, int, GLenum));
extern void tkoMouseUpFunc(GLenum (*)(int, int, GLenum));
extern void tkoMouseMoveFunc(GLenum (*)(int, int, GLenum));
extern void tkoIdleFunc(void (*)(void));

extern int tkoGetColorMapSize(void);
extern void tkoGetMouseLoc(int *, int *);
extern void tkoGetSystem(GLenum, void *);
extern GLenum tkoGetDisplayModePolicy(void);
extern GLint tkoGetDisplayModeID(void);
extern GLenum tkoGetDisplayMode(void);
extern long tkoGetContext(void);

extern void tkoSetOneColor(int, float, float, float);
extern void tkoSetFogRamp(int, int);
extern void tkoSetGreyRamp(void);
extern void tkoSetRGBMap(int, float *);
extern void tkoSetOverlayMap(int, float *);

extern void tkoNewCursor(GLint, GLubyte *, GLubyte *, GLenum, GLenum,
      GLint, GLint);
extern void tkoSetCursor(GLint);

extern TKO_RGBImageRec *tkoRGBImageLoad(char *);

extern GLenum tkoCreateStrokeFont(GLuint);
extern GLenum tkoCreateOutlineFont(GLuint);
extern GLenum tkoCreateFilledFont(GLuint);
extern GLenum tkoCreateBitmapFont(GLuint);
extern void tkoDrawStr(GLuint, char *);

extern void tkoWireSphere(GLuint, float);
extern void tkoSolidSphere(GLuint, float);
extern void tkoWireCube(GLuint, float);
extern void tkoSolidCube(GLuint, float);
extern void tkoWireBox(GLuint, float, float, float);
extern void tkoSolidBox(GLuint, float, float, float);
extern void tkoWireTorus(GLuint, float, float);
extern void tkoSolidTorus(GLuint, float, float);
extern void tkoWireCylinder(GLuint, float, float);
extern void tkoSolidCylinder(GLuint, float, float);
extern void tkoWireCone(GLuint, float, float);
extern void tkoSolidCone(GLuint, float, float);

/*###########################################################################*/
/* private.h */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <GL/glx.h>
/*#include <bstring.h>*/

#if defined(__cplusplus) || defined(c_plusplus)
#define class c_class
#endif


typedef struct _WINDOW_REC {
    int x, y, w, h;
    GLenum type;
    GLenum dmPolicy;
    Window wMain, wOverlay;
    XVisualInfo *vInfoMain, *vInfoOverlay;
    Colormap cMapMain, cMapOverlay;
    GLXContext cMain, cOverlay;
} WINDOW_REC;


extern Display *xDisplay;
extern int xScreen; 
extern Window wRoot;
extern WINDOW_REC w;
extern Atom deleteWindowAtom;

extern void (*ExposeFunc)(int, int);
extern void (*ReshapeFunc)(int, int);
extern void (*DisplayFunc)(void);
extern GLenum (*KeyDownFunc)(int, GLenum);
extern GLenum (*MouseDownFunc)(int, int, GLenum);
extern GLenum (*MouseUpFunc)(int, int, GLenum);
extern GLenum (*MouseMoveFunc)(int, int, GLenum);
extern void (*IdleFunc)(void);

extern GLenum drawAllowFlag;

extern int cursorNum;

extern WINDOW_REC w ;
extern Display *xDisplay ;

#endif
