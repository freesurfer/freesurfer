#ifndef xTypes_H
#define xTypes_H

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif


typedef unsigned char tBoolean;
typedef long          tSignature;

typedef struct {
  int mnX, mnY;
} xPoint2n, *xPoint2nRef;

typedef struct {
  float mfX, mfY;
} xPoint2f, *xPoint2fRef;

typedef struct {
  int mnRed, mnGreen, mnBlue;
} xColor3n, *xColor3nRef;

typedef struct {
  float mfRed, mfGreen, mfBlue;
} xColor3f, *xColor3fRef;

typedef enum {
  tAxis_X = 0,
  tAxis_Y,
  tAxis_Z
} tAxis;

#endif
