#ifndef xTypes_H
#define xTypes_H


typedef unsigned char tBoolean;
typedef long          tSignature;

typedef struct {
  int mnX, mnY;
} xPoint2n, *xPoint2nRef;

typedef struct {
  float mfX, mfY;
} xPoint2f, *xPoint2fRef;

#endif
