//
// ctrpoints.h
//
// purpose: read and write control points
//
#ifndef ctrpoints_h
#define ctrpoints_h

typedef struct { Real x; Real y; Real z ;} MPoint;

// reading control points 
// returning array of MGHPoint
// if useRealRAS = 1, then it is in scanner RAS
// if useRealRAS = 0, then it is in surface RAS
MPoint *MRIreadControlPoints(const char *fname, int *count, int *useRealRAS);

// writing control points 
// returning array of MGHPoint
// we will write whether they are in scannerRAS or surfaceRAS
// Note that use must tell which coordinate system it is.
int MRIwriteControlPoints(MPoint *pointArray, int count, int useRealRAS, char *fname);

#endif // inclusion guard
