//
// rgbutils.h
//
#ifndef c_rgbutils_h
#define c_rgbutils_h

int    RGBwrite(IMAGE *I, char *fname, int frame) ;
IMAGE *RGBReadImage(char *fname);
IMAGE *RGBReadHeader(char *fname, IMAGE *);

#endif
