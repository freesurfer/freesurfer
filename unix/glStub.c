#ifndef HAVE_OPENGL
// empty stubs to allow building against VTK libs
// on cluster platforms which do not have graphics libs
void glPushClientAttrib(){}
void glPushAttrib(){}
void glPixelStorei(){}
void glEnable(){}
void glBlendFunc(){}
void glDisable(){}
void glPopAttrib(){}
void glPopClientAttrib(){}
void glGetFloatv(){}
void glBitmap(){}
void glDrawPixels(){}
#endif
