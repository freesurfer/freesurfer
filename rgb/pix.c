/*
 *  getpix and putpix -
 *
 *        Paul Haeberli - 1984
 *
 */
#include  <stdio.h>
#include  <gl/image.h>

#undef getpix
#undef putpix

int getpix(IMAGE *image)
{
    if(--(image)->cnt>=0)
      return (int)(*(image)->ptr++);
    else
  return ifilbuf(image);
}

unsigned int putpix(IMAGE *image, unsigned int pix)
{
    if(--(image)->cnt>=0)
        return (unsigned int)(*(image)->ptr++ = pix);
    else
  return iflsbuf(image,pix);
}
