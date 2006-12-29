/*
 *  isetname and isetcolormap -
 *
 *        Paul Haeberli - 1984
 *
 */
#include  <stdio.h>
#include  <string.h>
#include  "rgb_image.h"

void isetname(RGB_IMAGE *image, char *name) {
  strncpy(image->name,name,80);
}

void isetcolormap(RGB_IMAGE *image, int colormap) {
  image->colormap = colormap;
}
