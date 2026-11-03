/**
 * @brief Misc utility class.
 *
 */
/*
 * Original Author: Ruopeng Wang
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */

#include "MyUtils.h"
#include <QFile>
#include "exif.h"

typedef struct
{
  int xl, xr, y, dy;
}
LINESEGMENT;

#define MAXDEPTH 10000

#define PUSH(XL, XR, Y, DY) \
  if( sp < stack+MAXDEPTH && Y+(DY) >= min_y && Y+(DY) <= max_y ) \
{ sp->xl = XL; sp->xr = XR; sp->y = Y; sp->dy = DY; sp++; }

#define POP(XL, XR, Y, DY) \
{ sp--; XL = sp->xl; XR = sp->xr; Y = sp->y + (DY = sp->dy); }

void MyUtils::FloodFill(char* data, int x, int y,
                        int dim_x, int dim_y,
                        int fill_value, int border_value)
{
  int left, x1, x2, dy;
  LINESEGMENT stack[MAXDEPTH], *sp = stack;

  if (data[y*dim_x+x] == border_value || data[y*dim_x+x] == fill_value)
  {
    return;
  }

  int min_x = 0, max_x = dim_x-1, min_y = 0, max_y = dim_y-1;
  if (x < min_x || x > max_x || y < min_y || y > max_y)
  {
    return;
  }

  PUSH(x, x, y, 1);        /* needed in some cases */
  PUSH(x, x, y+1, -1);    /* seed segment (popped 1st) */

  while (sp > stack )
  {
    POP(x1, x2, y, dy);

    for (x = x1; x >= min_x &&
         data[y*dim_x+x] != border_value && data[y*dim_x+x] != fill_value; x--)
    {
      data[y*dim_x+x] = fill_value;
    }

    if ( x >= x1 )
    {
      goto SKIP;
    }

    left = x+1;
    if ( left < x1 )
    {
      PUSH(left, x1-1, y, -dy);  /* leak on left? */
    }

    x = x1+1;

    do
    {
      for (; x<=max_x &&
           data[y*dim_x+x] != border_value && data[y*dim_x+x] != fill_value; x++)
      {
        data[y*dim_x+x] = fill_value;
      }

      PUSH(left, x-1, y, dy);

      if (x > x2+1)
      {
        PUSH(x2+1, x-1, y, -dy);  /* leak on right? */
      }

SKIP:
      for (x++; x <= x2 &&
           (data[y*dim_x+x] == border_value || data[y*dim_x+x] == fill_value); x++)
      {
        ;
      }

      left = x;
    }
    while (x <= x2);
  }
}

int MyUtils::GetExifOrientation(const QString& filename)
{
  int orientation = 0;
  QFile file(filename);
  if (file.open(QIODevice::ReadOnly))
  {
    Exif exif;
    exif.getExifOrientation(file, &orientation);
  }
  return orientation;
}

QImage MyUtils::ApplyExifOrientation(const QImage& image_in, int orientation)
{
  QImage image = image_in;
  if (!image.isNull() && orientation > 1)
  {
    QTransform t;
    if (orientation == 7 || orientation == 8)
      t.rotate(-90);
    else if (orientation == 3 || orientation == 4)
      t.rotate(180);
    else if (orientation == 5 || orientation == 6)
      t.rotate(90);
    image = image.transformed(t, Qt::SmoothTransformation);
  }
  return image;
}

QImage MyUtils::ReadImageWithExifAwareness(const QString& filename)
{
  QImage image(filename);
  if (!image.isNull())
    image = ApplyExifOrientation(image, GetExifOrientation(filename));
  return image;
}
