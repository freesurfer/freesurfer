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

#ifndef MyUtils_h
#define MyUtils_h

#include <QImage>

class MyUtils
{
public:
  static void FloodFill( char* data, int x, int y, int dim_x, int dim_y, int fill_value, int border_value );
  static QImage ReadImageWithExifAwareness(const QString& filename);

  // EXIF orientation (1-8, 0 if none) of filename. Separate from
  // ReadImageWithExifAwareness() so an image derived from a photo, such as a
  // mask with no EXIF of its own, can get the photo's rotation.
  static int GetExifOrientation(const QString& filename);
  // Applies the same rotation ReadImageWithExifAwareness() would for the
  // given orientation code -- see GetExifOrientation().
  static QImage ApplyExifOrientation(const QImage& image, int orientation);
};

#endif
