/*
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
 */


#ifndef MRI_CONFORM_H
#define MRI_CONFORM_H


MRI *MRIconform(MRI *mri);
MATRIX *MRIgetConformMatrix(MRI *mri);
MRI *MRIconformedTemplate(MRI *mri, int conform_width, double conform_size, int KeepDC);

/*  EOF  */

#endif
