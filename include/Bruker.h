/**
 * @file  Bruker.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
 *    $Revision: 1.5 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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


//
// Bruker.h
//
// created: y.tosa
// date   : Aug 25th, 2003
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: nicks $
// Revision Date  : $Date: 2011/03/02 00:04:09 $
// Revision       : $Revision: 1.5 $

#ifndef c_bruker_h
#define c_bruker_h

typedef struct
{
  int transposition; // 0 no , 1 x-y flipped
  double fov[3];
  int size[3];
  int ft_size[3];
  double vox_size[3]; // calculated from fov, size, ft_size
  double offset[3];   // calculated from vox, ft_size
  double read_offset;
  double phase1_offset;
  double slice_offset;
  double grad_matrix[9];
  int type;           // MRI data type
  int dim;            // record the dimension (for 2d case)
}
BrukerTransform;

// specify directory name
MRI *brukerRead(char *fname, int read_volume);

// utility routines
// check and create filenames for Bruker volume from directory name
int checkBrukerFiles(char *fname, char *methodFile, char *acqpFile, char *dataFile, char *d3procFile,
                     char *recoFile, int flag); // whether to say something(1) or not
// reconstructed image info
int readBrukerD3proc(char *d3procFile, int *pwidth, int *pheight, int *pdepth, int *ptype, int *pnframes);
int readBrukerReco(char *recoFile, BrukerTransform *bTran);
int readBrukerVolume(MRI *mri, char *dataFile);
// acquisition information
int readBrukerAcqp(char *acqpFile, double *TR, double *TE, double *TI, double *flip_angle, BrukerTransform *bTran);
// check whether fname is a bruker directory
int is_bruker(char *fname);

int buildVoxToRASTransform(MRI *mri, BrukerTransform *bTran);

#endif
