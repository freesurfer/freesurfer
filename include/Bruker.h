// 
// Bruker.h
//
// created: y.tosa
// date   : Aug 25th, 2003
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: tosa $
// Revision Date  : $Date: 2003/09/04 15:48:00 $
// Revision       : $Revision: 1.3 $

#ifndef c_bruker_h
#define c_bruker_h

typedef struct {
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
} BrukerTransform;
  
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
