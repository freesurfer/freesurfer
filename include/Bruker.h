// 
// Bruker.h
//
// created: y.tosa
// date   : Aug 25th, 2003
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: tosa $
// Revision Date  : $Date: 2003/08/29 16:43:55 $
// Revision       : $Revision: 1.1 $

#ifndef c_bruker_h
#define c_bruker_h

// specify directory name
MRI *brukerRead(char *fname, int read_volume);

// utility routines
// check and create filenames for Bruker volume from directory name
int checkBrukerFiles(char *fname, char *methodFile, char *acqpFile, char *dataFile, char *d3procFile, 
		     int flag); // whether to say something(1) or not
// gives width, height, depth, type, nframes
int readBrukerD3proc(char *d3procFile, int *pwidth, int *pheight, int *pdepth, int *ptype, int *pnframes);
// allocate MRI
MRI *readBrukerMethod(char *methodFile, char *d3procFile);
// modify MRI
int readBrukerAcqp(MRI *mri, char *acqpFile);
int readBrukerVolume(MRI *mri, char *dataFile);
// check whether fname is a bruker directory
int is_bruker(char *fname);

#endif
