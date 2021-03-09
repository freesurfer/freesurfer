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


#include "mri.h"
#include "gca.h"

#ifndef TOPO_INCLUDED
#define TOPO_INCLUDED

//TOPOLOGICAL / CONNECTIVITY STUFF
int connectivityNumber(int connectivity);
int associatedConnectivity(int connectivity);

typedef unsigned char Nbh[3][3][3];
typedef unsigned char LGNBH[5][5][5];

Nbh* loadNbh2(unsigned char*** im,int x,int y, int z,int label);
Nbh* loadNbh(MRI* mri,Nbh* nbh_dst,int i,int j, int k,int label);
Nbh* loadSNbh(MRI* mri,Nbh* nbh_dst,int i,int j, int k,int label);
Nbh* reverseNbh(Nbh* nbh_src,Nbh *nbh_dst);
Nbh *N_6_1(Nbh* nbh_src,Nbh* nbh_dst);
Nbh* N_6_2(Nbh* nbh_src,Nbh* nbh_dst);
Nbh* N_6_3(Nbh* nbh_src,Nbh* nbh_dst);
Nbh *N_18_1(Nbh* nbh_src,Nbh* nbh_dst);
Nbh* N_18_2(Nbh* nbh_src,Nbh* nbh_dst);
Nbh *N_26_1(Nbh* nbh_src,Nbh* nbh_dst);
int checkTn(Nbh *nbh_src,Nbh *nbh_dst,int connectivity);
int checkSimple(Nbh *nbh,int connectivity);
int checkSC(MRI* mri,int i0, int j0, int k0,int i1,int j1,int k1,int inside_label,int outside_label,int connectivity);
int checkWC(MRI* mri,int i0, int j0, int k0,int i1,int j1,int k1,int inside_label,int outside_label,int connectivity);
int checkNbh(MRI *mri,int i,int j,int k,int label,int connectivity);
int checkSNbh(MRI *mri,int i,int j,int k,int label,int connectivity);
Nbh* Nnk(Nbh* nbh_src,Nbh *nbh_dst,int connectivity);
int checkSP(Nbh *fgnbh_src,Nbh *fgnbh_dst,int *fgtp,Nbh *bgnbh_src,Nbh *bgnbh_dst,int *bgtp,int connectivity);

#define MAX_NUMBER_OF_LABELS 50

//different mode of corrections
#define NORMAL_MODE 0
#define VOXEL_MODE 1
#define PROB_MODE 2
#define PROB_MAP_MODE 3
#define MAP_MODE 4

typedef struct MRI_TOPOLOGY_PARMS
{
  //connectivity type (1, 2, 3 or 4)
  int connectivity;
  //label information
  int nlabels;
  int labels[MAX_NUMBER_OF_LABELS];

  //correction mode (VOXEL, PRIORS, ...)
  int mode;
  int background_priority;
  int only;
  int using_gca_maps;
  char *transform_fname,*gca_fname;
  char* prior_map_file;
  float alpha,beta;
  int guess_initial_segmentation;
  //just in case, they are already allocated...
  GCA *gca;
  TRANSFORM *transform;

  //surfaces information
  int generate_surface;
  char *initial_surface_file,*final_surface_file;
  int tesselation_mode;
  int MarchingCubes;

  //debugging information
  char *debugging_map_folder;

  int verbose_mode;

}
MRI_TOPOLOGY_PARMS;

MRI_TOPOLOGY_PARMS* MRI_TOPOLOGY_PARMSalloc(void);
MRI *MRIcorrectTopology(MRI *mri_orig,MRI *mri_seg, MRI*mri_output,MRI_TOPOLOGY_PARMS *mritopparms);

#endif



