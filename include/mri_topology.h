#include "mri.h"

#ifndef TOPO_INCLUDED
#define TOPO_INCLUDED

//TOPOLOGICAL / CONNECTIVITY STUFF
int connectivityNumber(int connectivity);
int associatedConnectivity(int connectivity);

typedef unsigned char Nbh[3][3][3];

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
int checkSC(MRI* mri,int i0, int j0, int k0,int i1,int j1,int k1,int inside_label,int outside_label,int connectivity);
int checkWC(MRI* mri,int i0, int j0, int k0,int i1,int j1,int k1,int inside_label,int outside_label,int connectivity);
int checkNbh(MRI *mri,int i,int j,int k,int label,int connectivity);
int checkSNbh(MRI *mri,int i,int j,int k,int label,int connectivity);
Nbh* Nnk(Nbh* nbh_src,Nbh *nbh_dst,int connectivity);
int checkSP(Nbh *fgnbh_src,Nbh *fgnbh_dst,int *fgtp,Nbh *bgnbh_src,Nbh *bgnbh_dst,int *bgtp,int connectivity);
#endif
