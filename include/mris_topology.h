// include guard
#ifndef MRIS_TOPOLOGY_INCLUDED_H
#define MRIS_TOPOLOGY_INCLUDED_H

// The following is usable from C
#ifdef __cplusplus
extern "C" {
#endif

#include "mrisurf.h"

	//functions
	bool MRIScorrectDefect(MRIS *mris, int defect_number,TOPOFIX_PARMS &parms);
  bool MRIScorrectPatchTopology(MRIS* &mris,TOPOFIX_PARMS &parms);
  bool MRISincreaseEuler(MRIS* &mris,TOPOFIX_PARMS &parms);
	int MRISgetEuler(MRIS *mris, int defect_number = -1);
	int MRISgetEulerNumber(const MRIS *mris, const int *list_of_faces, int nfs);
	MRIP* MRIPextractFromMRIS(MRIS *mris, int defect_number);
	void MRISinitSurface(MRIS *mris);
	bool MRISaddMRIP(MRIS *mris_dst, MRIP *mrip);
	MRIS *MRISduplicateOver(MRIS *mris,int mode = 0);
	void MRIScopyHeader(MRIS *mris_src,MRIS *mris_dst);
 
#ifdef __cplusplus
}
#endif



// C++ portion starts here
#ifdef __cplusplus

extern "C" {
#include "mrisurf.h"
#include "error.h"
}
#include "topology/surface.h"


MRIP *MRIPalloc(int nvertices, int nfaces);
void MRIPfree(MRIP **mrip);
MRIP *MRIPclone(MRIP *src);

Surface *MRIStoSurface(MRIS *mris);
MRIS * SurfaceToMRIS(Surface *surface, MRIS *mris);


#endif

#endif
