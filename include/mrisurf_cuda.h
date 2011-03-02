/**
 * @file  mrisurf_cuda.h
 * @brief headers for mrisurf related CUDA functions
 *
 */
/*
 * Original Author: Thomas Witzel
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.4 $
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

#ifndef MRI_CUDA_SURF_H
#define MRI_CUDA_SURF_H
#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

typedef struct MRI_CUDA_SURFACE_
{
	
	unsigned int nvertices;     /* number of vertices */
	float4 *h_V;                /* host memory for vertices */
	float4 *d_V;                /* device memory for vertices */

	unsigned int ntotalneighbors;     /* number of "total" neighborhood */
	unsigned int *h_TotalNeighborArray;   /* neighborhood array */
	unsigned int *d_TotalNeighborArray;
	
	unsigned int *h_TotalNeighborOffsets; /* offsets into the neighbor array */
	unsigned int *d_TotalNeighborOffsets;
	
	unsigned int *h_NTotalNeighbors;      /* number of neighbors for each vertex */
	unsigned int *d_NTotalNeighbors;
	
	unsigned int *h_NNeighbors;           /* number of immediate neighborhood members */
	unsigned int *d_NNeighbors;           /* device version of " " */
	
	float *h_Distances;              /* distance vector */
	float *d_Distances;
	
	float *d_DistancesOrig;          /* original distance vector */
	
	float4 *h_D;      /* host memory for gradients */
	float4 *d_D;      /* device memory for gradients */
	float4 *d_TD;     /* device memory for temporary gradients */
	
} MRI_CUDA_SURFACE;

#ifdef __cplusplus
extern "C" {
#endif
	/* utility routines */
	void MRISCinitSurface(MRI_CUDA_SURFACE *mrics);
	void MRISCcleanupSurface(MRI_CUDA_SURFACE *mrics);
	void MRISCallocVertices(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCuploadVertices(MRI_CUDA_SURFACE *mrics,MRI_SURFACE *mris);
	void MRISCdownloadDistances(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCallocDistances(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCuploadTotalNeighborArray(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCallocTotalNeighborArray(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCcheckCUDAError(const char *msg);
	void MRISCallocGradients(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCfreeGradients(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCuploadGradients(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCdownloadGradients(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCdeviceInfo();

	void MRISCallocDistOrigArray(MRI_CUDA_SURFACE *mrics);
	void MRISCuploadDistOrigArray(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	
	/* computational routines */
	void MRISCcomputeVertexDistances(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris);
	void MRISCaverageGradients(MRI_CUDA_SURFACE *mrisc, MRI_SURFACE *mris, unsigned int niter);
	float MRISCcomputeDistanceError(MRI_CUDA_SURFACE *mrisc, float dist_scale);
	
#ifdef __cplusplus
}
#endif

#endif /* MRI_CUDA_SURF_H */
