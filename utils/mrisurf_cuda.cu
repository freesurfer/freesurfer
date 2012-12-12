/**
 * @file mrisurf_cuda.cu
 * @brief CUDA based support functions for mrisurf.c
 *
 */
/*
 * Original Author: Thomas Witzel
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/12/12 21:18:24 $
 *    $Revision: 1.10 $
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

/* manually this compiles by

   nvcc -I ../include/ -DFS_CUDA -arch=sm_13 -c mrisurf_cuda.cu

*/
#ifdef FS_CUDA

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <unistd.h>

#include <sys/time.h>
#include <thrust/version.h>
#include <thrust/device_ptr.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include <cmath>

/* Note mrisurf.h cannot be before thrust libraries .... */
#include "mrisurf.h"
#include "mrisurf_cuda.h"
#include "devicemanagement.h"

#include <cuda.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>

#define DMALLOC 0

#if DMALLOC
#include "dmalloc.h"
#endif

extern const char* Progname;

// uncomment this to expose code which shows timings of gpu activities:
//#define FS_CUDA_TIMINGS

#ifdef FS_CUDA_TIMINGS
static int timeval_subtract (struct timeval *result,
                             struct timeval * x,
                             struct timeval * y)
{
  /* Perform the carry for the later subtraction by updating y. */
  if (x->tv_usec < y->tv_usec)
  {
    int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
    y->tv_usec -= 1000000 * nsec;
    y->tv_sec += nsec;
  }
  if (x->tv_usec - y->tv_usec > 1000000)
  {
    int nsec = (x->tv_usec - y->tv_usec) / 1000000;
    y->tv_usec += 1000000 * nsec;
    y->tv_sec -= nsec;
  }

  /* Compute the time remaining to wait.
    tv_usec is certainly positive. */
  result->tv_sec = x->tv_sec - y->tv_sec;
  result->tv_usec = x->tv_usec - y->tv_usec;

  /* Return 1 if result is negative. */
  return x->tv_sec < y->tv_sec;
}
#endif // FS_CUDA_TIMINGS

/*
   this function is purely for debug and should be replaced
   with a proper function
   from the freesurfer environment
*/
void MRISCdeviceInfo()
{
  /*
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, 0);
    printf ("      CUDA device %d: %s\n", 0,deviceProp.name);
  */
  AcquireCUDADevice();
}

void MRISCcheckCUDAError(const char *msg)
{
  cudaError_t err = cudaGetLastError();
  if ( err != cudaSuccess)
  {
    fprintf(stderr,"CUDA error: %s: %s.\n",msg,cudaGetErrorString(err));
    exit(-1);
  }
}

void MRISCinitSurface(MRI_CUDA_SURFACE *mrics)
{
  /* initialize the CUDA surface structure */
  mrics->nvertices = 0;
  mrics->h_V = NULL;
  mrics->d_V = NULL;

  mrics->ntotalneighbors = 0;
  mrics->h_TotalNeighborArray = NULL;
  mrics->d_TotalNeighborArray = NULL;
  mrics->h_TotalNeighborOffsets = NULL;
  mrics->d_TotalNeighborOffsets = NULL;
  mrics->h_NTotalNeighbors = NULL;
  mrics->d_NTotalNeighbors = NULL;
  mrics->h_NNeighbors = NULL;
  mrics->d_NNeighbors = NULL;

  mrics->h_Distances = NULL;
  mrics->d_Distances = NULL;

  mrics->d_DistancesOrig = NULL;

  mrics->h_D = NULL;
  mrics->d_D = NULL;
  mrics->d_TD = NULL;
}

void MRISCcleanupSurface(MRI_CUDA_SURFACE *mrics)
{

  if (mrics->h_D != NULL)
  {
    free(mrics->h_D);
    mrics->h_D = NULL;
  }

  if (mrics->d_D != NULL)
  {
    cudaFree(mrics->d_D);
    mrics->d_D = NULL;
  }

  if (mrics->d_TD != NULL)
  {
    cudaFree(mrics->d_TD);
    mrics->d_TD = NULL;
  }

  if (mrics->h_V != NULL)
  {
    free(mrics->h_V);
    mrics->h_V = NULL;
  }

  if (mrics->d_V != NULL)
  {
    cudaFree(mrics->h_V);
    mrics->h_V = NULL;
  }

  if (mrics->h_Distances != NULL)
  {
    free(mrics->h_Distances);
    mrics->h_Distances = NULL;
  }

  if (mrics->d_Distances != NULL)
  {
    cudaFree(mrics->d_Distances);
    mrics->d_Distances = NULL;
  }

  if (mrics->h_TotalNeighborArray)
  {
    free(mrics->h_TotalNeighborArray);
    mrics->h_TotalNeighborArray = NULL;
  }
  if (mrics->d_TotalNeighborArray != NULL)
  {
    cudaFree(mrics->d_TotalNeighborArray);
    mrics->d_TotalNeighborArray = NULL;
  }

  if (mrics->h_TotalNeighborOffsets != NULL)
  {
    free(mrics->h_TotalNeighborOffsets);
    mrics->h_TotalNeighborOffsets = NULL;
  }

  if (mrics->d_TotalNeighborOffsets != NULL)
  {
    cudaFree(mrics->d_TotalNeighborOffsets);
    mrics->d_TotalNeighborArray = NULL;
  }

  if (mrics->h_NTotalNeighbors != NULL)
  {
    free(mrics->h_NTotalNeighbors);
    mrics->h_NTotalNeighbors = NULL;
  }

  if (mrics->d_NTotalNeighbors != NULL)
  {
    cudaFree(mrics->d_NTotalNeighbors);
    mrics->d_NTotalNeighbors = NULL;
  }

  if (mrics->h_NNeighbors != NULL)
  {
    free(mrics->h_NNeighbors);
    mrics->h_NNeighbors = NULL;
  }

  if (mrics->d_NNeighbors != NULL)
  {
    cudaFree(mrics->d_NNeighbors);
    mrics->d_NNeighbors = NULL;
  }

  if (mrics->d_DistancesOrig != NULL)
  {
    cudaFree(mrics->d_DistancesOrig);
    mrics->d_DistancesOrig = NULL;
  }
}

void MRISCallocVertices(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  /* if there is some preexisting memory buffers, free them */
  if (mrics->h_V != NULL)
  {
    free(mrics->h_V);
    mrics->h_V = NULL;
  }

  if (mrics->d_V != NULL)
  {
    cudaFree(mrics->h_V);
    mrics->h_V = NULL;
  }

  /* there should be better error checks here */
  mrics->h_V = (float4 *)malloc(mris->nvertices*sizeof(float4));
  if (mrics->h_V == NULL)
  {
    fprintf(stderr,"MRISallocCudaVertices: cannot allocate host memory !\n");
    exit(-1);
  }
  cudaMalloc((void **)&(mrics->d_V),mris->nvertices*sizeof(float4));
  MRISCcheckCUDAError("MRISallocCudaVertices: "
                      "cannot allocate device memory !");

  mrics->nvertices = mris->nvertices;
}

void MRISCallocDistances(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  /* count the number of total neighbors */
  unsigned int ntotal = 0;
  int vno;

  for (vno=0; vno<mris->nvertices; vno++)
  {
    ntotal += mris->vertices[vno].vtotal;
  }

  /* then allocate */
  if (mrics->h_Distances != NULL)
  {
    free(mrics->h_Distances);
  }

  mrics->h_Distances = (float *)malloc(ntotal*sizeof(float));
  if (mrics->h_Distances == NULL)
  {
    fprintf(stderr,"MRISCallocDistances: cannot allocate host memory !\n");
    exit(-1);
  }
  /* allocate device memory */
  if (mrics->d_Distances != NULL)
  {
    cudaFree(mrics->d_Distances);
  }

  cudaMalloc((void **)&(mrics->d_Distances),ntotal*sizeof(float));
  MRISCcheckCUDAError("MRISCallocDistances: cannot allocate device memory !");
}

void MRISCuploadVertices(MRI_CUDA_SURFACE *mrics,MRI_SURFACE *mris)
{
  int vno;
  float4 *h_V;
#ifdef FS_CUDA_TIMINGS
  struct timeval tv1,tv2,result;
#endif // FS_CUDA_TIMINGS

  /* if it appears that either array is not allocated, reallocate */
  if (mrics->h_V == NULL
      || mrics->d_V == NULL
      || mrics->nvertices != mris->nvertices)
  {
    MRISCallocVertices(mrics,mris);
  }

  h_V = mrics->h_V;

#ifdef FS_CUDA_TIMINGS
  gettimeofday(&tv1,NULL);
#endif // FS_CUDA_TIMINGS

  // copy the data from the vertex structure into the host array
  for (vno=0; vno<mris->nvertices; vno++)
  {
    h_V[vno].x = mris->vertices[vno].x;
    h_V[vno].y = mris->vertices[vno].y;
    h_V[vno].z = mris->vertices[vno].z;
    h_V[vno].w = 0;
  }

  // finally transfer the data
  cudaMemcpy(mrics->d_V,
             mrics->h_V,
             mris->nvertices*sizeof(float4),
             cudaMemcpyHostToDevice);
  MRISCcheckCUDAError("MRISuploadVertices: cannot transfer data !");

#ifdef FS_CUDA_TIMINGS
  gettimeofday(&tv2,NULL);
  timeval_subtract(&result,&tv2,&tv1);
  printf("MRISuploadVertices: %ld ms\n",
         result.tv_sec*1000+result.tv_usec/1000);
  fflush(stdout);
#endif // FS_CUDA_TIMINGS
}

void MRISCdownloadDistances(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  VERTEX *v;
  int vno; /* n; */
#ifdef FS_CUDA_TIMINGS
  timeval tv1,tv2,result;
#endif // FS_CUDA_TIMINGS

  if (mrics->h_Distances == NULL || mrics->d_Distances == NULL)
  {
    fprintf(stderr,"MRISCdownloadDistances: "
            "Cannot operate on unallocated memory !\n");
  }
#ifdef FS_CUDA_TIMINGS
  gettimeofday(&tv1,NULL);
#endif // FS_CUDA_TIMINGS
  // download the data from the device
  cudaMemcpy(mrics->h_Distances,
             mrics->d_Distances,
             mrics->ntotalneighbors*sizeof(float),
             cudaMemcpyDeviceToHost);

  unsigned int ctr = 0;
  // now sort into structure
  for (vno=0; vno<mris->nvertices; vno++)
  {
    v = &(mris->vertices[vno]);
    /* Lets see whether we can speed this up with memcpy
      yep, memcpy seems about 35% faster here than the loop below
    */
    memcpy(v->dist,mrics->h_Distances+ctr,v->vtotal*sizeof(float));
    ctr += v->vtotal;
    /*
    for(n=0; n<v->vtotal; n++,ctr++)
     v->dist[n] = mrics->h_Distances[ctr];
    */
  }
#ifdef FS_CUDA_TIMINGS
  gettimeofday(&tv2,NULL);
  timeval_subtract(&result,&tv2,&tv1);
  printf("MRISdownloadDistances: %ld ms\n",
         result.tv_sec*1000+result.tv_usec/1000);
  fflush(stdout);
#endif // FS_CUDA_TIMINGS
}

void MRISCallocTotalNeighborArray(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  unsigned int ntotal = 0;
  int vno;

  for (vno = 0; vno < mris->nvertices; vno++)
  {
    ntotal += mris->vertices[vno].vtotal;
  }

  mrics->ntotalneighbors = ntotal;

  if (mrics->h_TotalNeighborArray != NULL)
  {
    free(mrics->h_TotalNeighborArray);
  }
  mrics->h_TotalNeighborArray =
    (unsigned int *)malloc(ntotal*sizeof(unsigned int));
  if (mrics->h_TotalNeighborArray == NULL)
  {
    fprintf(stderr,"MRISCallocTotalNeighborArray: "
            "Cannot allocate host memory for neighbor array !\n");
    exit(-1);
  }

  if (mrics->d_TotalNeighborArray != NULL)
  {
    cudaFree(mrics->d_TotalNeighborArray);
  }

  cudaMalloc((void **)&(mrics->d_TotalNeighborArray),
             ntotal*sizeof(unsigned int));
  if (mrics->d_TotalNeighborArray == NULL)
  {
    fprintf(stderr,"MRISCallocTotalNeighborArray: "
            "Cannot allocate device memory for neighbor array !\n");
    exit(-1);
  }

  if (mrics->h_NTotalNeighbors != NULL)
  {
    free(mrics->h_NTotalNeighbors);
  }

  mrics->h_NTotalNeighbors =
    (unsigned int *)malloc(mrics->nvertices*sizeof(unsigned int));
  if (mrics->h_NTotalNeighbors == NULL)
  {
    fprintf(stderr,"MRISCallocTotalNeighborArray: "
            "Cannot allocate host memory for n-neighbor array !\n");
    exit(-1);
  }

  if (mrics->d_NTotalNeighbors != NULL)
  {
    cudaFree(mrics->d_NTotalNeighbors);
  }
  cudaMalloc((void **)&(mrics->d_NTotalNeighbors),
             mrics->nvertices*sizeof(unsigned int));
  if (mrics->d_NTotalNeighbors == NULL)
  {
    fprintf(stderr,"MRISCallocTotalNeighborArray: "
            "Cannot allocate device memory for n-neighbor array !\n");
    exit(-1);
  }

  if (mrics->h_NNeighbors != NULL)
  {
    free(mrics->h_NNeighbors);
  }

  mrics->h_NNeighbors =
    (unsigned int *)malloc(mrics->nvertices*sizeof(unsigned int));
  if (mrics->h_NNeighbors == NULL)
  {
    fprintf(stderr,"MRISCallocTotalNeighborArray: "
            "Cannot allocate host memory for n-neighbor array !\n");
    exit(-1);
  }

  if (mrics->d_NNeighbors != NULL)
  {
    cudaFree(mrics->d_NNeighbors);
  }
  cudaMalloc((void **)&(mrics->d_NNeighbors),
             mrics->nvertices*sizeof(unsigned int));
  if (mrics->d_NNeighbors == NULL)
  {
    fprintf(stderr,"MRISCallocTotalNeighborArray: "
            "Cannot allocate device memory for n-neighbor array !\n");
    exit(-1);
  }


  if (mrics->h_TotalNeighborOffsets != NULL)
  {
    free(mrics->h_TotalNeighborOffsets);
  }
  mrics->h_TotalNeighborOffsets =
    (unsigned int *)malloc(mrics->nvertices*sizeof(unsigned int));
  if (mrics->h_TotalNeighborOffsets == NULL)
  {
    fprintf(stderr,"MRISCallocTotalNeighborArray: "
            "Cannot allocate host memory for neighbor-offsets array !\n");
    exit(-1);
  }

  if (mrics->d_TotalNeighborOffsets != NULL)
  {
    cudaFree(mrics->d_TotalNeighborOffsets);
  }

  cudaMalloc((void **)&(mrics->d_TotalNeighborOffsets),
             mrics->nvertices*sizeof(unsigned int));
  if (mrics->d_TotalNeighborOffsets == NULL)
  {
    fprintf(stderr,"MRISCallocTotalNeighborArray: "
            "Cannot allocate device memory for neighbor-offsets array !\n");
    exit(-1);
  }
}

/* allocate a dist_orig vector if it is not already allocated */
void MRISCallocDistOrigArray(MRI_CUDA_SURFACE *mrics)
{
  unsigned int ntotal = mrics->ntotalneighbors;
  /* allocate device memory */
  if (mrics->d_DistancesOrig != NULL)
  {
    cudaFree(mrics->d_DistancesOrig);
  }

  cudaMalloc((void **)&(mrics->d_DistancesOrig),ntotal*sizeof(float));
  MRISCcheckCUDAError("MRISCallocDistOrigArray: "
                      "cannot allocate device memory !");
}

/* upload the dist_orig vector to the GPU device.
   This assumes that the neighborhood structure is
   already uploaded to some extent
   a download routine is not needed because this is never modified by the GPU
*/
void MRISCuploadDistOrigArray(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  if (mrics->d_DistancesOrig == NULL)
  {
    MRISCallocDistOrigArray(mrics);
  }

  unsigned int ntotal = mrics->ntotalneighbors;
  unsigned int ctr = 0;

  float *tmp = new float[ntotal];
  if (tmp == NULL)
  {
    fprintf(stderr,"MRISCuploadDistOrigArray: "
            "cannot allocate temporary buffer !");
    exit(-1);
  }

  for (int vno = 0; vno < mris->nvertices; vno++)
  {
    VERTEX *v = &(mris->vertices[vno]);
    if (v->ripflag)
    {
      fprintf(stderr,"PANIC: RIPFLAG IS NOT SUPPORTED !\n");
    }

    for (int n=0; n<v->vtotal; n++)
    {
      tmp[ctr++] = v->dist_orig[n];
    }
  }
  // transfer the data
  cudaMemcpy(mrics->d_DistancesOrig,
             tmp,
             mrics->ntotalneighbors*sizeof(float),
             cudaMemcpyHostToDevice);
  MRISCcheckCUDAError("MRISCuploadDistOrigArray: cannot transfer the data !");
  // free the tmp vector
  delete [] tmp;
}

/* upload the neighbor structure to the GPU
  a download routine is not needed because this is never modifed by the GPU
*/
void MRISCuploadTotalNeighborArray(MRI_CUDA_SURFACE *mrics,
                                   MRI_SURFACE *mris)
{
  unsigned int ntotal = 0;
  unsigned int ctr=0;
  int vno,n;
  VERTEX *v;

  for (vno = 0; vno < mris->nvertices; vno++)
  {
    ntotal += mris->vertices[vno].vtotal;
  }

  if (ntotal != mrics->ntotalneighbors
      || mrics->h_TotalNeighborArray == NULL
      || mrics->d_TotalNeighborArray == NULL
      || mrics->h_TotalNeighborOffsets == NULL
      || mrics->d_TotalNeighborOffsets == NULL
      || mrics->h_NTotalNeighbors == NULL
      || mrics->d_NTotalNeighbors == NULL)
  {
    MRISCallocTotalNeighborArray(mrics,mris);
  }

  for (vno = 0; vno < mris->nvertices; vno++)
  {

    v = &(mris->vertices[vno]);
    if (v->ripflag)
    {
      fprintf(stderr,"PANIC: RIPFLAG IS NOT SUPPORTED !\n");

    }

    mrics->h_NTotalNeighbors[vno] = v->vtotal;
    mrics->h_NNeighbors[vno] = v->vnum;

    if (vno == 0)
    {
      mrics->h_TotalNeighborOffsets[vno] = 0;
    }
    else
      mrics->h_TotalNeighborOffsets[vno] =
        mrics->h_TotalNeighborOffsets[vno-1] + mrics->h_NTotalNeighbors[vno-1];

    for (n=0; n<v->vtotal; n++)
    {
      mrics->h_TotalNeighborArray[ctr++] = v->v[n];
    }
  }

  cudaMemcpy(mrics->d_TotalNeighborArray,
             mrics->h_TotalNeighborArray,
             mrics->ntotalneighbors*sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  MRISCcheckCUDAError("MRISCuploadTotalNeighborArray: "
                      "cannot transfer neighbor array !");

  cudaMemcpy(mrics->d_TotalNeighborOffsets,
             mrics->h_TotalNeighborOffsets,
             mrics->nvertices*sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  MRISCcheckCUDAError("MRISCuploadTotalNeighborArray: "
                      "cannot transfer neighbor offsets !");

  cudaMemcpy(mrics->d_NTotalNeighbors,
             mrics->h_NTotalNeighbors,
             mrics->nvertices*sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  MRISCcheckCUDAError("MRISCuploadTotalNeighborArray: "
                      "cannot transfer ntotalneighbors !");

  cudaMemcpy(mrics->d_NNeighbors,
             mrics->h_NNeighbors,
             mrics->nvertices*sizeof(unsigned int),
             cudaMemcpyHostToDevice);
  MRISCcheckCUDAError("MRISCuploadTotalNeighborArray: "
                      "cannot transfer nneighbors !");
}

void MRISCallocGradients(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{

  if (mrics->nvertices != mris->nvertices)
  {
    mrics->nvertices = mris->nvertices;
    MRISCallocVertices(mrics,mris);
  }

  if (mrics->h_D != NULL)
  {
    free(mrics->h_D);
  }

  mrics->h_D = (float4 *)malloc(mrics->nvertices*sizeof(float4));
  if (mrics->h_D == NULL)
  {
    fprintf(stderr,"MRISCallocGradients: "
            "Can't allocate host memory for gradients !\n");
    exit(-1);
  }

  if (mrics->d_D != NULL)
  {
    cudaFree(mrics->d_D);
  }

  cudaMalloc((void **)&(mrics->d_D),mrics->nvertices*sizeof(float4));
  if (mrics->d_D == NULL)
  {
    fprintf(stderr,"MRISCallocGradients: "
            "Can't allocate device memory for gradients !\n");
    exit(-1);
  }

  if (mrics->d_TD != NULL)
  {
    cudaFree(mrics->d_TD);
  }

  cudaMalloc((void **)&(mrics->d_TD),mrics->nvertices*sizeof(float4));
  if (mrics->d_D == NULL)
  {
    fprintf(stderr,"MRISCallocGradients: "
            "Can't allocate device memory for gradients !\n");
    exit(-1);
  }

}

void MRISCfreeGradients(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  if (mrics->h_D != NULL)
  {
    free(mrics->h_D);
    mrics->h_D = NULL;
  }

  if (mrics->d_D != NULL)
  {
    cudaFree(mrics->d_D);
    mrics->d_D = NULL;
  }

  if (mrics->d_TD != NULL)
  {
    cudaFree(mrics->d_TD);
    mrics->d_TD = NULL;
  }
}

void MRISCuploadGradients(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  int vno;
  float4 *h_D;
#ifdef FS_CUDA_TIMINGS
  struct timeval tv1,tv2,result;
#endif // FS_CUDA_TIMINGS

  /* if it appears that either array is not allocated, reallocate */
  if (mrics->h_D == NULL
      || mrics->d_D == NULL
      || mrics->nvertices != mris->nvertices)
  {
    MRISCallocGradients(mrics,mris);
  }

  h_D = mrics->h_D;

#ifdef FS_CUDA_TIMINGS
  gettimeofday(&tv1,NULL);
#endif // FS_CUDA_TIMINGS

  // copy the data from the vertex structure into the host array
  for (vno=0; vno<mris->nvertices; vno++)
  {
    h_D[vno].x = mris->vertices[vno].dx;
    h_D[vno].y = mris->vertices[vno].dy;
    h_D[vno].z = mris->vertices[vno].dz;
    h_D[vno].w = 0;
  }

  // finally transfer the data
  cudaMemcpy(mrics->d_D,
             mrics->h_D,mris->nvertices*sizeof(float4),
             cudaMemcpyHostToDevice);
  MRISCcheckCUDAError("MRISuploadGradients: cannot transfer data !");

#ifdef FS_CUDA_TIMINGS
  gettimeofday(&tv2,NULL);
  timeval_subtract(&result,&tv2,&tv1);
  printf("MRISuploadGradients: %ld ms\n",
         result.tv_sec*1000+result.tv_usec/1000);
  fflush(stdout);
#endif // FS_CUDA_TIMINGS
}

void MRISCdownloadGradients(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  int vno;
  float4 *h_D;

  /* if it appears that either array is not allocated, reallocate */
  if (mrics->h_D == NULL
      || mrics->d_D == NULL
      || mrics->nvertices != mris->nvertices)
  {
    MRISCallocGradients(mrics,mris);
  }

  h_D = mrics->h_D;

  cudaMemcpy(h_D,mrics->d_D,
             mris->nvertices*sizeof(float4),
             cudaMemcpyDeviceToHost);
  MRISCcheckCUDAError("MRISdownloadGradients: cannot transfer data !");

  // copy the data from the host array into the vertex structure
  for (vno=0; vno<mris->nvertices; vno++)
  {
    mris->vertices[vno].dx = h_D[vno].x;
    mris->vertices[vno].dy = h_D[vno].y;
    mris->vertices[vno].dz = h_D[vno].z;
  }
}

/* Computational Routines */

#define BLOCK_SIZE_CVD 128

__global__ void computeSphereVertexDistancesKernel(float4 *V, float *dist,
    unsigned int *NEIGHBOR,
    unsigned int *NBOFFSETS,
    unsigned int *nNeighbors,
    unsigned int nVertices,
    float circumference)
{
  int n,N;
  int offset,soffset;

  // since we are using multiple threads per blocks as well as multiple blocks
  int vidxb = 4*(blockIdx.x * blockDim.x) + threadIdx.x;
  int basevert = 4*(blockIdx.x * blockDim.x);

  int vidx,tab;
  float4 nv,tv;
  float dot,n1,n2,norm;

  // create a cache for 4 elements per block (4*BLOCK_SIZE elements)
  __shared__ float4 SI[4*BLOCK_SIZE_CVD];

  int bidx = threadIdx.x;
  // this means we have 128 neighboring vertices cached
  for (vidx=vidxb; vidx<vidxb+4*BLOCK_SIZE_CVD; vidx+=BLOCK_SIZE_CVD)
  {
    if (vidx < nVertices)
    {
      SI[bidx] = V[vidx];
      bidx+=BLOCK_SIZE_CVD;
    }
  }

  __syncthreads();

  bidx = threadIdx.x;
  // preload the current BLOCK_SIZE vertices
  for (vidx=vidxb; vidx<vidxb+4*BLOCK_SIZE_CVD; vidx+=BLOCK_SIZE_CVD)
  {
    if (vidx < nVertices)
    {
      offset = NBOFFSETS[ vidx ];
      N = nNeighbors[ vidx ];
      tv = SI[bidx];

      bidx += BLOCK_SIZE_CVD;

      for (n = 0; n < N; n++)
      {
        soffset = NEIGHBOR[offset+n];

        /* There seems to be little to NO benefit of this local caching,
          either because we have no hits, or reading from the shared memory
          is just as slow as reading from global memory
        */
        tab = soffset - basevert;
        if (tab > 0 && tab < 4*BLOCK_SIZE_CVD)
        {
          nv = SI[tab];
        }
        else
        {
          nv = V[soffset];
        }

        // avoid FMADS
        //dot = tv.x*nv.x + tv.y*nv.y + tv.z*nv.z;

        dot = __fmul_rn(tv.x,nv.x);
        dot = __fadd_rn(dot,__fmul_rn(tv.y,nv.y));
        dot = __fadd_rn(dot,__fmul_rn(tv.z,nv.z));

        //n1 = tv.x*tv.x + tv.y*tv.y + tv.z*tv.z;

        n1 = __fmul_rn(tv.x,tv.x);
        n1 = __fadd_rn(n1,__fmul_rn(tv.y,tv.y));
        n1 = __fadd_rn(n1,__fmul_rn(tv.z,tv.z));

        //n2 = nv.x*nv.x + nv.y*nv.y + nv.z*nv.z;

        n2 = __fmul_rn(nv.x,nv.x);
        n2 = __fadd_rn(n2,__fmul_rn(nv.y,nv.y));
        n2 = __fadd_rn(n2,__fmul_rn(nv.z,nv.z));

        norm = __fmul_rn(__fsqrt_rn(n1),__fsqrt_rn(n2));

        //norm = __fsqrt_rn(n1) * __fsqrt_rn(n2);

        // this seems to be a quell of numerical error here
        if (norm < 1.0e-7f)
        {
          dist[offset+n] = 0.0f;
        }
        else if (fabsf(dot) > norm)
        {
          dist[offset+n] = 0.0f;
        }
        else
        {
          dist[offset+n] = __fmul_rn(circumference,fabsf(acosf(dot/norm)));
        }
      }
    }
  }
}

void MRISCcomputeVertexDistances(MRI_CUDA_SURFACE *mrics, MRI_SURFACE *mris)
{
  float circumference;
  float4 v = mrics->h_V[0];

  circumference = sqrtf(v.x*v.x + v.y*v.y + v.z*v.z);
  dim3 dimBlock(BLOCK_SIZE_CVD,1);
  dim3 dimGrid(mrics->nvertices/(4*BLOCK_SIZE_CVD)+1,1);

  computeSphereVertexDistancesKernel<<<dimGrid, dimBlock>>>(mrics->d_V,
      mrics->d_Distances,
      mrics->d_TotalNeighborArray,
      mrics->d_TotalNeighborOffsets,
      mrics->d_NTotalNeighbors,
      mrics->nvertices,
      circumference);
  MRISCcheckCUDAError("MRISCcomputeVertexDistances");
  cudaThreadSynchronize();
}

#define BLOCK_SIZE_AVGG 128

__global__ void GradientAverageKernel(float4 *D,
                                      float4 *TD,
                                      unsigned int *NEIGHBOR,
                                      unsigned int *NBOFFSETS,
                                      unsigned int *nNeighbors,
                                      unsigned int nVertices)
{
  int n,N;
  int offset,soffset;

  // since we are using multiple threads per blocks as well as multiple blocks
  int vidxb = 4*(blockIdx.x * blockDim.x) + threadIdx.x;
  //int basevert = 4*(blockIdx.x * blockDim.x);

  int vidx; //,tab;
  float4 nbd,td;

  // create a cache for 4 elements per block (4*BLOCK_SIZE elements)
  __shared__ float4 SI[4*BLOCK_SIZE_AVGG];

  int bidx = 4*threadIdx.x;
  // this means we have 128 neighboring vertices cached
  for (vidx=vidxb; vidx<vidxb+4*BLOCK_SIZE_AVGG; vidx+=BLOCK_SIZE_AVGG)
  {
    if (vidx < nVertices)
    {
      SI[bidx] = D[vidx];
      bidx++;
    }
  }

  __syncthreads();

  bidx = 4*threadIdx.x;
  // preload the current BLOCK_SIZE vertices
  for (vidx=vidxb; vidx<vidxb+4*BLOCK_SIZE_AVGG; vidx+=BLOCK_SIZE_AVGG)
  {
    if (vidx < nVertices)
    {

      offset = NBOFFSETS[ vidx ];
      N = nNeighbors[ vidx ];

      td = SI[bidx++];

      for (n = 0; n < N; n++)
      {
        soffset = NEIGHBOR[offset+n];
        /*
         tab = soffset - basevert;
        if(tab > 0 && tab < 4*BLOCK_SIZE)
        nbd = SI[tab];
         else
        */
        nbd = D[soffset];

        td.x += nbd.x;
        td.y += nbd.y;
        td.z += nbd.z;
      }

      td.x /= (float)(N+1);
      td.y /= (float)(N+1);
      td.z /= (float)(N+1);

      TD[vidx] = td;
    }
  }
}

__global__ void updateGradientsKernel(float4 *D,
                                      float4 *TD,
                                      unsigned int nVertices)
{
  int vidx = 4*(blockIdx.x * blockDim.x) + threadIdx.x;
  int idx;
  for (idx=0; idx<4*BLOCK_SIZE_AVGG; idx+=BLOCK_SIZE_AVGG)
  {
    D[vidx+idx] = TD[vidx+idx];
  }
}

void MRISCaverageGradients(MRI_CUDA_SURFACE *mrisc,
                           MRI_SURFACE *mris,
                           unsigned int niter)
{
  unsigned int idx;

  dim3 dimBlock(BLOCK_SIZE_AVGG,1);
  dim3 dimGrid(mrisc->nvertices/(4*BLOCK_SIZE_AVGG)+1,1);

  for (idx=0; idx<niter; idx++)
  {
    // run the kernel
    GradientAverageKernel<<<dimGrid, dimBlock>>>(mrisc->d_D,
        mrisc->d_TD,
        mrisc->d_TotalNeighborArray,
        mrisc->d_TotalNeighborOffsets,
        mrisc->d_NNeighbors,
        mrisc->nvertices);
    cudaThreadSynchronize();
    // run the update kernel
    updateGradientsKernel<<<dimGrid, dimBlock>>>(mrisc->d_D,
        mrisc->d_TD,
        mrisc->nvertices);
    cudaThreadSynchronize();
  }
  MRISCcheckCUDAError("MRISCaverageGradients");
}

/* MRISCcomputeDistanceError()

   this function is based on the mrisComputeDistanceError
   function used in mrisurf.c
   Based on observations this function will have the following
   limitations compared to
   the original function
   it does NOT support parms->dist_error NOR parms->vsmoothness
   it also IGNORES the .neg property

   This is a reduction function. It can only be done in
   double precision on SM_13 hardware
*/
template<typename T>
struct square
{
  __host__ __device__
  T operator()(const T& x) const
  {
    return x * x;
  }
};

struct ssxpy_functor
{
  const float a;

  ssxpy_functor(float _a) : a(_a) {}

  __host__ __device__
  float operator()(const float& x, const float& y) const
  {
    return x - a * y;
  }
};

float MRISCcomputeDistanceError(MRI_CUDA_SURFACE *mrisc, float dist_scale)
{
  // allocate a temporary difference vector
  thrust::device_vector<float> Diff(mrisc->ntotalneighbors);

  // generate thrust device pointers from the device pointers we already have
  thrust::device_ptr<float> dist_ptr(mrisc->d_Distances);
  thrust::device_ptr<float> dist_origptr(mrisc->d_DistancesOrig);

  if (dist_scale == 1.0)
  {
    thrust::transform(dist_origptr,
                      dist_origptr+mrisc->ntotalneighbors,
                      dist_ptr,
                      Diff.begin(),
                      thrust::minus<float>());
  }
  else
  {
    thrust::transform(dist_origptr,
                      dist_origptr+mrisc->ntotalneighbors,
                      dist_ptr,
                      Diff.begin(),
                      ssxpy_functor(dist_scale));
  }
  square<float>        unary_op;
  thrust::plus<float> binary_op;
  float init = 0;

  float norm = thrust::transform_reduce(Diff.begin(),
                                        Diff.end(),
                                        unary_op,
                                        init,
                                        binary_op);

  return norm;
}

#endif /* FS_CUDA */
