const char *dijkstraStr = "\
// \n\
// \n\
//  Description: \n\
//      Implementation of Dijkstra's Single-Source Shortest Path (SSSP) algorithm on the GPU. \n\
//      The basis of this implementation is the paper: \n\
// \n\
//          \"Accelerating large graph algorithms on the GPU using CUDA\" by \n\
//          Parwan Harish and P.J. Narayanan \n\
// \n\
// \n\
//  Author: \n\
//      Dan Ginsburg \n\
//      <daniel.ginsburg@childrens.harvard.edu> \n\
// \n\
//  Children's Hospital Boston \n\
// \n\
 \n\
 \n\
/// \n\
/// This is part 1 of the Kernel from Algorithm 4 in the paper \n\
/// \n\
__kernel  void OCL_SSSP_KERNEL1(__global int *vertexArray, __global int *edgeArray, __global float *weightArray, \n\
                               __global int *maskArray, __global float *costArray, __global float *updatingCostArray, \n\
                               int vertexCount, int edgeCount ) \n\
{ \n\
    // access thread id \n\
    int tid = get_global_id(0); \n\
 \n\
    if ( maskArray[tid] != 0 ) \n\
    { \n\
        maskArray[tid] = 0; \n\
 \n\
        int edgeStart = vertexArray[tid]; \n\
        int edgeEnd; \n\
        if (tid + 1 < (vertexCount)) \n\
        { \n\
            edgeEnd = vertexArray[tid + 1]; \n\
        } \n\
        else \n\
        { \n\
            edgeEnd = edgeCount; \n\
        } \n\
 \n\
        for(int edge = edgeStart; edge < edgeEnd; edge++) \n\
        { \n\
            int nid = edgeArray[edge]; \n\
 \n\
            // One note here: whereas the paper specified weightArray[nid], I \n\
            //  found that the correct thing to do was weightArray[edge].  I think \n\
            //  this was a typo in the paper.  Either that, or I misunderstood \n\
            //  the data structure. \n\
            if (updatingCostArray[nid] > (costArray[tid] + weightArray[edge])) \n\
            { \n\
                updatingCostArray[nid] = (costArray[tid] + weightArray[edge]); \n\
            } \n\
        } \n\
    } \n\
} \n\
 \n\
/// \n\
/// This is part 2 of the Kernel from Algorithm 5 in the paper. \n\
/// \n\
__kernel  void OCL_SSSP_KERNEL2(__global int *vertexArray, __global int *edgeArray, __global float *weightArray, \n\
                                __global int *maskArray, __global float *costArray, __global float *updatingCostArray, \n\
                                int vertexCount) \n\
{ \n\
    // access thread id \n\
    int tid = get_global_id(0); \n\
 \n\
 \n\
    if (costArray[tid] > updatingCostArray[tid]) \n\
    { \n\
        costArray[tid] = updatingCostArray[tid]; \n\
        maskArray[tid] = 1; \n\
    } \n\
 \n\
    updatingCostArray[tid] = costArray[tid]; \n\
} \n\
 \n\
 \n\
/// \n\
/// Kernel to initialize buffers \n\
/// \n\
__kernel void initializeBuffers( __global int *maskArray, __global float *costArray, __global float *updatingCostArray, \n\
                                 int sourceVertex, int vertexCount ) \n\
{ \n\
    // access thread id \n\
    int tid = get_global_id(0); \n\
 \n\
 \n\
    if (sourceVertex == tid) \n\
    { \n\
        maskArray[tid] = 1; \n\
        costArray[tid] = 0.0; \n\
        updatingCostArray[tid] = 0.0; \n\
    } \n\
    else \n\
    { \n\
        maskArray[tid] = 0; \n\
        costArray[tid] = FLT_MAX; \n\
        updatingCostArray[tid] = FLT_MAX; \n\
    } \n\
 \n\
} \n\
 \n\
";
