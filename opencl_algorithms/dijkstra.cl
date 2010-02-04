//
//
//  Description:
//      Implementation of Dijkstra's Single-Source Shortest Path (SSSP) algorithm on the GPU.
//      The basis of this implementation is the paper:
//
//          \"Accelerating large graph algorithms on the GPU using CUDA\" by
//          Parwan Harish and P.J. Narayanan
//
//
//  Author:
//      Dan Ginsburg
//      <daniel.ginsburg@childrens.harvard.edu>
//
//  Children's Hospital Boston
//


///
/// This is part 1 of the Kernel from Algorithm 4 in the paper
///
__kernel  void OCL_SSSP_KERNEL1(__global int *vertexArray, __global int *edgeArray, __global float *weightArray,
                               __global int *maskArray, __global float *costArray, __global float *updatingCostArray,
                               int vertexCount, int edgeCount )
{
    // access thread id
    int tid = get_global_id(0);

    if ( maskArray[tid] != 0 )
    {
        maskArray[tid] = 0;

        int edgeStart = vertexArray[tid];
        int edgeEnd;
        if (tid + 1 < (vertexCount))
        {
            edgeEnd = vertexArray[tid + 1];
        }
        else
        {
            edgeEnd = edgeCount;
        }

        for(int edge = edgeStart; edge < edgeEnd; edge++)
        {
            int nid = edgeArray[edge];

            // One note here: whereas the paper specified weightArray[nid], I
            //  found that the correct thing to do was weightArray[edge].  I think
            //  this was a typo in the paper.  Either that, or I misunderstood
            //  the data structure.
            if (updatingCostArray[nid] > (costArray[tid] + weightArray[edge]))
            {
                updatingCostArray[nid] = (costArray[tid] + weightArray[edge]);
            }
        }
    }
}

///
/// This is part 2 of the Kernel from Algorithm 5 in the paper.  
///
__kernel  void OCL_SSSP_KERNEL2(__global int *vertexArray, __global int *edgeArray, __global float *weightArray,
                                __global int *maskArray, __global float *costArray, __global float *updatingCostArray,
                                int vertexCount)
{
    // access thread id
    int tid = get_global_id(0);


    if (costArray[tid] > updatingCostArray[tid])
    {
        costArray[tid] = updatingCostArray[tid];
        maskArray[tid] = 1;
    }

    updatingCostArray[tid] = costArray[tid];
}


///
/// Kernel to initialize buffers
///
__kernel void initializeBuffers( __global int *maskArray, __global float *costArray, __global float *updatingCostArray,
                                 int sourceVertex, int vertexCount )
{
    // access thread id
    int tid = get_global_id(0);


    if (sourceVertex == tid)
    {
        maskArray[tid] = 1;
        costArray[tid] = 0.0;
        updatingCostArray[tid] = 0.0;
    }
    else
    {
        maskArray[tid] = 0;
        costArray[tid] = FLT_MAX;
        updatingCostArray[tid] = FLT_MAX;
    }

}

