//
//
//  Description:
//      Implementation of Dijkstra's Single-Source Shortest Path (SSSP) algorithm on the GPU.
//      The basis of this implementation is the paper:
//
//          "Accelerating large graph algorithms on the GPU using CUDA" by
//          Parwan Harish and P.J. Narayanan
//
//
//  Author:
//      Dan Ginsburg
//      <daniel.ginsburg@childrens.harvard.edu>
//
//  Children's Hospital Boston
//
#include <float.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include "oclDijkstraKernel.h"

///
//  Macro Options
//
#define KERNEL_AS_HEADER

///
//  Function prototypes
//
bool maskArrayEmpty(int *maskArray, int count);

///
//  Namespaces
//
using namespace std;

#ifdef FS_OPENCL
#include "oclCommon.h"
#include <pthread.h>

#ifdef KERNEL_AS_HEADER
    #include "dijkstra.cl.h"
#endif


///
//  Types
//

// This structure is used in the multi-GPU implementation of the algorithm.
// This structure defines the workload for each GPU.  The code chunks up
// the work on a per-GPU basis.
typedef struct
{
    // Context
    cl_context context;

    // Device number to run algorithm on
    cl_device_id deviceId;

    // Pointer to graph data
    GraphData *graph;

    // Source vertex indices to process
    int *sourceVertices;

    // End vertex indices to process
    int *endVertices;

    // Results of processing
    float *outResultCosts;

    // Number of results
    int numResults;

} DevicePlan;

///
//  Globals
//
pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

///////////////////////////////////////////////////////////////////////////////
//
//  Private Functions
//
//

///
/// Load and build an OpenCL program from source file
/// \param gpuContext GPU context on which to load and build the program
/// \param fileName File name of source file that holds the kernels
/// \return Handle to the program
///
cl_program loadAndBuildProgram( cl_context gpuContext, const char *fileName )
{
    pthread_mutex_lock(&mutex);

    size_t programLength;
    cl_int errNum;
    cl_program program;

    // Load the OpenCL source code from the .cl file
#ifndef KERNEL_AS_HEADER
    const char* sourcePath = "dijkstra.cl";
    const char *source = oclLoadProgSource(sourcePath, "", &programLength);
#else
    const char *source = dijkstraStr;
    programLength = strlen(source) + 1;
#endif

    oclCheckError(source != NULL, true);

    // Create the program for all GPUs in the context
    program = clCreateProgramWithSource(gpuContext, 1, (const char **)&source, &programLength, &errNum);
    oclCheckError(errNum, CL_SUCCESS);
    // build the program for all devices on the context
    errNum = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (errNum != CL_SUCCESS)
    {
        char cBuildLog[10240];
        clGetProgramBuildInfo(program, oclGetFirstDev(gpuContext), CL_PROGRAM_BUILD_LOG,
                              sizeof(cBuildLog), cBuildLog, NULL );

        cerr << cBuildLog << endl;
        oclCheckError(errNum, CL_SUCCESS);
    }

    pthread_mutex_unlock(&mutex);
    return program;
}

///
///  Allocate memory for input CUDA buffers and copy the data into device memory
///
void allocateOCLBuffers(cl_context gpuContext, cl_command_queue commandQueue, GraphData *graph,
                        cl_mem *vertexArrayDevice, cl_mem *edgeArrayDevice, cl_mem *weightArrayDevice,
                        cl_mem *maskArrayDevice, cl_mem *costArrayDevice, cl_mem *updatingCostArrayDevice,
                        size_t globalWorkSize)
{
    cl_int errNum;
    cl_mem hostVertexArrayBuffer;
    cl_mem hostEdgeArrayBuffer;
    cl_mem hostWeightArrayBuffer;

    // First, need to create OpenCL Host buffers that can be copied to device buffers
    hostVertexArrayBuffer = clCreateBuffer(gpuContext, CL_MEM_COPY_HOST_PTR | CL_MEM_ALLOC_HOST_PTR,
                                           sizeof(int) * graph->vertexCount, graph->vertexArray, &errNum);
    oclCheckError(errNum, CL_SUCCESS);

    hostEdgeArrayBuffer = clCreateBuffer(gpuContext, CL_MEM_COPY_HOST_PTR | CL_MEM_ALLOC_HOST_PTR,
                                           sizeof(int) * graph->edgeCount, graph->edgeArray, &errNum);
    oclCheckError(errNum, CL_SUCCESS);

    hostWeightArrayBuffer = clCreateBuffer(gpuContext, CL_MEM_COPY_HOST_PTR | CL_MEM_ALLOC_HOST_PTR,
                                           sizeof(float) * graph->edgeCount, graph->weightArray, &errNum);
    oclCheckError(errNum, CL_SUCCESS);

    // Now create all of the GPU buffers
    *vertexArrayDevice = clCreateBuffer(gpuContext, CL_MEM_READ_ONLY, sizeof(int) * globalWorkSize, NULL, &errNum);
    oclCheckError(errNum, CL_SUCCESS);
    *edgeArrayDevice = clCreateBuffer(gpuContext, CL_MEM_READ_ONLY, sizeof(int) * graph->edgeCount, NULL, &errNum);
    oclCheckError(errNum, CL_SUCCESS);
    *weightArrayDevice = clCreateBuffer(gpuContext, CL_MEM_READ_ONLY, sizeof(float) * graph->edgeCount, NULL, &errNum);
    oclCheckError(errNum, CL_SUCCESS);
    *maskArrayDevice = clCreateBuffer(gpuContext, CL_MEM_READ_WRITE, sizeof(int) * globalWorkSize, NULL, &errNum);
    oclCheckError(errNum, CL_SUCCESS);
    *costArrayDevice = clCreateBuffer(gpuContext, CL_MEM_READ_WRITE, sizeof(float) * globalWorkSize, NULL, &errNum);
    oclCheckError(errNum, CL_SUCCESS);
    *updatingCostArrayDevice = clCreateBuffer(gpuContext, CL_MEM_READ_WRITE, sizeof(float) * globalWorkSize, NULL, &errNum);
    oclCheckError(errNum, CL_SUCCESS);

    // Now queue up the data to be copied to the device
    errNum = clEnqueueCopyBuffer(commandQueue, hostVertexArrayBuffer, *vertexArrayDevice, 0, 0,
                                 sizeof(int) * graph->vertexCount, 0, NULL, NULL);
    oclCheckError(errNum, CL_SUCCESS);

    errNum = clEnqueueCopyBuffer(commandQueue, hostEdgeArrayBuffer, *edgeArrayDevice, 0, 0,
                                 sizeof(int) * graph->edgeCount, 0, NULL, NULL);
    oclCheckError(errNum, CL_SUCCESS);

    errNum = clEnqueueCopyBuffer(commandQueue, hostWeightArrayBuffer, *weightArrayDevice, 0, 0,
                                 sizeof(float) * graph->edgeCount, 0, NULL, NULL);
    oclCheckError(errNum, CL_SUCCESS);

    clReleaseMemObject(hostVertexArrayBuffer);
    clReleaseMemObject(hostEdgeArrayBuffer);
    clReleaseMemObject(hostWeightArrayBuffer);
}

///
/// Initialize OpenCL buffers for single run of Dijkstra
///
void initializeOCLBuffers(cl_command_queue commandQueue, cl_kernel initializeKernel, GraphData *graph,
                          size_t maxWorkGroupSize)
{
    cl_int errNum;
    // Set # of work items in work group and total in 1 dimensional range
    size_t localWorkSize = maxWorkGroupSize;
    size_t globalWorkSize = oclRoundWorkSizeUp(localWorkSize, graph->vertexCount);

    errNum = clEnqueueNDRangeKernel(commandQueue, initializeKernel, 1, NULL, &globalWorkSize, &localWorkSize,
                                    0, NULL, NULL);
    oclCheckError(errNum, CL_SUCCESS);
}

///
/// Worker thread for running the algorithm on one of the compute devices
///
void dijkstraThread(DevicePlan *plan)
{
    runDijkstra( plan->context, plan->deviceId, plan->graph, plan->sourceVertices,
                 plan->outResultCosts, plan->numResults );
}

///////////////////////////////////////////////////////////////////////////////
//
//  Public Functions
//
//

///
/// Run Dijkstra's shortest path on the GraphData provided to this function.  This
/// function will compute the shortest path distance from sourceVertices[n] ->
/// endVertices[n] and store the cost in outResultCosts[n].  The number of results
/// it will compute is given by numResults.
///
/// This version of the function will run the algorithm on either just the CPU,
/// CPU + GPU, GPU, or Multi GPU depending on what compute resources are available
/// on the system.
///
/// \param graph Structure containing the vertex, edge, and weight arra
///              for the input graph
/// \param startVertices Indices into the vertex array from which to
///                      start the search
/// \param outResultsCosts A pre-allocated array where the results for
///                        each shortest path search will be written.
///                        This must be sized numResults * graph->numVertices.
/// \param numResults Should be the size of all three passed inarrays
///
void runDijkstraOpenCL( GraphData* graph, int *sourceVertices,
                        float *outResultCosts, int numResults )
{
    // See what kind of devices are available
    cl_int errNum;
    cl_context cpuContext;
    cl_context gpuContext;

    // create the OpenCL context on available GPU devices
    gpuContext = clCreateContextFromType(0, CL_DEVICE_TYPE_GPU, NULL, NULL, &errNum);

    // Create an OpenCL context on available CPU devices
    cpuContext = clCreateContextFromType(0, CL_DEVICE_TYPE_CPU, NULL, NULL, &errNum);

    if (cpuContext == 0 && gpuContext == 0)
    {
        cerr << "ERROR: could not create any OpenCL context on CPU or GPU" << endl;
        return;
    }

    // For just a single result, just use multi-threaded CPU or single GPU
    if (numResults == 1)
    {
        if (gpuContext != 0)
        {
            cout << "Dijkstra OpenCL: Running single GPU version." << endl;
            runDijkstra(gpuContext, oclGetMaxFlopsDev(gpuContext), graph, sourceVertices,
                        outResultCosts, numResults);
        }
        else
        {
            cout << "Dijkstra OpenCL: Running multithreaded CPU version." << endl;
            runDijkstra(cpuContext, oclGetMaxFlopsDev(cpuContext), graph, sourceVertices,
                        outResultCosts, numResults);
        }
    }
    // For multiple results, prefer multi-GPU and fallback to CPU
    else
    {
        // Prefer Multi-GPU if multiple GPUs are available
        if (gpuContext != 0)
        {
            cout << "Dijkstra OpenCL: Running multi-GPU version." << endl;
            runDijkstraMultiGPU( gpuContext, graph, sourceVertices,
                                 outResultCosts, numResults );
        }
        // For now, fallback to CPU in this case.  I have a multi GPU+CPU path
        // but it does not seem to perform well because of the CPU overhead of
        // running the GPU version slows down the CPU version.
        else
        {
            cout << "Dijkstra OpenCL: Running multithreaded CPU version." << endl;
            runDijkstra(cpuContext, oclGetMaxFlopsDev(cpuContext), graph, sourceVertices,
                        outResultCosts, numResults);
        }
    }

    clReleaseContext(cpuContext);
    clReleaseContext(gpuContext);
}

///
/// Run Dijkstra's shortest path on the GraphData provided to this function.  This
/// function will compute the shortest path distance from sourceVertices[n] ->
/// endVertices[n] and store the cost in outResultCosts[n].  The number of results
/// it will compute is given by numResults.
///
/// This function will run the algorithm on a single GPU.
///
/// \param gpuContext Current GPU context, must be created by caller
/// \param deviceId The device ID on which to run the kernel.  This can
///                 be determined externally by the caller or the multi
///                 GPU version will automatically split the work across
///                 devices
/// \param graph Structure containing the vertex, edge, and weight arra
///              for the input graph
/// \param startVertices Indices into the vertex array from which to
///                      start the search
/// \param outResultsCosts A pre-allocated array where the results for
///                        each shortest path search will be written
/// \param numResults Should be the size of all three passed inarrays
///
void runDijkstra( cl_context gpuContext, cl_device_id deviceId, GraphData* graph,
                  int *sourceVertices, float *outResultCosts, int numResults)
{
    // Create command queue
    cl_int errNum;
    cl_command_queue commandQueue;
    commandQueue = clCreateCommandQueue( gpuContext, deviceId, 0, &errNum );
    oclCheckError(errNum, CL_SUCCESS);

    // Program handle
    cl_program program = loadAndBuildProgram( gpuContext, "dijkstra.cl" );
    if (program <= 0 )
    {
        return;
    }

    // Get the max workgroup size
    size_t maxWorkGroupSize;
    clGetDeviceInfo(deviceId, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &maxWorkGroupSize, NULL);
    oclCheckError(errNum, CL_SUCCESS);
    cout << "MAX_WORKGROUP_SIZE: " << maxWorkGroupSize << endl;

    // Set # of work items in work group and total in 1 dimensional range
    size_t localWorkSize = maxWorkGroupSize;
    size_t globalWorkSize = oclRoundWorkSizeUp(localWorkSize, graph->vertexCount);

    cl_mem vertexArrayDevice;
    cl_mem edgeArrayDevice;
    cl_mem weightArrayDevice;
    cl_mem maskArrayDevice;
    cl_mem costArrayDevice;
    cl_mem updatingCostArrayDevice;

    // Allocate buffers in Device memory
    allocateOCLBuffers( gpuContext, commandQueue, graph, &vertexArrayDevice, &edgeArrayDevice, &weightArrayDevice,
                        &maskArrayDevice, &costArrayDevice, &updatingCostArrayDevice, globalWorkSize);


    // Create the Kernels
    cl_kernel initializeBuffersKernel;
    initializeBuffersKernel = clCreateKernel(program, "initializeBuffers", &errNum);
    oclCheckError(errNum, CL_SUCCESS);

    // Set the args values and check for errors
    errNum |= clSetKernelArg(initializeBuffersKernel, 0, sizeof(cl_mem), &maskArrayDevice);
    errNum |= clSetKernelArg(initializeBuffersKernel, 1, sizeof(cl_mem), &costArrayDevice);
    errNum |= clSetKernelArg(initializeBuffersKernel, 2, sizeof(cl_mem), &updatingCostArrayDevice);

    // 3 set below in loop
    errNum |= clSetKernelArg(initializeBuffersKernel, 4, sizeof(int), &graph->vertexCount);
    oclCheckError(errNum, CL_SUCCESS);

    // Kernel 1
    cl_kernel ssspKernel1;
    ssspKernel1 = clCreateKernel(program, "OCL_SSSP_KERNEL1", &errNum);
    oclCheckError(errNum, CL_SUCCESS);
    errNum |= clSetKernelArg(ssspKernel1, 0, sizeof(cl_mem), &vertexArrayDevice);
    errNum |= clSetKernelArg(ssspKernel1, 1, sizeof(cl_mem), &edgeArrayDevice);
    errNum |= clSetKernelArg(ssspKernel1, 2, sizeof(cl_mem), &weightArrayDevice);
    errNum |= clSetKernelArg(ssspKernel1, 3, sizeof(cl_mem), &maskArrayDevice);
    errNum |= clSetKernelArg(ssspKernel1, 4, sizeof(cl_mem), &costArrayDevice);
    errNum |= clSetKernelArg(ssspKernel1, 5, sizeof(cl_mem), &updatingCostArrayDevice);
    errNum |= clSetKernelArg(ssspKernel1, 6, sizeof(int), &graph->vertexCount);
    errNum |= clSetKernelArg(ssspKernel1, 7, sizeof(int), &graph->edgeCount);
    oclCheckError(errNum, CL_SUCCESS);

    // Kernel 2
    cl_kernel ssspKernel2;
    ssspKernel2 = clCreateKernel(program, "OCL_SSSP_KERNEL2", &errNum);
    oclCheckError(errNum, CL_SUCCESS);
    errNum |= clSetKernelArg(ssspKernel2, 0, sizeof(cl_mem), &vertexArrayDevice);
    errNum |= clSetKernelArg(ssspKernel2, 1, sizeof(cl_mem), &edgeArrayDevice);
    errNum |= clSetKernelArg(ssspKernel2, 2, sizeof(cl_mem), &weightArrayDevice);
    errNum |= clSetKernelArg(ssspKernel2, 3, sizeof(cl_mem), &maskArrayDevice);
    errNum |= clSetKernelArg(ssspKernel2, 4, sizeof(cl_mem), &costArrayDevice);
    errNum |= clSetKernelArg(ssspKernel2, 5, sizeof(cl_mem), &updatingCostArrayDevice);
    errNum |= clSetKernelArg(ssspKernel2, 6, sizeof(int), &graph->vertexCount);

    oclCheckError(errNum, CL_SUCCESS);

    int *maskArrayHost = (int*) malloc(sizeof(int) * graph->vertexCount);

    for ( int i = 0 ; i < numResults; i++ )
    {

        errNum |= clSetKernelArg(initializeBuffersKernel, 3, sizeof(int), &sourceVertices[i]);
        oclCheckError(errNum, CL_SUCCESS);

        // Initialize mask array to false, C and U to infiniti
        initializeOCLBuffers( commandQueue, initializeBuffersKernel, graph, maxWorkGroupSize );

        // Read mask array from device -> host
        cl_event readDone;
        errNum = clEnqueueReadBuffer( commandQueue, maskArrayDevice, CL_FALSE, 0, sizeof(int) * graph->vertexCount,
                                      maskArrayHost, 0, NULL, &readDone);
        oclCheckError(errNum, CL_SUCCESS);
        clWaitForEvents(1, &readDone);

        while(!maskArrayEmpty(maskArrayHost, graph->vertexCount))
        {
            size_t localWorkSize = maxWorkGroupSize;
            size_t globalWorkSize = oclRoundWorkSizeUp(localWorkSize, graph->vertexCount);

            // execute the kernel
            errNum = clEnqueueNDRangeKernel(commandQueue, ssspKernel1, 1, 0, &globalWorkSize, &localWorkSize,
                                           0, NULL, NULL);
            oclCheckError(errNum, CL_SUCCESS);

            errNum = clEnqueueNDRangeKernel(commandQueue, ssspKernel2, 1, 0, &globalWorkSize, &localWorkSize,
                                           0, NULL, NULL);
            oclCheckError(errNum, CL_SUCCESS);

            errNum = clEnqueueReadBuffer(commandQueue, maskArrayDevice, CL_FALSE, 0, sizeof(int) * graph->vertexCount,
                                         maskArrayHost, 0, NULL, &readDone);
            oclCheckError(errNum, CL_SUCCESS);
            clWaitForEvents(1, &readDone);
        }

        // Copy the result back
        errNum = clEnqueueReadBuffer(commandQueue, costArrayDevice, CL_FALSE, 0, sizeof(float) * graph->vertexCount,
                                     &outResultCosts[i * graph->vertexCount], 0, NULL, &readDone);
        oclCheckError(errNum, CL_SUCCESS);
        clWaitForEvents(1, &readDone);
    }

    free (maskArrayHost);

    clReleaseMemObject(vertexArrayDevice);
    clReleaseMemObject(edgeArrayDevice);
    clReleaseMemObject(weightArrayDevice);
    clReleaseMemObject(maskArrayDevice);
    clReleaseMemObject(costArrayDevice);
    clReleaseMemObject(updatingCostArrayDevice);

    clReleaseKernel(initializeBuffersKernel);
    clReleaseKernel(ssspKernel1);
    clReleaseKernel(ssspKernel2);

    clReleaseCommandQueue(commandQueue);
    clReleaseProgram(program);
}



///
/// Run Dijkstra's shortest path on the GraphData provided to this function.  This
/// function will compute the shortest path distance from sourceVertices[n] ->
/// endVertices[n] and store the cost in outResultCosts[n].  The number of results
/// it will compute is given by numResults.
///
/// This function will run the algorithm on as many GPUs as is available.  It will
/// create N threads, one for each GPU, and chunk the workload up to perform
/// (numResults / N) searches per GPU.
///
/// \param gpuContext Current GPU context, must be created by caller
/// \param graph Structure containing the vertex, edge, and weight arra
///              for the input graph
/// \param startVertices Indices into the vertex array from which to
///                      start the search
/// \param endVertices Indices into the vertex array from which to end
///                    the search.
/// \param outResultsCosts A pre-allocated array where the results for
///                        each shortest path search will be written
/// \param numResults Should be the size of all three passed inarrays
///
///
void runDijkstraMultiGPU( cl_context gpuContext, GraphData* graph, int *sourceVertices,
                          float *outResultCosts, int numResults )
{

    // Find out how many GPU's to compute on all available GPUs
    cl_int errNum;
    size_t deviceBytes;
    cl_uint deviceCount;

    errNum = clGetContextInfo(gpuContext, CL_CONTEXT_DEVICES, 0, NULL, &deviceBytes);
    oclCheckError(errNum, CL_SUCCESS);
    deviceCount = (cl_uint)deviceBytes/sizeof(cl_device_id);

    if (deviceCount == 0)
    {
        cerr << "ERROR: no GPUs present!" << endl;
        return;
    }

    DevicePlan *devicePlans = (DevicePlan*) malloc(sizeof(DevicePlan) * deviceCount);
    pthread_t *threadIDs = (pthread_t*) malloc(sizeof(pthread_t) * deviceCount);

    // Divide the workload out per device
    int resultsPerDevice = numResults / deviceCount;

    int offset = 0;

    for (unsigned int i = 0; i < deviceCount; i++)
    {
        devicePlans[i].context = gpuContext;
        devicePlans[i].deviceId = oclGetDev(gpuContext, i);;
        devicePlans[i].graph = graph;
        devicePlans[i].sourceVertices = &sourceVertices[offset];
        devicePlans[i].outResultCosts = &outResultCosts[offset * graph->vertexCount];
        devicePlans[i].numResults = resultsPerDevice;

        oclPrintDevInfo(devicePlans[i].deviceId);
        offset += resultsPerDevice;
    }

    // Add any remaining work to the last GPU
    if (offset < numResults)
    {
        devicePlans[deviceCount - 1].numResults += (numResults - offset);
    }

    // Launch all the threads
    for (unsigned int i = 0; i < deviceCount; i++)
    {
        pthread_create(&threadIDs[i], NULL, (void* (*)(void*))dijkstraThread, (void*)(devicePlans + i));
    }

    // Wait for the results from all threads
    for (unsigned int i = 0; i < deviceCount; i++)
    {
        pthread_join(threadIDs[i], NULL);
    }

    free (devicePlans);
    free (threadIDs);
}

///
/// Run Dijkstra's shortest path on the GraphData provided to this function.  This
/// function will compute the shortest path distance from sourceVertices[n] ->
/// endVertices[n] and store the cost in outResultCosts[n].  The number of results
/// it will compute is given by numResults.
///
/// This function will run the algorithm on as many GPUs as is available along with
/// the CPU.  It will create N threads, one for each device, and chunk the workload up to perform
/// (numResults / N) searches per device.
///
/// \param gpuContext Current GPU context, must be created by caller
/// \param cpuContext Current CPU context, must be created by caller
/// \param graph Structure containing the vertex, edge, and weight arra
///              for the input graph
/// \param startVertices Indices into the vertex array from which to
///                      start the search
/// \param outResultsCosts A pre-allocated array where the results for
///                        each shortest path search will be written
/// \param numResults Should be the size of all three passed inarrays
///
///
void runDijkstraMultiGPUandCPU( cl_context gpuContext, cl_context cpuContext, GraphData* graph,
                                int *sourceVertices,
                                float *outResultCosts, int numResults )
{
    float ratioCPUtoGPU = 0.65; // CPU seems to run it at 2.26X on GT120 GPU

    // Find out how many GPU's to compute on all available GPUs
    cl_int errNum;
    size_t deviceBytes;
    cl_uint gpuDeviceCount;
    cl_uint cpuDeviceCount;

    errNum = clGetContextInfo(gpuContext, CL_CONTEXT_DEVICES, 0, NULL, &deviceBytes);
    oclCheckError(errNum, CL_SUCCESS);
    gpuDeviceCount = (cl_uint)deviceBytes/sizeof(cl_device_id);

    if (gpuDeviceCount == 0)
    {
        cerr << "ERROR: no GPUs present!" << endl;
        return;
    }

    errNum = clGetContextInfo(cpuContext, CL_CONTEXT_DEVICES, 0, NULL, &deviceBytes);
    oclCheckError(errNum, CL_SUCCESS);
    cpuDeviceCount = (cl_uint)deviceBytes/sizeof(cl_device_id);

    if (cpuDeviceCount == 0)
    {
        cerr << "ERROR: no CPUs present!" << endl;
        return;
    }

    cl_uint totalDeviceCount = gpuDeviceCount + cpuDeviceCount;

    DevicePlan *devicePlans = (DevicePlan*) malloc(sizeof(DevicePlan) * totalDeviceCount);
    pthread_t *threadIDs = (pthread_t*) malloc(sizeof(pthread_t) * totalDeviceCount);

    int gpuResults = numResults / (ratioCPUtoGPU);
    int cpuResults = numResults - gpuResults;

    // Divide the workload out per device
    int resultsPerGPU = gpuResults / totalDeviceCount;

    int offset = 0;

    int curDevice = 0;
    for (unsigned int i = 0; i < gpuDeviceCount; i++)
    {
        devicePlans[curDevice].context = gpuContext;
        devicePlans[curDevice].deviceId = oclGetDev(gpuContext, i);;
        devicePlans[curDevice].graph = graph;
        devicePlans[curDevice].sourceVertices = &sourceVertices[offset];
        devicePlans[curDevice].outResultCosts = &outResultCosts[offset * graph->vertexCount];
        devicePlans[curDevice].numResults = resultsPerGPU;

        oclPrintDevInfo(devicePlans[curDevice].deviceId);
        offset += resultsPerGPU;
        curDevice++;
    }

    int resultsPerCPU = cpuResults;

    for (unsigned int i = 0; i < cpuDeviceCount; i++)
    {
        devicePlans[curDevice].context = cpuContext;
        devicePlans[curDevice].deviceId = oclGetDev(cpuContext, i);;
        devicePlans[curDevice].graph = graph;
        devicePlans[curDevice].sourceVertices = &sourceVertices[offset];
        devicePlans[curDevice].outResultCosts = &outResultCosts[offset * graph->vertexCount];
        devicePlans[curDevice].numResults = resultsPerCPU;

        oclPrintDevInfo(devicePlans[curDevice].deviceId);
        offset += resultsPerCPU;
        curDevice++;
    }

    // Add any remaining work to the last GPU
    if (offset < numResults)
    {
        devicePlans[totalDeviceCount - 1].numResults += (numResults - offset);
    }

    // Launch all the threads
    for (unsigned int i = 0; i < totalDeviceCount; i++)
    {
        pthread_create(&threadIDs[i], NULL, (void* (*)(void*))dijkstraThread, (void*)(devicePlans + i));
    }

    // Wait for the results from all threads
    for (unsigned int i = 0; i < totalDeviceCount; i++)
    {
        pthread_join(threadIDs[i], NULL);
    }

    free (devicePlans);
    free (threadIDs);
}

#endif // FS_OPENCL

///
/// Check whether the mask array is empty.  This tells the algorithm whether
/// it needs to continue running or not.
///
bool maskArrayEmpty(int *maskArray, int count)
{
    for(int i = 0; i < count; i++ )
    {
        if (maskArray[i] == 1)
        {
            return false;
        }
    }

    return true;
}

///
/// Run Dijkstra's shortest path on the GraphData provided to this function.  This
/// function will compute the shortest path distance from sourceVertices[n] ->
/// endVertices[n] and store the cost in outResultCosts[n].  The number of results
/// it will compute is given by numResults.
///
/// This is a CPU *REFERENCE* implementation for use as a fallback.
///
/// \param graph Structure containing the vertex, edge, and weight arra
///              for the input graph
/// \param startVertices Indices into the vertex array from which to
///                      start the search
/// \param outResultsCosts A pre-allocated array where the results for
///                        each shortest path search will be written.
///                        This must be sized numResults * graph->numVertices.
/// \param numResults Should be the size of all three passed inarrays
///
void runDijkstraRef( GraphData* graph, int *sourceVertices,
                     float *outResultCosts, int numResults )
{

    // Create the arrays needed for processing the algorithm
    float *costArray = new float[graph->vertexCount];
    float *updatingCostArray = new float[graph->vertexCount];
    int *maskArray = new int[graph->vertexCount];

    for (int i = 0; i < numResults; i++)
    {
        // Initialize the buffer for this run
        for (int v = 0; v < graph->vertexCount; v++)
        {
            if (v == sourceVertices[i])
            {
                maskArray[v] = 1;
                costArray[v] = 0.0;
                updatingCostArray[v] = 0.0;
            }
            else
            {
                maskArray[v] = 0;
                costArray[v] = FLT_MAX;
                updatingCostArray[v] = FLT_MAX;
            }
        }

        while(!maskArrayEmpty(maskArray, graph->vertexCount))
        {
            // Equivalent of OCL_SSSP_KERNEL1()
            for (int tid = 0; tid < graph->vertexCount; tid++)
            {
                if ( maskArray[tid] != 0 )
                {
                    maskArray[tid] = 0;

                    int edgeStart = graph->vertexArray[tid];
                    int edgeEnd;
                    if (tid + 1 < (graph->vertexCount))
                    {
                        edgeEnd = graph->vertexArray[tid + 1];
                    }
                    else
                    {
                        edgeEnd = graph->edgeCount;
                    }

                    for(int edge = edgeStart; edge < edgeEnd; edge++)
                    {
                        int nid = graph->edgeArray[edge];

                        // One note here: whereas the paper specified weightArray[nid], I
                        //  found that the correct thing to do was weightArray[edge].  I think
                        //  this was a typo in the paper.  Either that, or I misunderstood
                        //  the data structure.
                        if (updatingCostArray[nid] > (costArray[tid] + graph->weightArray[edge]))
                        {
                            updatingCostArray[nid] = (costArray[tid] + graph->weightArray[edge]);
                        }
                    }
                }
            }

            // Equivalent of OCL_SSSP_KERNEL2()
            for (int tid = 0; tid < graph->vertexCount; tid++)
            {
                if (costArray[tid] > updatingCostArray[tid])
                {
                    costArray[tid] = updatingCostArray[tid];
                    maskArray[tid] = 1;
                }

                updatingCostArray[tid] = costArray[tid];
            }
        }

        // Copy the result back
        memcpy(&outResultCosts[i * graph->vertexCount], costArray, sizeof(float) * graph->vertexCount);
    }

    // Free temporary computation buffers
    delete [] costArray;
    delete [] updatingCostArray;
    delete [] maskArray;
}
