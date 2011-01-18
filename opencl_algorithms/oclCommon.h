//
//
//  Description:
//      This contains a set of utility functions used for common interfacing
//      with OpenCL such as error detection and kernel loading.  Much of this
//      code was adapted from the NVIDIA GPU Computing SDK.
//
//  Author:
//      Dan Ginsburg
//      <daniel.ginsburg@childrens.harvard.edu>
//
//  Children's Hospital Boston
//
#ifndef OCL_COMMON_H
#define OCL_COMMON_H

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#ifdef __APPLE__
    #include <OpenCL/cl.h>
#else
    #include <CL/cl.h>
#endif


///
//  Macros
//
#define oclCheckError(a, b) oclCheckErrorFileLine(a, b, __FILE__ , __LINE__)

///
//  Public Functions
//

///
/// Print info about the device to stdout (modified from NVIDIA SDK)
/// @param device         OpenCL id of the device
///
void oclPrintDevInfo(cl_device_id device);

///
/// Loads a Program file and prepends the cPreamble to the code. (from the NVIDIA SDK)
///
/// @return the source string if succeeded, 0 otherwise
/// @param cFilename        program filename
/// @param cPreamble        code that is prepended to the loaded file, typically a set of #defines or a header
/// @param szFinalLength    returned length of the code string
///
char* oclLoadProgSource(const char* cFilename, const char* cPreamble, size_t* szFinalLength);

///
/// Gets the id of the nth device from the context (from the NVIDIA SDK)
///
/// @return the id or -1 when out of range
/// @param cxGPUContext         OpenCL context
/// @param device_idx            index of the device of interest
///
cl_device_id oclGetDev(cl_context cxGPUContext, unsigned int nr);

///
/// Gets the id of the first device from the context (from the NVIDIA SDK)
///
/// @return the id
/// @param cxGPUContext         OpenCL context
///
cl_device_id oclGetFirstDev(cl_context cxGPUContext);

///
/// Gets the id of device with maximal FLOPS from the context (from NVIDIA SDK)
///
/// @return the id
/// @param cxGPUContext         OpenCL context
///
cl_device_id oclGetMaxFlopsDev(cl_context cxGPUContext);

///
/// Check for error condition and exit if found.  Print file and line number
/// of error.
///
/// @param errNum Error value to check
/// @param expected Expected value
/// @param file File name (filled in by macro)
/// @param lineNumber Line number (filled in by macro)
///
void oclCheckErrorFileLine(int errNum, int expected, const char* file, const int lineNumber);

///
/// Round the local work size up to the next multiple of the size
/// @param groupSize Size of the group
/// @param globalSize Global size
///
int oclRoundWorkSizeUp(int groupSize, int globalSize);

#endif // OCL_COMMON_H
