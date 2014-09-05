//
//
//  Description:
//      This contains a set of utility functions used for common interfacing
//      with OpenCL such as error detection and kernel loading.
//
//
//  Author:
//      Dan Ginsburg
//      <daniel.ginsburg@childrens.harvard.edu>
//
//  Children's Hospital Boston
//
#include "oclCommon.h"
#include <iostream>
#include <stdlib.h>
#include <string.h>

///
//  Namespaces
//
using namespace std;

///
//  Public Functions
//

///
/// Print info about the device to stdout (modified from NVIDIA SDK)
/// @param device         OpenCL id of the device
///
void oclPrintDevInfo(cl_device_id device)
{
    char device_string[1024];
    bool nv_device_attibute_query = false;

    // CL_DEVICE_VENDOR
    clGetDeviceInfo(device, CL_DEVICE_VENDOR, sizeof(device_string), &device_string, NULL);
    fprintf(stdout, "  CL_DEVICE_VENDOR: \t\t\t%s\n", device_string);

    // CL_DEVICE_NAME
    clGetDeviceInfo(device, CL_DEVICE_NAME, sizeof(device_string), &device_string, NULL);
    fprintf(stdout, "  CL_DEVICE_NAME: \t\t\t%s\n", device_string);

    // CL_DRIVER_VERSION
    clGetDeviceInfo(device, CL_DRIVER_VERSION, sizeof(device_string), &device_string, NULL);
    fprintf(stdout, "  CL_DRIVER_VERSION: \t\t\t%s\n", device_string);

    // CL_DEVICE_INFO
    cl_device_type type;
    clGetDeviceInfo(device, CL_DEVICE_TYPE, sizeof(type), &type, NULL);
    if( type & CL_DEVICE_TYPE_CPU )
        fprintf(stdout, "  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_CPU");
    if( type & CL_DEVICE_TYPE_GPU )
        fprintf(stdout, "  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_GPU");
    if( type & CL_DEVICE_TYPE_ACCELERATOR )
        fprintf(stdout, "  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_ACCELERATOR");
    if( type & CL_DEVICE_TYPE_DEFAULT )
        fprintf(stdout, "  CL_DEVICE_TYPE:\t\t\t%s\n", "CL_DEVICE_TYPE_DEFAULT");

    // CL_DEVICE_MAX_COMPUTE_UNITS
    cl_uint compute_units;
    clGetDeviceInfo(device, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(compute_units), &compute_units, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_COMPUTE_UNITS:\t\t%d\n", compute_units);

    // CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS
    cl_uint max_dims;
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS, sizeof(max_dims), &max_dims, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS:\t%d\n", max_dims);

    // CL_DEVICE_MAX_WORK_ITEM_SIZES
    size_t workitem_size[3];
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(workitem_size), &workitem_size, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_WORK_ITEM_SIZES:\t%d / %d / %d \n", (int)workitem_size[0], (int)workitem_size[1], (int)workitem_size[2]);

    // CL_DEVICE_MAX_WORK_GROUP_SIZE
    size_t workgroup_size;
    clGetDeviceInfo(device, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(workgroup_size), &workgroup_size, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_WORK_GROUP_SIZE:\t%d\n", (int)workgroup_size);

    // CL_DEVICE_MAX_CLOCK_FREQUENCY
    cl_uint clock_frequency;
    clGetDeviceInfo(device, CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clock_frequency), &clock_frequency, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_CLOCK_FREQUENCY:\t%u MHz\n", clock_frequency);

    // CL_DEVICE_ADDRESS_BITS
    cl_uint addr_bits;
    clGetDeviceInfo(device, CL_DEVICE_ADDRESS_BITS, sizeof(addr_bits), &addr_bits, NULL);
    fprintf(stdout, "  CL_DEVICE_ADDRESS_BITS:\t\t%d\n", addr_bits);

    // CL_DEVICE_IMAGE_SUPPORT
    cl_bool image_support;
    clGetDeviceInfo(device, CL_DEVICE_IMAGE_SUPPORT, sizeof(image_support), &image_support, NULL);
    fprintf(stdout, "  CL_DEVICE_IMAGE_SUPPORT:\t\t%d\n", image_support);

    // CL_DEVICE_MAX_READ_IMAGE_ARGS
    cl_uint max_read_image_args;
    clGetDeviceInfo(device, CL_DEVICE_MAX_READ_IMAGE_ARGS, sizeof(max_read_image_args), &max_read_image_args, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_READ_IMAGE_ARGS:\t%d\n", max_read_image_args);

    // CL_DEVICE_MAX_WRITE_IMAGE_ARGS
    cl_uint max_write_image_args;
    clGetDeviceInfo(device, CL_DEVICE_MAX_WRITE_IMAGE_ARGS, sizeof(max_write_image_args), &max_write_image_args, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_WRITE_IMAGE_ARGS:\t%d\n", max_write_image_args);

    // CL_DEVICE_IMAGE2D_MAX_WIDTH
    size_t image2d_max_width;
    clGetDeviceInfo(device, CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(image2d_max_width), &image2d_max_width, NULL);
    size_t image2d_max_height;
    clGetDeviceInfo(device, CL_DEVICE_IMAGE2D_MAX_HEIGHT, sizeof(image2d_max_height), &image2d_max_height, NULL);
    size_t image3d_max_width;
    clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_WIDTH, sizeof(image3d_max_width), &image3d_max_width, NULL);
    size_t image3d_max_height;
    clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_HEIGHT, sizeof(image3d_max_height), &image3d_max_height, NULL);
    size_t image3d_max_depth;
    clGetDeviceInfo(device, CL_DEVICE_IMAGE3D_MAX_DEPTH, sizeof(image3d_max_depth), &image3d_max_depth, NULL);
    fprintf(stdout, "  CL_DEVICE_IMAGE_MAX_WIDTH:\t\t2d width %d, 2d height %d, 3d width %d, 3d height %d, 3d depth %d\n", (int)image2d_max_width, (int)image2d_max_height, (int)image3d_max_width, (int)image3d_max_height, (int)image3d_max_depth);

    // CL_DEVICE_MAX_MEM_ALLOC_SIZE
    cl_ulong max_mem_alloc_size;
    clGetDeviceInfo(device, CL_DEVICE_MAX_MEM_ALLOC_SIZE, sizeof(max_mem_alloc_size), &max_mem_alloc_size, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_MEM_ALLOC_SIZE:\t\t%d MByte\n", (unsigned int)(max_mem_alloc_size / (1024 * 1024)));

    // CL_DEVICE_GLOBAL_MEM_SIZE
    cl_ulong mem_size;
    clGetDeviceInfo(device, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
    fprintf(stdout, "  CL_DEVICE_GLOBAL_MEM_SIZE:\t\t%d MByte\n", (unsigned int)(mem_size / (1024 * 1024)));

    // CL_DEVICE_ERROR_CORRECTION_SUPPORT
    cl_bool error_correction_support;
    clGetDeviceInfo(device, CL_DEVICE_ERROR_CORRECTION_SUPPORT, sizeof(error_correction_support), &error_correction_support, NULL);
    fprintf(stdout, "  CL_DEVICE_ERROR_CORRECTION_SUPPORT:\t%s\n", error_correction_support == CL_TRUE ? "yes" : "no");

    // CL_DEVICE_LOCAL_MEM_TYPE
    cl_device_local_mem_type local_mem_type;
    clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_TYPE, sizeof(local_mem_type), &local_mem_type, NULL);
    fprintf(stdout, "  CL_DEVICE_LOCAL_MEM_TYPE:\t\t%s\n", local_mem_type == 1 ? "local" : "global");

    // CL_DEVICE_LOCAL_MEM_SIZE
    clGetDeviceInfo(device, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
    fprintf(stdout, "  CL_DEVICE_LOCAL_MEM_SIZE:\t\t%d KByte\n", (unsigned int)(mem_size / 1024));

    // CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE
    clGetDeviceInfo(device, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(mem_size), &mem_size, NULL);
    fprintf(stdout, "  CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE:\t%d KByte\n", (unsigned int)(mem_size / 1024));

    // CL_DEVICE_QUEUE_PROPERTIES
    cl_command_queue_properties queue_properties;
    clGetDeviceInfo(device, CL_DEVICE_QUEUE_PROPERTIES, sizeof(queue_properties), &queue_properties, NULL);
    if( queue_properties & CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE )
        fprintf(stdout, "  CL_DEVICE_QUEUE_PROPERTIES:\t\t%s\n", "CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE");
    if( queue_properties & CL_QUEUE_PROFILING_ENABLE )
        fprintf(stdout, "  CL_DEVICE_QUEUE_PROPERTIES:\t\t%s\n", "CL_QUEUE_PROFILING_ENABLE");

    // CL_DEVICE_EXTENSIONS: get device extensions, and if any then parse & log the string onto separate lines
    clGetDeviceInfo(device, CL_DEVICE_EXTENSIONS, sizeof(device_string), &device_string, NULL);
    if (device_string != 0)
    {
        fprintf(stdout, "  CL_DEVICE_EXTENSIONS:\n");
        std::string stdDevString;
        stdDevString = std::string(device_string);
        size_t szOldPos = 0;
        size_t szSpacePos = stdDevString.find(' ', szOldPos); // extensions string is space delimited
        while (szSpacePos != stdDevString.npos && (szSpacePos - szOldPos) > 0)
        {
            if( strcmp("cl_nv_device_attribute_query", stdDevString.substr(szOldPos, szSpacePos - szOldPos).c_str()) == 0 )
                nv_device_attibute_query = true;

            fprintf(stdout, "\t\t\t\t\t%s\n", stdDevString.substr(szOldPos, szSpacePos - szOldPos).c_str());

            szOldPos = szSpacePos + 1;
            szSpacePos = stdDevString.find(' ', szOldPos);
        }
    }
    else
    {
        fprintf(stdout, "  CL_DEVICE_EXTENSIONS: None\n");
    }

    // CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR
    cl_uint vec_width_char;
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_CHAR, sizeof(vec_width_char), &vec_width_char, NULL);
    cl_uint vec_width_short;
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_SHORT, sizeof(vec_width_short), &vec_width_short, NULL);
    cl_uint vec_width_int;
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_INT, sizeof(vec_width_int), &vec_width_int, NULL);
    cl_uint vec_width_long;
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_LONG, sizeof(vec_width_long), &vec_width_long, NULL);
    cl_uint vec_width_float;
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_FLOAT, sizeof(vec_width_float), &vec_width_float, NULL);
    cl_uint vec_width_double;
    clGetDeviceInfo(device, CL_DEVICE_PREFERRED_VECTOR_WIDTH_DOUBLE, sizeof(vec_width_double), &vec_width_double, NULL);
    fprintf(stdout, "  CL_DEVICE_PREFERRED_VECTOR_WIDTH:\tchar %d, short %d, int %d, long %d, float %d, double %d\n", vec_width_char, vec_width_short, vec_width_int, vec_width_long, vec_width_float, vec_width_double);

    fprintf(stdout, "\n");
}

///
/// Loads a Program file and prepends the cPreamble to the code. (from the NVIDIA SDK)
///
/// @return the source string if succeeded, 0 otherwise
/// @param cFilename        program filename
/// @param cPreamble        code that is prepended to the loaded file, typically a set of #defines or a header
/// @param szFinalLength    returned length of the code string
///
char* oclLoadProgSource(const char* cFilename, const char* cPreamble, size_t* szFinalLength)
{
    // locals
    FILE* pFileStream = NULL;
    size_t szSourceLength;

    // open the OpenCL source code file
    pFileStream = fopen(cFilename, "rb");
    if(pFileStream == 0)
    {
        return NULL;
    }

    size_t szPreambleLength = strlen(cPreamble);

    // get the length of the source code
    fseek(pFileStream, 0, SEEK_END);
    szSourceLength = ftell(pFileStream);
    fseek(pFileStream, 0, SEEK_SET);

    // allocate a buffer for the source code string and read it in
    char* cSourceString = (char *)malloc(szSourceLength + szPreambleLength + 1);
    memcpy(cSourceString, cPreamble, szPreambleLength);
    if (fread((cSourceString) + szPreambleLength, szSourceLength, 1, pFileStream) != 1)
    {
        fclose(pFileStream);
        free(cSourceString);
        return 0;
    }

    // close the file and return the total length of the combined (preamble + source) string
    fclose(pFileStream);
    if(szFinalLength != 0)
    {
        *szFinalLength = szSourceLength + szPreambleLength;
    }
    cSourceString[szSourceLength + szPreambleLength] = '\0';

    return cSourceString;
}

///
/// Gets the id of the nth device from the context (from the NVIDIA SDK)
///
/// @return the id or -1 when out of range
/// @param cxGPUContext         OpenCL context
/// @param device_idx            index of the device of interest
///
cl_device_id oclGetDev(cl_context cxGPUContext, unsigned int nr)
{
    size_t szParmDataBytes;
    cl_device_id* cdDevices;

    // get the list of GPU devices associated with context
    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, 0, NULL, &szParmDataBytes);

    if( szParmDataBytes / sizeof(cl_device_id) < nr )
    {
        return (cl_device_id)-1;
    }

    cdDevices = (cl_device_id*) malloc(szParmDataBytes);

    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, szParmDataBytes, cdDevices, NULL);

    cl_device_id device = cdDevices[nr];
    free(cdDevices);

    return device;
}

///
/// Gets the id of the first device from the context (from the NVIDIA SDK)
///
/// @return the id
/// @param cxGPUContext         OpenCL context
///
cl_device_id oclGetFirstDev(cl_context cxGPUContext)
{
    size_t szParmDataBytes;
    cl_device_id* cdDevices;

    // get the list of GPU devices associated with context
    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, 0, NULL, &szParmDataBytes);
    cdDevices = (cl_device_id*) malloc(szParmDataBytes);

    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, szParmDataBytes, cdDevices, NULL);

    cl_device_id first = cdDevices[0];
    free(cdDevices);

    return first;
}

///
/// Gets the id of device with maximal FLOPS from the context (from NVIDIA SDK)
///
/// @return the id
/// @param cxGPUContext         OpenCL context
///
cl_device_id oclGetMaxFlopsDev(cl_context cxGPUContext)
{
    size_t szParmDataBytes;
    cl_device_id* cdDevices;

    // get the list of GPU devices associated with context
    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, 0, NULL, &szParmDataBytes);
    cdDevices = (cl_device_id*) malloc(szParmDataBytes);
    size_t device_count = szParmDataBytes / sizeof(cl_device_id);

    clGetContextInfo(cxGPUContext, CL_CONTEXT_DEVICES, szParmDataBytes, cdDevices, NULL);

    cl_device_id max_flops_device = cdDevices[0];
	int max_flops = 0;

	size_t current_device = 0;

    // CL_DEVICE_MAX_COMPUTE_UNITS
    cl_uint compute_units;
    clGetDeviceInfo(cdDevices[current_device], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(compute_units), &compute_units, NULL);

    // CL_DEVICE_MAX_CLOCK_FREQUENCY
    cl_uint clock_frequency;
    clGetDeviceInfo(cdDevices[current_device], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clock_frequency), &clock_frequency, NULL);

	max_flops = compute_units * clock_frequency;
	++current_device;

	while( current_device < device_count )
	{
        // CL_DEVICE_MAX_COMPUTE_UNITS
        cl_uint compute_units;
        clGetDeviceInfo(cdDevices[current_device], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(compute_units), &compute_units, NULL);

        // CL_DEVICE_MAX_CLOCK_FREQUENCY
        cl_uint clock_frequency;
        clGetDeviceInfo(cdDevices[current_device], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(clock_frequency), &clock_frequency, NULL);

        int flops = compute_units * clock_frequency;
		if( flops > max_flops )
		{
			max_flops        = flops;
			max_flops_device = cdDevices[current_device];
		}
		++current_device;
	}

    free(cdDevices);

	return max_flops_device;
}

///
/// Check for error condition and exit if found.  Print file and line number
/// of error.
///
/// @param errNum Error value to check
/// @param expected Expected value
/// @param file File name (filled in by macro)
/// @param lineNumber Line number (filled in by macro)
///
void oclCheckErrorFileLine(int errNum, int expected, const char* file, const int lineNumber)
{
    if (errNum != expected)
    {
        cerr << "Line " << lineNumber << " in File " << file << endl;
        exit(1);
    }
}


///
/// Round the local work size up to the next multiple of the size
/// @param groupSize Size of the group
/// @param globalSize Global size
/// @return The rounded local work size
///
int oclRoundWorkSizeUp(int groupSize, int globalSize)
{
    int remainder = globalSize % groupSize;
    if (remainder == 0)
    {
        return globalSize;
    }
    else
    {
        return globalSize + groupSize - remainder;
    }
}
