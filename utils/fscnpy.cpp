#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "fscnpy.h"

fscnpy::fscnpy()
{
  __major_version = 0;
  __minor_version  = 0;
  __header_len = 0;

  /* get endianness of current architecture */
  __archendian = __getarchendian();

  __npyfp = NULL;

  __dtype = -1;
  __data = NULL;

  __squeezed = false;

  __verbose = false;
  __counter = 0;
}


fscnpy::~fscnpy()
{
  if (__data != NULL)
  {
    free(__data);
    __data = NULL;
  }
}


unsigned char *fscnpy::read(const char *npy)
{
  __npyfp = fopen(npy, "rb");
  if (!__npyfp)
  {
    printf("Error reading numpy file %s\n", npy);
    exit(1);
  }

  __parse_header_dictionary(npy);

  // read the data
  size_t ndata = __shape[0] * __shape[1] * __shape[2];
  size_t datasize = __dtype_len * ndata;
  __data = new unsigned char[datasize];
  memset(__data, 0, datasize);
  size_t bytesread = fread(__data, __dtype_len, ndata, __npyfp);
  if (bytesread != datasize)
  {
    printf("ERROR failed to read %ld bytes data (%ld read)\n", datasize, bytesread);
  }

  fclose(__npyfp);

  return __data;
}


MRI* fscnpy::npy2mri(const char *npy, bool verbose)
{
#if __FSCNPY_DEBUG
  ++__counter;
  printf("fscnpy::npy2mri() called: %d\n", __counter);
#endif

  __npyfp = fopen(npy, "rb");
  if (!__npyfp)
  {
    printf("Error reading numpy file %s\n", npy);
    exit(1);
  }

  __verbose = verbose;
  __parse_header_dictionary(npy);

  if (__ndims > 4)
  {
    printf("numpy file contains %dD data. Cannot fit them in MRI structure.\n", __ndims);
    exit(1);
  }

  // create MRI structure
  __getmritype();
  MRI *mri = new MRI(__shape, __dtype);

  if (__verbose)
  {
    printf("MRI dimensions: %d x %d x %d x %d\n", mri->width, mri->height, mri->depth, mri->nframes);
    printf("          type: %s (%d)\n",
           mri->type == MRI_UCHAR   ? "UCHAR" :
           mri->type == MRI_SHORT   ? "SHORT" :
           mri->type == MRI_USHRT   ? "USHRT" :
           mri->type == MRI_INT     ? "INT" :
           mri->type == MRI_LONG    ? "LONG" :
           mri->type == MRI_BITMAP  ? "BITMAP" :
           mri->type == MRI_TENSOR  ? "TENSOR" :
           mri->type == MRI_FLOAT   ? "FLOAT" : "UNKNOWN", mri->type);
    fflush(stdout);
  }

  size_t bytesread = 0;
#if 1
    size_t bytes_per_slice = __dtype_len * mri->width * mri->height;
    for (int f = 0; f < mri->nframes; f++)
    {
      for (int s = 0; s < mri->depth; s++)
      {
        void *buf = &MRIseq_vox(mri, 0, 0, s, f);

        bytesread += fread(buf, 1, bytes_per_slice, __npyfp);
        if (__archendian != __endian)
        {
          // swap bytes
          if (__dtype_len == 2) 
            byteswapbufshort(__data, bytes_per_slice);
          if (__dtype_len == 4) 
            byteswapbuffloat(__data, bytes_per_slice);
          if (__dtype_len == 8) 
            byteswapbuffloat(__data, bytes_per_slice);
	}

        // freeview progress bar
        exec_progress_callback(s, mri->depth, f, mri->nframes);
      }
    }
#else
  // load data into MRI structure
  size_t bytes_per_col = __dtype_len * mri->width;
  for (int f = 0; f < mri->nframes; f++)
  {
    for (int s = 0; s < mri->depth; s++)
    {
      for (int r = 0; r < mri->height; r++)
      {
        void *buf = &MRIseq_vox(mri, 0, r, s, f);

        bytesread += fread(buf, 1, bytes_per_col, __npyfp);
        if (__archendian != __endian)
        {
          // swap bytes
          if (__dtype_len == 2) 
            byteswapbufshort(buf, bytes_per_col);
          if (__dtype_len == 4) 
            byteswapbuffloat(buf, bytes_per_col);
          if (__dtype_len == 8) 
            byteswapbuffloat(buf, bytes_per_col);
        }
      }
    }
  }
#endif

  fclose(__npyfp);

  if (__verbose)
  {
    printf("%ld bytes read\n", bytesread);
    fflush(stdout);
  }

  // need to swap the bytes if numpy array is in C-order
  if (!__fortran_order)
    __contiguousreorder(mri);

  return mri;
}


void fscnpy::write(const char *npy, unsigned char *data)
{
}



/*
 * https://docs.scipy.org/doc/numpy-1.14.0/neps/npy-format.html
 *
 * npy format specification V1.0:
 * The first 6 bytes are a magic string: exactly “x93NUMPY”.
 * The next 1 byte is an unsigned byte: the major version number of the file format, e.g. x01.
 * The next 1 byte is an unsigned byte: the minor version number of the file format, e.g. x00. Note: the version of the file format is not tied to the version of the numpy package.
 * The next 2 bytes form a little-endian unsigned short int: the length of the header data HEADER_LEN.
 * The next HEADER_LEN bytes form the header data describing the array’s format. It is an ASCII string which contains a Python literal expression of a dictionary. 
 *     It is terminated by a newline (‘\n’) and padded with spaces (‘x20’) to make the total length of the magic string + 4 + HEADER_LEN be evenly divisible by 16 for alignment purposes.
 * The dictionary contains three keys:
 *     “descr” : dtype.descr       An object that can be passed as an argument to the numpy.dtype() constructor to create the array’s dtype.
 *     “fortran_order” : bool      Whether the array data is Fortran-contiguous or not. Since Fortran-contiguous arrays are a common form of non-C-contiguity, 
 *                                 we allow them to be written directly to disk for efficiency.
 *     “shape” : tuple of int      The shape of the array.
 *
 *
 * npy format specification V2.0:
 * The description of the fourth element of the header therefore has become:
 *     The next 4 bytes form a little-endian unsigned int: the length of the header data HEADER_LEN.
 *
 *
 * npy format specification V3.0:    (https://numpy.org/devdocs/reference/generated/numpy.lib.format.html)
 *       This version replaces the ASCII string (which in practice was latin1) with a utf8-encoded string, so supports structured types with any unicode field names.
 *
 *
 * sample dictionary:
 *  {   'descr' : '<i8', 
 *      'fortran_order' : False,
 *      'shape' : (10   ,),       
 *  }
 */
void fscnpy::__parse_header_dictionary(const char *npy)
{
  // read magic string : 6 bytes
  char magic_string[6]={'\0'};
  fread(magic_string, 1, __magic_string_len, __npyfp);
  if (strcmp(magic_string, __magic_string) != 0)
  {
    printf("Error %s is not npy file\n", npy);
    exit(1);
  }

  // read major  version : 1 bytes
  // read mminor version : 1 bytes
  fread(&__major_version, 1, 1, __npyfp);
  fread(&__minor_version, 1, 1, __npyfp);

  // read header_len : 2 bytes for version 1, and 4 bytes for version 2
  // 2 bytes form a little-endian unsigned short int: the length of the header data HEADER_LEN
  // ???on big endian machines, the bytes need to be swapped???
  fread(&__header_len, 1, __header_len_bytes[__major_version-1], __npyfp);

  // read dictionary
  char dict[__header_len];
  fread(dict, 1, __header_len, __npyfp);

  // skip the rest of padded with spaces
  if ( ((__magic_string_len + 2 + __header_len_bytes[__major_version-1] + __header_len) % 16) != 0 )
  {
    size_t restheader_len = 16 - (__magic_string_len + 2 + __header_len_bytes[__major_version-1] + __header_len) % 16;
    unsigned char restheader[restheader_len];
    fread(restheader, 1, restheader_len, __npyfp);
  }

  char *currpos = dict;
  for (; *currpos != '\n'; currpos++)
  {
    if (*currpos == '{' || *currpos == '}' || *currpos == ',' || *currpos == ' ')
      continue;

    char *key   = currpos;
    while (*currpos != ':' && *currpos != '\n')
      currpos++;
    
    if (*currpos == ':')
    {
      *currpos = '\0';
      currpos++;
      while (*currpos == ' ')
        currpos++;

      if (strcmp(key, "\'descr\'") == 0)
      {
        char *value = ++currpos;  // skip '
        __endian = value[0];
        __dtype_char  = value[1];
        currpos += 2;

        value = currpos;
        while (isdigit(*currpos))
	  currpos++;

        *currpos = '\0';
        __dtype_len   = atoi(value);
      } // descr
      else if (strcmp(key, "\'fortran_order\'") == 0)
      {
        char *value = currpos;
        while (*currpos != ',')
          currpos++;

        if (*currpos == ',')
          *currpos = '\0';
       
        __fortran_order = !strcmp(value, "False") ? false : true;
      } // fortran_order
      else if (strcmp(key, "\'shape\'") == 0)
      {
	while (!isdigit(*currpos))
          currpos++;

        char *value = currpos;
        while (true)
	{
          if (isdigit(*currpos) && value == NULL)
            value = currpos;
          else if (*currpos == ',' || *currpos == ')')
	  {
            bool done = false;
            if (*currpos == ')')
              done = true;

            *currpos = '\0';
            __shape.push_back(atoi(value));
            value = NULL;
            if (done)
	    {
              currpos++;
              break;
	    }
	  }

          currpos++;
        }
      } // shape
    } 
  } // for

  if (__verbose)
  {
    printf("archendian : %c (%s)\n", __archendian, (__archendian == '<') ? "little endian" : "big endian");
    printf("NPY version: %d.%d\n", __major_version, __minor_version);
    printf("NPY header_len: %ld\n", __header_len);
    printf("NPY dict.descr : %c%c%d\n", __endian, __dtype_char, __dtype_len);
    printf("NPY dict.fortran_order: %s\n", (__fortran_order) ? "True" : "False");
    printf("NPY dict.shape: {");
    for (int i = 0; i < __shape.size(); i++)
      printf(" %d, ", __shape[i]);
    printf("}\n");
  }


  __ndims = __shape.size();
  // remove the last dimensions of length 1
  while (__shape[__ndims-1] == 1)
  {
    __squeezed = true;
    __shape.pop_back();
    __ndims--;
  }

  // remove the first dimensions of length 1
  while (__shape[0] == 1)
  {
    __shape.erase(__shape.begin());
    __ndims--;
    __squeezed = true;
  }

  if (__verbose)
  {
    printf("NPY shape dimension: %d (%s) ( ", __ndims, (__squeezed) ? "squeezed" : "");
    for (int n = 0; n < __shape.size(); n++)
      printf("%d, ", __shape[n]);
    printf(")\n");

    if (__archendian != __endian)
      printf("Data is in different endian. Need to swap bytes\n");
    else
      printf("Data is in same endian\n");
  
    fflush(stdout);
  }
}


// 
int fscnpy::__getmritype()
{
  switch (__dtype_char)
  {
    case 'f':
    {
      if (__dtype_len == sizeof(float))
        __dtype = MRI_FLOAT;
      else
        printf("WARN: !!! unsupported __dtype_char=%c, __dtype_len=%d\n", __dtype_char, __dtype_len); 
      break;
    }
    case 'i':
    {
      if (__dtype_len == sizeof(short))
        __dtype = MRI_SHORT;
      else if (__dtype_len == sizeof(int))
        __dtype = MRI_INT;
      else if (__dtype_len == sizeof(long))
        __dtype = MRI_LONG;
      else
        printf("WARN: !!! unsupported __dtype_char=%c, __dtype_len=%d\n", __dtype_char, __dtype_len); 
      break;
    }
    case 'u':
    {
      if (__dtype_len == sizeof(unsigned char))
        __dtype = MRI_UCHAR;
      else if (__dtype_len == sizeof(unsigned short))
        __dtype = MRI_USHRT;
      else
        printf("WARN: !!! unsupported __dtype_char=%c, __dtype_len=%d\n", __dtype_char, __dtype_len); 
      break;
    }
    default:
      printf("WARN: !!! unsupported __dtype_char=%c!!!\n", __dtype_char);
      break;
  }

  return __dtype;
}


// return '<' for little endian;
// return '>' for big endian 
char fscnpy::__getarchendian() 
{
    int x = 1;
    return (((char *)&x)[0]) ? '<' : '>';
}


void fscnpy::__contiguousreorder(MRI *mri)
{
  if (__shape.size() > 3)
  {
    printf("__indexreorder() not implemented for 4d mri\n");
    return;
  }

  if (__dtype != MRI_UCHAR && __dtype != MRI_SHORT && __dtype != MRI_USHRT && __dtype != MRI_INT && __dtype != MRI_LONG && __dtype != MRI_FLOAT)
  {
    printf("Unsupported datatype, no index reordering\n");
    return;
  }

  if (__verbose)
  {
    printf("array contiguous reordering ...\n");
    fflush(stdout);
  }

  for (int s = 0; s < mri->depth; s++)
  {
    for (int r = 0; r < mri->height; r++)
    {
      for (int c = 0; c < mri->width; c++)
      {
        if (s <= c)
          continue;

#if __FSCNPY_DEBUG
        printf("swap [%d, %d, %d] with [%d, %d, %d]\n", c, r, s, s, r, c);

        printf("(before) loc1 = %.7f\n", MRIFvox(mri, c, r, s));
        printf("(before) loc2 = %.7f\n", MRIFvox(mri, s, r, c));
        fflush(stdout);
#endif

        if (__dtype == MRI_UCHAR)
	{
          unsigned char tmp = MRIvox(mri, c, r, s);
          MRIvox(mri, c, r, s) = MRIvox(mri, s, r, c);
          MRIvox(mri, s, r, c) = tmp;
        }
        else if (__dtype == MRI_SHORT)
	{
          short tmp = MRISvox(mri, c, r, s);
          MRISvox(mri, c, r, s) = MRISvox(mri, s, r, c);
          MRISvox(mri, s, r, c) = tmp;
        }
        else if (__dtype == MRI_USHRT)
	{
          unsigned short tmp = MRIUSvox(mri, c, r, s);
          MRIUSvox(mri, c, r, s) = MRIUSvox(mri, s, r, c);
          MRIUSvox(mri, s, r, c) = tmp;
        }
        else if (__dtype == MRI_INT)
	{
          int tmp = MRIIvox(mri, c, r, s);
          MRIIvox(mri, c, r, s) = MRIIvox(mri, s, r, c);
          MRIIvox(mri, s, r, c) = tmp;
	}
        else if (__dtype == MRI_LONG)
	{
          long tmp = MRILvox(mri, c, r, s);
          MRILvox(mri, c, r, s) = MRILvox(mri, s, r, c);
          MRILvox(mri, s, r, c) = tmp;
	}
        else if (__dtype == MRI_FLOAT)
	{
          float tmp = MRIFvox(mri, c, r, s);
          MRIFvox(mri, c, r, s) = MRIFvox(mri, s, r, c);
          MRIFvox(mri, s, r, c) = tmp;
	}

#if __FSCNPY_DEBUG
        printf("(after) loc1 = %.7f\n", MRIFvox(mri, c, r, s));
        printf("(after) loc2 = %.7f\n", MRIFvox(mri, s, r, c));
        fflush(stdout);
        printf("\n");
#endif
      }
    }
  }
}




