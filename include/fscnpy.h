#ifndef FSCNPY_H
#define FSCNPY_H

#include <vector>

#include "mri.h"
#include "bfileio.h"

#define __FSCNPY_DEBUG 0

class fscnpy
{
public:
  fscnpy();
  ~fscnpy();

  unsigned char *read(const char *npy);
  MRI *npy2mri(const char *npy, bool verbose=false);
  void write(const char *npy, unsigned char *data);

private:
  void __parse_header_dictionary(const char *npy);
  int  __getmritype();
  char __getarchendian();
  void __contiguousreorder(MRI *mri);

  const char *__magic_string = "\x93NUMPY";
  static const int  __magic_string_len = 6;
  static const char __little_endian_char = '<';
  static const char __big_endian_char = '>';
  static const char __no_endian_char = '|';
  const int  __header_len_bytes[2] = {2, 4};

  unsigned char __major_version;
  unsigned char __minor_version;
  //unsigned char __header_len_str;
  size_t __header_len;

  char __endian;
  unsigned char __dtype_char;
  int __dtype_len;
  bool __fortran_order;
  std::vector<int> __shape;

  unsigned char *__data;

  int __dtype; // this is converted to fs mri types
  int __ndims;
  bool __squeezed;

  FILE *__npyfp;
  char __archendian;

  bool __verbose;
  int  __counter;
};

#endif
