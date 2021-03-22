/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

/*-----------------------------------------------------
                FFT Filter

We have implemented a FFT-filtering. Which means we apply
first a FFT to the mri, then we apply a filter; and we
finally apply en inverse FFT.
Only one filter has been implemented so far: The Gaussian

We have first some basic functions used in the FFT Filter.

-----------------------------------------------------*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define FourierForward 1
#define FourierBackward -1

const int cMaxLength = 4096;
const int cMinLength = 1;
static int _rev = 0;
int *reversedBits;
int **_reverseBits = NULL;
static int _lookupTabletLength = -1;
typedef struct complex
{
  double a, b;
} complex;
typedef struct complexF
{
  float a, b;
} complexF;

static complex **_uRLookup = NULL;
static complex **_uILookup = NULL;
static complexF **_uRLookupF = NULL;
static complexF **_uILookupF = NULL;

static float *uRLookup = NULL;
static float *uILookup = NULL;

static void FFTerror(const char *string)
{
  fprintf(stdout, "\nFFT Error: %s\n", string);
  exit(1);
}

void FFTdebugAssert(int b, const char *string)
{
  if (!b) FFTerror(string);
}

static void Swap(float *a, float *b)
{
  float temp = *a;
  *a = *b;
  *b = temp;
}
int FFTisPowerOf2(int x) { return (x & (x - 1)) == 0; }
int FFTpow2(int exponent)
{
  if (exponent >= 0 && exponent < 31) return 1 << exponent;
  return 0;
}
int FFTlog2(int x)
{  // ceiling of FFTlog2
  if (x <= 65536) {
    if (x <= 256) {
      if (x <= 16) {
        if (x <= 4) {
          if (x <= 2) {
            if (x <= 1) {
              return 0;
            }
            return 1;
          }
          return 2;
        }
        if (x <= 8) return 3;
        return 4;
      }
      if (x <= 64) {
        if (x <= 32) return 5;
        return 6;
      }
      if (x <= 128) return 7;
      return 8;
    }
    if (x <= 4096) {
      if (x <= 1024) {
        if (x <= 512) return 9;
        return 10;
      }
      if (x <= 2048) return 11;
      return 12;
    }
    if (x <= 16384) {
      if (x <= 8192) return 13;
      return 14;
    }
    if (x <= 32768) return 15;
    return 16;
  }
  if (x <= 16777216) {
    if (x <= 1048576) {
      if (x <= 262144) {
        if (x <= 131072) return 17;
        return 18;
      }
      if (x <= 524288) return 19;
      return 20;
    }
    if (x <= 4194304) {
      if (x <= 2097152) return 21;
      return 22;
    }
    if (x <= 8388608) return 23;
    return 24;
  }
  if (x <= 268435456) {
    if (x <= 67108864) {
      if (x <= 33554432) return 25;
      return 26;
    }
    if (x <= 134217728) return 27;
    return 28;
  }
  if (x <= 1073741824) {
    if (x <= 536870912) return 29;
    return 30;
  }
  return 31;
}

static void ReorderArray(float *data, int data_length)
{
  FFTdebugAssert(data != NULL, "ReorderArray : Data = Null");
  int length = data_length / 2;
  int numberOfBits = FFTlog2(length);
  int maxBits = FFTpow2(numberOfBits);
  int i, j;
  FFTdebugAssert(FFTisPowerOf2(length) == 1, "length no power of 2");
  FFTdebugAssert(length >= cMinLength, "length < cMinLength");
  FFTdebugAssert(length <= cMaxLength, "length > cMaxLength");
  if (_rev == 0) {
    reversedBits = (int *)malloc(512 * sizeof(int));
    for (i = 0; i < maxBits; i++) {
      int oldBits = i;
      int newBits = 0;
      for (j = 0; j < numberOfBits; j++) {
        newBits = (newBits << 1) | (oldBits & 1);
        oldBits = (oldBits >> 1);
      }
      reversedBits[i] = newBits;
    }
    _rev = 1;
  }
  for (i = 0; i < length; i++) {
    int swap = reversedBits[i];
    if (swap > i) {
      Swap(&data[(i << 1)], &data[(swap << 1)]);
      Swap(&data[(i << 1) + 1], &data[(swap << 1) + 1]);
    }
  }
}

static int _ReverseBits(int bits, int n)
{
  int bitsReversed = 0;
  int i;
  for (i = 0; i < n; i++) {
    bitsReversed = (bitsReversed << 1) | (bits & 1);
    bits = (bits >> 1);
  }
  return bitsReversed;
}
static void InitializeReverseBits(int levels)
{
  int i, j;
  if (_reverseBits) free(_reverseBits);
  _reverseBits = (int **)malloc((levels + 1) * sizeof(int *));
  for (j = 0; j < (levels + 1); j++) {
    int count = (int)FFTpow2(j);
    _reverseBits[j] = (int *)malloc(count * sizeof(int));
    for (i = 0; i < count; i++) {
      _reverseBits[j][i] = _ReverseBits(i, j);
    }
  }
}

static void InitializeComplexRotations(int levels)
{
  int ln = levels;

  if (_uRLookup) free(_uRLookup);
  if (_uILookup) free(_uILookup);
  if (_uRLookupF) free(_uRLookupF);
  if (_uILookupF) free(_uILookupF);
  _uRLookup = (complex **)malloc((levels + 1) * sizeof(complex *));
  _uILookup = (complex **)malloc((levels + 1) * sizeof(complex *));
  _uRLookupF = (complexF **)malloc((levels + 1) * sizeof(complexF *));
  _uILookupF = (complexF **)malloc((levels + 1) * sizeof(complexF *));
  int level, j, N = 1;
  for (level = 1; level <= ln; level++) {
    int M = N;
    N <<= 1;
    {
      // positive sign ( i.e. [M,0] )
      double uR = 1;
      double uI = 0;
      double angle = (double)M_PI / M * 1;
      double wR = (double)cos(angle);
      double wI = (double)sin(angle);

      _uRLookup[level] = (complex *)malloc(M * sizeof(complex));
      _uILookup[level] = (complex *)malloc(M * sizeof(complex));
      _uRLookupF[level] = (complexF *)malloc(M * sizeof(complexF));
      _uILookupF[level] = (complexF *)malloc(M * sizeof(complexF));

      for (j = 0; j < M; j++) {
        _uRLookupF[level][j].a = (float)(_uRLookup[level][j].a = uR);
        _uILookupF[level][j].a = (float)(_uILookup[level][j].a = uI);
        double uwI = uR * wI + uI * wR;
        uR = uR * wR - uI * wI;
        uI = uwI;
      }
    }
    {
      // negative sign ( i.e. [M,1] )
      double uR = 1;
      double uI = 0;
      double angle = (double)M_PI / M * -1;
      double wR = (double)cos(angle);
      double wI = (double)sin(angle);

      for (j = 0; j < M; j++) {
        _uRLookupF[level][j].b = (float)(_uRLookup[level][j].b = uR);
        _uILookupF[level][j].b = (float)(_uILookup[level][j].b = uI);
        double uwI = uR * wI + uI * wR;
        uR = uR * wR - uI * wI;
        uI = uwI;
      }
    }
  }
}

static void SyncLookupTableLength(int length)
{
  FFTdebugAssert(length < 1024 * 10, "SyncLookupTableLength : length too big");
  FFTdebugAssert(length >= 0, "SyncLookupTableLength : length<0");
  if (length > _lookupTabletLength) {
    int level = FFTlog2(length);
    InitializeReverseBits(level);
    InitializeComplexRotations(level);
    _lookupTabletLength = length;
  }
}

static void copy_vect(float *vect, complexF **mat, int level, int signIndex, int M)
{
  int j;
  if (signIndex)
    for (j = 0; j < M; j++) vect[j] = mat[level][j].b;
  else
    for (j = 0; j < M; j++) {
      vect[j] = mat[level][j].a;
    }
}
/*-----------------------------------------------------
 FFT performs a complex FFT where the complex numbers are
 represented as pairs in the vector data. That is wfy we
 should have data_length >= length*2.
 The direction means either from time to frequency, or
 from frequency to time.
------------------------------------------------------*/
static void FFT(float *data, int data_length, int length, int direction)
{
  FFTdebugAssert(data != NULL, "FFT : DATA = NULL");
  FFTdebugAssert(data_length >= length * 2, "FFT : data_length < length*2");
  FFTdebugAssert(FFTisPowerOf2(length) == 1, "FFT : length is not a power of 2");

  SyncLookupTableLength(length);
  int ln = FFTlog2(length);

  // reorder array
  ReorderArray(data, data_length);

  // successive doubling
  int N = 1;
  int signIndex = (direction == FourierForward) ? 0 : 1;

  int level, j, evenT;
  for (level = 1; level <= ln; level++) {
    int M = N;
    N <<= 1;

    uRLookup = (float *)malloc(M * sizeof(float));
    uILookup = (float *)malloc(M * sizeof(float));
    copy_vect(uRLookup, _uRLookupF, level, signIndex, M);
    copy_vect(uILookup, _uILookupF, level, signIndex, M);

    for (j = 0; j < M; j++) {
      float uR = uRLookup[j];
      float uI = uILookup[j];

      for (evenT = j; evenT < length; evenT += N) {
        int even = evenT << 1;
        int odd = (evenT + M) << 1;
        float r = data[odd];
        float i = data[odd + 1];

        float odduR = r * uR - i * uI;
        float odduI = r * uI + i * uR;

        r = data[even];
        i = data[even + 1];

        data[even] = r + odduR;
        data[even + 1] = i + odduI;

        data[odd] = r - odduR;
        data[odd + 1] = i - odduI;
      }
    }
  }
  free(uRLookup);
  free(uILookup);
}

/*!
 \fn void RFFT( float* data, int data_length, int length, int direction )
 \param data - input (and output) vector
 \param data_length - length of input
 \brief RFFT performs a real FFT: data is a vector with the real
 input with n nb. The complex solution (n complexes) is
 represented in the same vector data, this way :
 if the solution is :
 A	H	G	F	E	F	G	H
        -h*i	-c*i	-f*i		+f*i	+g*i	+h*i
 The vector data will be :
 A	E	H	h	G	g	F	f
*/
void RFFT(float *data, int data_length, int length, int direction)
{
  FFTdebugAssert(data != NULL, "RFFT :data");
  FFTdebugAssert(data_length >= length, "RFFT : length must be at least as large as data_length parameter");
  FFTdebugAssert(FFTisPowerOf2(length) == 1, "RFFT : length must be a power of 2");
  float c1 = 0.5f, c2;
  float theta = (float)M_PI / (length / 2);
  int i;

  if (direction == FourierForward) {
    c2 = -0.5f;
    FFT(data, data_length, length / 2, direction);
  }
  else {
    c2 = 0.5f;
    theta = -theta;
  }

  float wtemp = (float)sin(0.5 * theta);
  float wpr = -2 * wtemp * wtemp;
  float wpi = (float)sin(theta);
  float wr = 1 + wpr;
  float wi = wpi;

  // do / undo packing
  for (i = 1; i < length / 4; i++) {
    int a = 2 * i;
    int b = length - 2 * i;
    float h1r = c1 * (data[a] + data[b]);
    float h1i = c1 * (data[a + 1] - data[b + 1]);
    float h2r = -c2 * (data[a + 1] + data[b + 1]);
    float h2i = c2 * (data[a] - data[b]);
    data[a] = h1r + wr * h2r - wi * h2i;
    data[a + 1] = h1i + wr * h2i + wi * h2r;
    data[b] = h1r - wr * h2r + wi * h2i;
    data[b + 1] = -h1i + wr * h2i + wi * h2r;
    wr = (wtemp = wr) * wpr - wi * wpi + wr;
    wi = wi * wpr + wtemp * wpi + wi;
  }

  if (direction == FourierForward) {
    float hir = data[0];
    data[0] = hir + data[1];
    data[1] = hir - data[1];
  }
  else {
    float hir = data[0];
    data[0] = c1 * (hir + data[1]);
    data[1] = c1 * (hir - data[1]);
    FFT(data, data_length, length / 2, direction);
  }
}
/*-----------------------------------------------------
 RFFTforward performs a real FFT. Here, the result is given
 in the two vectors : re and im
------------------------------------------------------*/
void RFFTforward(float *data, int length, float *re, float *im)
{
  int j;
  RFFT(data, length, length, FourierForward);
  re[0] = data[0];
  im[0] = 0;
  re[(length + 2) / 2 - 1] = data[1];
  im[(length + 2) / 2 - 1] = 0;
  for (j = length - 1; j > length / 2; j--) {
    re[j] = data[2 * (length - j)];
    im[j] = data[2 * (length - j) + 1];
    re[length - j] = re[j];
    im[length - j] = -im[j];
  }
}

/*-----------------------------------------------------
 CFFTforward performs a cplx FFT. Here, the input is given with
 one vector for the real parts and one for the imag ones.
 The results are given the same way.
 ------------------------------------------------------*/
void CFFTforward(float *re, float *im, int length)
{
  float *rec, *imc;
  int j;
  rec = (float *)malloc(length * sizeof(float));
  imc = (float *)malloc(length * sizeof(float));

  RFFT(re, length, length, FourierForward);
  RFFT(im, length, length, FourierForward);

  for (j = 0; j < length; j++) {
    rec[j] = re[j];
    imc[j] = im[j];
  }
  re[0] = rec[0];
  im[0] = imc[0];
  re[(length + 2) / 2 - 1] = rec[1];
  im[(length + 2) / 2 - 1] = imc[1];

  for (j = length - 1; j > length / 2; j--) {
    re[j] = rec[2 * (length - j)] - imc[2 * (length - j) + 1];
    im[j] = rec[2 * (length - j) + 1] + imc[2 * (length - j)];

    re[length - j] = rec[2 * (length - j)] + imc[2 * (length - j) + 1];
    im[length - j] = -rec[2 * (length - j) + 1] + imc[2 * (length - j)];
  }
  free(rec);
  free(imc);
}

/*-----------------------------------------------------
 CFFTbackward performs a cplx FFT inverse .
 ------------------------------------------------------*/
void CFFTbackward(float *re, float *im, int length)
{
  float *a, *b;
  int j;
  a = (float *)malloc(length * sizeof(float));
  b = (float *)malloc(length * sizeof(float));

  // fft(a+ib) -> fft(a) , fft(b)
  //  re, im   ->   a        b
  // then         ifft(a)   ifft(b)
  //               re         im

  a[0] = re[0];
  b[0] = im[0];
  a[1] = re[(length + 2) / 2 - 1];
  b[1] = im[(length + 2) / 2 - 1];

  for (j = length - 1; j > length / 2; j--) {
    a[2 * (length - j)] = 0.5f * (re[j] + re[length - j]);
    a[2 * (length - j) + 1] = 0.5f * (im[j] - im[length - j]);

    b[2 * (length - j)] = 0.5f * (im[j] + im[length - j]);
    b[2 * (length - j) + 1] = 0.5f * (-re[j] + re[length - j]);
  }

  RFFT(a, length, length, FourierBackward);
  RFFT(b, length, length, FourierBackward);

  for (j = 0; j < length; j++) {
    re[j] = a[j] / (float)length * 2.0f;
    im[j] = b[j] / (float)length * 2.0f;
  }
  free(a);
  free(b);
}

/*-----------------------------------------------------
 switch_with_z  switch the z coords with either the x (if
 is_y == 0) or the y (is_y ==1)one of vect.
 This function is used to perform a FFT of a 3D image with
 only 1D-FFT. We perform one 1D-FFT, then we swith 2 coords
 and perform another 1D-FFT and we repeat this process once.
 ------------------------------------------------------*/
void FFTswitch_with_z(float ***vect, int dimension, int is_y)
{
  float ***res;
  int x, y, z;
  res = (float ***)malloc(dimension * sizeof(float **));

  for (z = 0; z < dimension; z++) {
    res[z] = (float **)malloc(dimension * sizeof(float *));
    for (y = 0; y < dimension; y++) res[z][y] = (float *)malloc(dimension * sizeof(float));
  }
  for (z = 0; z < dimension; z++)
    for (y = 0; y < dimension; y++)
      for (x = 0; x < dimension; x++) {
        if (is_y)
          res[x][y][z] = vect[x][z][y];
        else
          res[x][y][z] = vect[z][y][x];
      }
  for (z = 0; z < dimension; z++)
    for (y = 0; y < dimension; y++)
      for (x = 0; x < dimension; x++) vect[x][y][z] = res[x][y][z];
  free(res);
}

float ***FFTinv_quarter(float ***vect, int dimension)
{
  int transl = dimension / 2;
  int x, y, z, k, j;
  float ***res;
  res = (float ***)malloc(dimension * sizeof(float **));
  for (k = 0; k < dimension; k++) {
    res[k] = (float **)malloc(dimension * sizeof(float *));
    for (j = 0; j < dimension; j++) {
      res[k][j] = (float *)malloc(dimension * sizeof(float));
    }
  }
  for (z = 0; z < transl; z++) {
    for (y = 0; y < transl; y++) {
      for (x = 0; x < transl; x++) res[x][y][z] = vect[x + transl][y + transl][z + transl];
      for (x = transl; x < dimension; x++) res[x][y][z] = vect[x - transl][y + transl][z + transl];
    }
    for (y = transl; y < dimension; y++) {
      for (x = 0; x < transl; x++) res[x][y][z] = vect[x + transl][y - transl][z + transl];
      for (x = transl; x < dimension; x++) res[x][y][z] = vect[x - transl][y - transl][z + transl];
    }
  }
  for (z = transl; z < dimension; z++) {
    for (y = 0; y < transl; y++) {
      for (x = 0; x < transl; x++) res[x][y][z] = vect[x + transl][y + transl][z - transl];
      for (x = transl; x < dimension; x++) res[x][y][z] = vect[x - transl][y + transl][z - transl];
    }
    for (y = transl; y < dimension; y++) {
      for (x = 0; x < transl; x++) res[x][y][z] = vect[x + transl][y - transl][z - transl];
      for (x = transl; x < dimension; x++) res[x][y][z] = vect[x - transl][y - transl][z - transl];
    }
  }
  return (res);
}
float argument(float re, float im)
{
  if (re == 0 && im == 0) return 0;
  if (re > 0) return (atan(im / re));
  if (re < 0) {
    if (im >= 0)
      return (atan(im / re) + M_PI);
    else
      return (atan(im / re) - M_PI);
  }
  if (im > 0)  // re == 0
    return (M_PI / 2);
  else  // im<0
    return (-M_PI / 2);
}

void FFTreim_to_modarg(float ***re_mod, float ***im_arg, int l)
{
  int x, y, z;
  float a, b;
  for (z = 0; z < l; z++)
    for (y = 0; y < l; y++)
      for (x = 0; x < l; x++) {
        a = re_mod[x][y][z];
        b = im_arg[x][y][z];
        re_mod[x][y][z] = sqrt(a * a + b * b);
        im_arg[x][y][z] = argument(a, b);
      }
}

void FFTmodarg_to_reim(float ***re_mod, float ***im_arg, int l)
{
  int x, y, z;
  float a, b;
  for (z = 0; z < l; z++)
    for (y = 0; y < l; y++)
      for (x = 0; x < l; x++) {
        a = re_mod[x][y][z];
        b = im_arg[x][y][z];
        re_mod[x][y][z] = a * cos(b);
        im_arg[x][y][z] = a * sin(b);
      }
}

float FFTdist(int x, int y, int z, int len)
{
  return ((x - (float)len / 2) * (x - (float)len / 2) + (y - (float)len / 2) * (y - (float)len / 2) +
          (z - (float)len / 2) * (z - (float)len / 2));
}
