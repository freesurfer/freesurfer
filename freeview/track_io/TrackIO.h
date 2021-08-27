/**
 * @brief I/O interface for track file (.trk)
 *
 */
/*
 * Original Author: Ruopeng Wang
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

///////////////////////////////////////////////////////////////////////////////
//
// Usage:
//      Sample lines to read track file -
//
//          include "TrackIO.h"
//
//      ...
//
//      CTrackReader reader;
//          TRACK_HEADER header;
//      if (!reader.Open("foo.trk", &header))
//      {
//        printf(reader.GetLastErrorMessage());
//        return;
//      }
//      ...
//
//      int cnt;
//      if (ignore_scalars && ignore_properties)
//      {
//        while (reader.GetNextPointCount(&cnt))
//        {
//          float* pts = new float[cnt*3];
//          reader.GetNextTrackData(cnt, pts);
//          ...
//                process_point_data(...);
//          ...
//          delete[] pts;
//        }
//      }
//      else
//      {
//        while (reader.GetNextPointCount(&cnt))
//        {
//          float* pts = new float[cnt*3];
//          float* scalars = new float[cnt*header.n_scalars];
//          float* properties = new float[header.n_properties];
//          reader.GetNextTrackData(cnt, pts, scalars, properties);
//          ...
//                process_point_and_scalar_data_etc.(...);
//          ...
//          delete[] pts;
//          delete[] scalars;
//          delete[] properties;
//        }
//      }
//
//      reader.Close();
//
///////////////////////////////////////////////////////////////////////////////

#ifndef _TrackIO_H_
#define _TrackIO_H_

#include <stdio.h>
//#include "Common.h"
#include "ByteSwap.h"
#include "ErrorCode.h"
#include <string.h>
#include <vector>

#ifndef DEFAULT_VOXEL_ORDER
#define DEFAULT_VOXEL_ORDER "LPS"
#endif

#ifndef HEADER_VERSION
#define HEADER_VERSION    2
#endif

struct TRACK_HEADER
{
  char      id_string[6]; // first 5 chars must be "TRACK"
  short int   dim[3];     // dimensions
  float     voxel_size[3];  // voxel size
  float     origin[3];    // origin. default are 0,0,0.
  short int   n_scalars;    // number of scalars saved per point besides xyz coordinates.
  char      scalar_name[10][20]; // name of the scalars

  short int   n_properties; // number of properties
  char      property_name[10][20]; // name of the properties

  float     vox_to_ras[4][4]; // voxel to ras (ijk to xyz) matrix, this is used for coordinate transformation
  // if vox_to_ras[3][3] is 0, it means v2r matrix is not recorded
  // this field is added from version 2.
  char      reserved[444];
  char      voxel_order[4]; // voxel order for this track space
  // if there was no reorientation, this should be the same
  // as voxel_order_original
  char      voxel_order_original[4];
  // voxel order of the original image data
  float     image_orientation_patient[6];
  // image orientation patient info as usually recorded in dicom tags
  // this will help display program to determine the correct
  // image orientation
  char      pad1[2];
  unsigned char invert_x;   // inversion/rotation flags used to generate this track file
  unsigned char invert_y;   // value is 0 or 1. can be ignored (for private use only).
  unsigned char invert_z;
  unsigned char swap_xy;
  unsigned char swap_yz;
  unsigned char swap_zx;
  int       n_count;    // total number of tracks. if 0, number of tracks was not recorded.
  // call GetNumberOfTracks(...) to get it
  int       version;    // version number
  int       hdr_size;   // size of the header. used to determine byte swap

  TRACK_HEADER()
  {
    Initialize();
  }

  TRACK_HEADER(int* d, float* vs, float* o, int n = 0)
  {
    Initialize();

    for (int i = 0; i < 3; i++)
    {
      dim[i] = d[i];
      voxel_size[i] = vs[i];
      origin[i] = o[i];
    }
    n_scalars = n;
  }

  TRACK_HEADER(short int* d, float* vs, float* o, int n = 0)
  {
    Initialize();

    for (int i = 0; i < 3; i++)
    {
      dim[i] = d[i];
      voxel_size[i] = vs[i];
      origin[i] = o[i];
    }
    n_scalars = n;
  }

  inline void ByteSwap()
  {
    SWAP_SHORT(dim, 3);
    SWAP_FLOAT(voxel_size, 3);
    SWAP_FLOAT(origin, 3);
    SWAP_SHORT(n_scalars);
    SWAP_SHORT(n_properties);
    SWAP_FLOAT(image_orientation_patient, 6);
    for ( int i = 0; i < 4; i++ )
    {
      SWAP_FLOAT(vox_to_ras[i], 4);
    }
    SWAP_INT(version);
    SWAP_INT(n_count);
    SWAP_INT(hdr_size);
  }

  inline void Initialize()
  {
    memset(this, 0, sizeof(struct TRACK_HEADER));
    strcpy(id_string, "TRACK");
    hdr_size = sizeof(struct TRACK_HEADER);

    // default image orientation is axial, voxel order is LPS
    image_orientation_patient[0] = 1;
    image_orientation_patient[4] = 1;
    strcpy(voxel_order, DEFAULT_VOXEL_ORDER);
//    strcpy(voxel_order_original, DEFAULT_VOXEL_ORDER);

    version = HEADER_VERSION;
  }

};


class CTrackIO
{
public:
  CTrackIO()
  {
    m_nErrorCode = 0;
    m_pFile = 0;
  }
  virtual ~CTrackIO()
  {
    Close();
  }

  bool GetHeader(TRACK_HEADER* header);
  virtual bool Close();
  const char* GetLastErrorMessage();
  int GetLastErrorCode()
  {
    return m_nErrorCode;
  }

protected:
  TRACK_HEADER  m_header;
  FILE* m_pFile;

  int   m_nErrorCode;
};

class CTrackReader : public CTrackIO
{
public:
  CTrackReader();

  bool Open(const char* filename, TRACK_HEADER* header = NULL);
  bool Open(const char* filename, int* dim, float* voxel_size);
  bool GetNextPointCount(int* ncount);
  bool GetNextRawData(int ncount, float* data);
  bool GetNextTrackData(int nCount, float* pt_data, float* scalars = NULL, float* properties = NULL);
  int GetProgress();
  int  GetNumberOfTracks();
  bool GetNumberOfTracks(int* cnt);
  bool ByteSwapped()
  {
    return m_bByteSwap;
  }
  bool IsOldFormat()
  {
    return m_bOldFormat;
  }
  void AllowOldFormat(bool bAllow)
  {
    m_bAllowOldFormat = bAllow;
  }

  static bool GetHeader(const char* filename, TRACK_HEADER* header);

protected:
  bool      m_bByteSwap;
  bool      m_bOldFormat;
  long      m_nSize;
  bool      m_bAllowOldFormat;
};

class CTrackWriter : public CTrackIO
{
public:
  bool Initialize(const char* filename, short int* dim, float* voxel_size, float* origin = NULL,
                  short int n_scalars = 0);
  bool Initialize(const char* filename, int* dim, float* voxel_size, float* origin = NULL,
                  short int n_scalars = 0);
  bool Initialize(const char* filename, TRACK_HEADER header);
  bool WriteNextTrack(int ncount, float* data);
  bool WriteNextTrack(int ncount, float* pts, float* scalars, float* properties);
  bool UpdateHeader(TRACK_HEADER header);

  virtual bool Close();

  static bool UpdateHeader(const char* filename, TRACK_HEADER header);
};

#endif
