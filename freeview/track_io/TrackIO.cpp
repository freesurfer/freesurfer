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

#include "TrackIO.h"

///// CTrackIO reference //////////////////
const char* error_message[] =
{
  "No error",
  "Can not open file",
  "Can not close file",
  "Can not read from file",
  "Can not write to file",
  "Not a compatible track file",
  "I/O not initialized. Call 'Open' or 'Initialize' first"
};

const char* CTrackIO::GetLastErrorMessage()
{
  return error_message[m_nErrorCode];
}

bool CTrackIO::GetHeader(TRACK_HEADER* header)
{
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    return false;
  }
  *header = m_header;

  return true;
}

bool CTrackIO::Close()
{
  bool ret = true;
  if (m_pFile)
  {
    if (fclose(m_pFile) == EOF)
    {
      m_nErrorCode = TE_CAN_NOT_CLOSE;
      ret = false;
    }
    m_pFile = NULL;
  }
//  if (ret)
//    m_nErrorCode = TE_NO_ERROR;

  return ret;
}

/////////////////////////////////////////

///// CTrackReader reference //////////////////

CTrackReader::CTrackReader()
{
  m_bByteSwap = false;
  m_bAllowOldFormat = false;
}


bool CTrackReader::Open(const char* filename, TRACK_HEADER* header)
{
  Close();
  m_nErrorCode = TE_NO_ERROR;
  m_bOldFormat = false;
  m_pFile = fopen(filename, "rb");
  if (!m_pFile)
  {
    m_nErrorCode = TE_CAN_NOT_OPEN;
    return false;
  }

  fseek(m_pFile, 0, SEEK_END);
  m_nSize = ftell(m_pFile);
  fseek(m_pFile, 0, SEEK_SET);

  if (fread(&m_header, sizeof(TRACK_HEADER), 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_NOT_TRACK_FILE;
  }

  if (m_header.voxel_order[0] == 0)
  {
    strcpy(m_header.voxel_order, DEFAULT_VOXEL_ORDER);
  }
  if (m_header.voxel_order_original[0] == 0)
  {
    strcpy(m_header.voxel_order_original, m_header.voxel_order);
  }

  // could be old preliminary format from old MGH data.
  // in most cases ignored
  //////////////////////////////////////
  bool bOldFormat = false;
  bool bBigEndian = IS_BIG_ENDIAN();
  if (!m_nErrorCode && strcmp(m_header.id_string, "TRACK") != 0)
  {
    if (!m_bAllowOldFormat)
    {
      m_nErrorCode = TE_NOT_TRACK_FILE;
      Close();
      return false;
    }
    bOldFormat = true;
    m_header.Initialize();
    fseek(m_pFile, 0, SEEK_SET);
    int dim[3];
    fread(dim, sizeof(int)*3, 1, m_pFile);

    if (bBigEndian)
    {
      SWAP_INT(dim, 3);
    }

    for (int i = 0; i < 3; i++)
    {
      m_header.dim[i] = dim[i];
    }

    if (bBigEndian)
    {
      SWAP_SHORT(m_header.dim, 3);
    }

    fread(m_header.voxel_size, sizeof(float)*3, 1, m_pFile);
    m_header.n_scalars = 0;
    m_header.n_properties = 0;

    m_bByteSwap = bBigEndian;

  }
  else
  {
    if (m_header.hdr_size == 0)
    {

      m_bByteSwap = bBigEndian;
    }
    else
    {
      m_bByteSwap = (m_header.hdr_size != sizeof(struct TRACK_HEADER));
    }
  }
  ///////////////////////////////////////////////////////////

  if (m_bByteSwap)
  {
    m_header.ByteSwap();
    m_header.hdr_size = sizeof(struct TRACK_HEADER);
  }

  /*
  if (bOldFormat)
  {
    if (mmin(m_header.dim, 3) <= 0 || mmin(m_header.voxel_size, 3) <= 0.00001 ||
      mmax(m_header.dim, 3) > 10000 || mmax(m_header.voxel_size, 3) > 100)
      m_nErrorCode = TE_NOT_TRACK_FILE;
  }
  */

  // if no image orientation info, assign axial standard
  if ( m_header.image_orientation_patient[0] == 0 &&
       m_header.image_orientation_patient[1] == 0 &&
       m_header.image_orientation_patient[2] == 0 )
  {
    m_header.image_orientation_patient[0] = m_header.image_orientation_patient[4] = 1;
  }

  if (header)
  {
    *header = m_header;
  }

  m_bOldFormat = bOldFormat;

  if (m_nErrorCode != TE_NO_ERROR)
  {
    Close();
    return false;
  }

  return true;
}


bool CTrackReader::Open(const char* filename, int* dim, float* voxel_size)
{
  bool ret = Open(filename);
  for (int i = 0; i < 3; i++)
  {
    dim[i] = m_header.dim[i];
    voxel_size[i] = m_header.voxel_size[i];
  }

  return ret;
}

// Get progess in percentage
int CTrackReader::GetProgress()
{
  if (m_pFile && m_nSize > 0)
  {
    return (int)(ftell(m_pFile) * 100.0 / m_nSize);
  }
  return 0;
}

bool CTrackReader::GetNextPointCount(int* n)
{
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    n = 0;
    return false;
  }
  if (fread(n, sizeof(int), 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_CAN_NOT_READ;
  }
  else
  {
    m_nErrorCode = TE_NO_ERROR;
  }

  if (m_bByteSwap)
  {
    SWAP_INT(*n);
  }

  return m_nErrorCode == TE_NO_ERROR;
}


// Attention: data buffer has to been pre-allocated. The following two routines
// do not allocate memory.
// GetNextRawData(...) read in together track properties, point coordinates and scalars, etc.
// GetNextTrackData(...) read in track properties, point coords and scalars in seperated buffers
// if scalars buffer is not given, only point coords are read.
// These two routines can only be called once after calling GetNextPointCount().
bool CTrackReader::GetNextRawData(int nCount, float* data)
{
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    return false;
  }
  int nSize = (3+m_header.n_scalars)*nCount + m_header.n_properties;
  if (fread(data, sizeof(float)*nSize, 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_CAN_NOT_READ;
  }
  else
  {
    m_nErrorCode = TE_NO_ERROR;
  }

  if (!m_nErrorCode && m_bByteSwap)
  {
    SWAP_FLOAT(data, nSize);
  }

  return m_nErrorCode == TE_NO_ERROR;
}

bool CTrackReader::GetNextTrackData(int nCount, float* pt_data, float* scalars, float* properties)
{
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    return false;
  }
  m_nErrorCode = TE_NO_ERROR;

  if (m_header.n_scalars == 0 && m_header.n_properties == 0)
  {
    return GetNextRawData(nCount, pt_data);
  }

  // read point and scalar data
  for (int i = 0; i < nCount; i++)
  {
    if (fread(pt_data+i*3, sizeof(float)*3, 1, m_pFile) != 1)
    {
      m_nErrorCode = TE_CAN_NOT_READ;
      break;
    }
    if (m_header.n_scalars && scalars)
    {
      if (fread(scalars+i*m_header.n_scalars, sizeof(float)*m_header.n_scalars, 1, m_pFile) != 1)
      {
        m_nErrorCode = TE_CAN_NOT_READ;
        break;
      }
    }
    else
    {
      fseek(m_pFile, sizeof(float)*m_header.n_scalars, SEEK_CUR);
    }
  }
  if (!m_nErrorCode && m_bByteSwap)
  {
    SWAP_FLOAT(pt_data, nCount*3);
    if (scalars)
    {
      SWAP_FLOAT(scalars, nCount*m_header.n_scalars);
    }
  }

  // read property data
  if (m_header.n_properties && properties &&
      fread(properties, sizeof(float)*m_header.n_properties, 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_CAN_NOT_READ;
  }

  if (!m_nErrorCode && m_bByteSwap)
  {
    if (properties)
    {
      SWAP_FLOAT(properties, m_header.n_properties);
    }
  }

  return m_nErrorCode == TE_NO_ERROR;
}

// Static function. Get header info from a given track file directly
bool CTrackReader::GetHeader(const char* filename, TRACK_HEADER *header)
{
  CTrackReader reader;
  bool ret = reader.Open(filename, header);
  reader.Close();

  return ret;
}

// if number of tracks was not recorded in the header, this routine will
// take longer to excute as it will go through the whole file
bool CTrackReader::GetNumberOfTracks(int* cnt)
{
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    *cnt = 0;
    return false;
  }
  if (m_header.n_count == 0)
  {
    long pos = ftell(m_pFile);
    fseek(m_pFile, m_bOldFormat ? 3*(sizeof(int)+sizeof(float)) : sizeof(TRACK_HEADER), SEEK_SET);

    int n;
    *cnt = 0;
    while (GetNextPointCount(&n))
    {
      *cnt += 1;
      fseek(m_pFile, sizeof(float)*(n*(3+m_header.n_scalars)+m_header.n_properties), SEEK_CUR);
    }
    fseek(m_pFile, pos, SEEK_SET);
    m_header.n_count = *cnt;
  }
  else
  {
    *cnt = m_header.n_count;
  }
  m_nErrorCode = TE_NO_ERROR;

  return true;
}

int CTrackReader::GetNumberOfTracks()
{
  int n;
  GetNumberOfTracks(&n);
  return n;
}

///////////////////////////////////////////////

///// CTrackReader reference //////////////////

// One of the Initializers must be called before WriteNextTrackData()
bool CTrackWriter::Initialize(const char* filename, short int* dim, float* voxel_size, float* origin,
                              short int n_scalars)
{
  float org[3] = { 0, 0, 0 };
  if (origin)
  {
    memcpy(org, origin, 3*sizeof(float));
  }

  TRACK_HEADER header(dim, voxel_size, org, n_scalars);

  return Initialize(filename, header);
}

bool CTrackWriter::Initialize(const char* filename, int* dim, float* voxel_size, float* origin,
                              short int n_scalars)
{
  short int ndim[3];
  for (int i = 0; i < 3; i++)
  {
    ndim[i] = dim[i];
  }

  return Initialize(filename, ndim, voxel_size, origin, n_scalars);
}

bool CTrackWriter::Initialize(const char* filename, TRACK_HEADER header)
{
  Close();
  m_pFile = fopen(filename, "wb");
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    return false;
  }

  if (header.hdr_size != sizeof(TRACK_HEADER))
  {
    header.ByteSwap();
  }
  m_header = header;
  m_header.version = HEADER_VERSION;
  m_header.hdr_size = sizeof(struct TRACK_HEADER);

  if (fwrite(&m_header, sizeof(struct TRACK_HEADER), 1, m_pFile) == 1)
  {
    m_nErrorCode = TE_NO_ERROR;
  }
  else
  {
    m_nErrorCode = TE_CAN_NOT_WRITE;
  }

  m_header.n_count = 0;
  return m_nErrorCode == TE_NO_ERROR;
}


// data is raw track data!! Must include scalars if n_scalars is not 0
bool CTrackWriter::WriteNextTrack(int ncount, float* data)
{
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    return false;
  }
  m_nErrorCode = TE_NO_ERROR;

  if (fwrite(&ncount, sizeof(int), 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_CAN_NOT_WRITE;
  }
  long nSize = (ncount*(3+m_header.n_scalars)+m_header.n_properties)*sizeof(float);

  if (fwrite(data, nSize, 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_CAN_NOT_WRITE;
  }

  if (!m_nErrorCode)
  {
    m_header.n_count ++;
  }
  return m_nErrorCode == TE_NO_ERROR;
}

bool CTrackWriter::WriteNextTrack(int ncount, float* pts, float* scalars, float* properties)
{
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    return false;
  }
  m_nErrorCode = TE_NO_ERROR;

  if (fwrite(&ncount, sizeof(int), 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_CAN_NOT_WRITE;
  }

  for (int i = 0; i < ncount; i++)
  {
    if (fwrite(pts+i*3, sizeof(float)*3, 1, m_pFile) != 1)
    {
      m_nErrorCode = TE_CAN_NOT_WRITE;
    }
    if (m_header.n_scalars
        && fwrite(scalars+i*m_header.n_scalars, sizeof(float)*m_header.n_scalars, 1, m_pFile) != 1)
    {
      m_nErrorCode = TE_CAN_NOT_WRITE;
    }
  }

  if (m_header.n_properties && fwrite(properties, sizeof(float)*m_header.n_properties, 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_CAN_NOT_WRITE;
  }

  if (!m_nErrorCode)
  {
    m_header.n_count ++;
  }
  return m_nErrorCode == TE_NO_ERROR;
}


bool CTrackWriter::UpdateHeader(TRACK_HEADER header)
{
  if (!m_pFile)
  {
    m_nErrorCode = TE_NOT_INITIALIZED;
    return false;
  }

  if (header.hdr_size != sizeof(TRACK_HEADER))
  {
    header.ByteSwap();
  }

  m_header = header;
  m_header.hdr_size = sizeof(struct TRACK_HEADER);

  long pos = ftell(m_pFile);
  fseek(m_pFile, 0, SEEK_SET);

  if (fwrite(&m_header, sizeof(struct TRACK_HEADER), 1, m_pFile) != 1)
  {
    m_nErrorCode = TE_CAN_NOT_WRITE;
  }
  else
  {
    m_nErrorCode = TE_NO_ERROR;
  }
  fseek(m_pFile, pos, SEEK_SET);

  return m_nErrorCode == TE_NO_ERROR;
}

bool CTrackWriter::Close()
{
  return UpdateHeader(m_header) && CTrackIO::Close();
}

// Write the given header to a existing track file.
bool CTrackWriter::UpdateHeader(const char* filename, TRACK_HEADER header)
{
  TRACK_HEADER hdr;
  FILE* fp = fopen(filename, "r+b");
  if (!fp)
  {
    return false;
  }

  if (fread(&hdr, sizeof(TRACK_HEADER), 1, fp) != 1)
  {
    fclose(fp);
    return false;
  }

  bool bswap = (hdr.hdr_size != sizeof(struct TRACK_HEADER));
  hdr = header;
  hdr.hdr_size = sizeof(struct TRACK_HEADER);

  if (bswap)
  {
    hdr.ByteSwap();
  }

  fseek(fp, 0, SEEK_SET);
  if (fwrite(&hdr, sizeof(TRACK_HEADER), 1, fp) != 1)
  {
    fclose(fp);
    return false;
  }
  fclose(fp);

  return true;
}

