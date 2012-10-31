/**
 * @file  FSSurface.h
 * @brief Base surface class that takes care of I/O and data conversion.
 *
 */
/*
 * Original Author: Ruopeng Wang
 * CVS Revision Info:
 *    $Author: rpwang $
 *    $Date: 2012/10/31 20:10:11 $
 *    $Revision: 1.40 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 *
 */

#ifndef FSSurface_h
#define FSSurface_h

#include <QObject>
#include "vtkSmartPointer.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkMatrix4x4.h"

#include <vector>
#include <string>

extern "C"
{
#include "mrisurf.h"
#include "mri.h"
}

#define NUM_OF_VSETS 5

class vtkTransform;
class FSVolume;

class FSSurface : public QObject
{
  Q_OBJECT
public:
  FSSurface( FSVolume* ref, QObject* parent = NULL );
  virtual ~FSSurface();

  enum ACTIVE_SURFACE { SurfaceMain = 0, SurfaceInflated, SurfaceWhite, SurfacePial, SurfaceOriginal };

  bool MRISRead( const QString& filename,
                 const QString& vector_filename = QString(),
                 const QString& patch_filename = QString(),
                 const QString& target_filename = QString(),
                 bool bAllSurfaces = false);

  bool MRISWrite( const QString& filename );

  bool MRISReadVectors( const QString& filename );

  void GetBounds ( float oRASBounds[6] );

  void ConvertSurfaceToRAS ( float iSurfX, float iSurfY, float iSurfZ,
                             float& oRASX, float& oRASY, float& oRASZ ) const;
  void ConvertSurfaceToRAS ( double iSurfX, double iSurfY, double iSurfZ,
                             double& oRASX, double& oRASY, double& oRASZ ) const;
  void ConvertRASToSurface ( float iRASX, float iRASY, float iRASZ,
                             float& oSurfX, float& oSurfY, float& oSurfZ) const;
  void ConvertRASToSurface ( double iRASX, double iRASY, double iRASZ,
                             double& oSurfX, double& oSurfY, double& oSurfZ) const;
  void ConvertSurfaceToRAS ( float const iSurf[3], float oRAS[3] ) const;
  void ConvertSurfaceToRAS ( double const iSurf[3], double oRAS[3] ) const;
  void ConvertRASToSurface ( float const iRAS[3], float oSurf[3] ) const;
  void ConvertRASToSurface ( double const iRAS[3], double oSurf[3] ) const;

  // Description:
  // Get the vertex number from a RAS or surface RAS point. This uses
  // the hash table and finds only the closest vertex point. If
  // oDistance is not NULL, the distance to the found point will be
  // returned there.
  int FindVertexAtRAS        ( float  const iRAS[3],       float*  oDistance );
  int FindVertexAtRAS        ( double const iRAS[3],       double* oDistance );
  int FindVertexAtSurfaceRAS ( float  const iSurfaceRAS[3],float*  oDistance );
  int FindVertexAtSurfaceRAS ( double const iSurfaceRAS[3],double* oDistance );

  // Description:
  // Get the RAS or surface RAS coords at a vertex index.
  bool GetRASAtVertex        ( int inVertex, float  ioRAS[3] );
  bool GetRASAtVertex        ( int inVertex, double ioRAS[3] );
  bool GetSurfaceRASAtVertex ( int inVertex, float  ioRAS[3] );
  bool GetSurfaceRASAtVertex ( int inVertex, double ioRAS[3] );

  int GetNumberOfVertices () const;

  bool LoadSurface    ( const QString& filename, int nSet );
  bool LoadCurvature  ( const QString& filename = NULL );
  bool LoadOverlay    ( const QString& filename, const QString& fn_reg,
                        float** data_out, int* nvertices_out, int* nframes_out );

  bool IsSurfaceLoaded( int nSet )
  {
    return m_bSurfaceLoaded[nSet];
  }

  bool IsCurvatureLoaded()
  {
    return m_bCurvatureLoaded;
  }

  double GetCurvatureValue( int nVertex );

  bool SetActiveSurface( int nIndex );

  int GetActiveSurface()
  {
    return m_nActiveSurface;
  }

  bool HasVectorSet()
  {
    return m_vertexVectors.size() > 0;
  }

  int GetNumberOfVectorSets()
  {
    return m_vertexVectors.size();
  }

  QString GetVectorSetName( int nSet );

  int GetActiveVector()
  {
    return m_nActiveVector;
  }

  bool SetActiveVector( int nIndex );

  vtkPolyData* GetPolyData()
  {
    return m_polydata;
  }

  vtkPolyData* GetTargetPolyData()
  {
    return m_polydataTarget;
  }

  vtkPolyData* GetVectorPolyData()
  {
    return m_polydataVector;
  }

  vtkPolyData* GetVector2DPolyData( int n )
  {
    return m_polydataVector2D[n];
  }

  vtkPolyData* GetVertexPolyData()
  {
    return m_polydataVertices;
  }

  void GetVectorAtVertex( int nVertex, double* vec_out, int nVector = -1 );

  vtkPolyData* GetWireframePolyData()
  {
    return m_polydataWireframes;
  }

  MRIS* GetMRIS()
  {
    return m_MRIS;
  }

  void GetNormalAtVertex( int nVertex, double* vec_out );

  void UpdateVector2D( int nPlane, double slice_pos,
                       vtkPolyData* contour_polydata = NULL );

  void Reposition( FSVolume* volume, int target_vnos, double target_val, int nsize, double sigma, int flags = 0 );

  void Reposition( FSVolume* volume, int target_vnos, double* coord, int nsize, double sigma, int flags = 0 );

  void RepositionVertex( int vno, double* coord );

  bool Smooth(int nMethod, int niters, double lambda, double k_cutoff);

  void RemoveIntersections();

  void UndoReposition();

  bool HasValidVolumeGeometry()
  {
    return m_bValidVolumeGeometry;
  }

protected:
  void UpdatePolyData();
  void UpdatePolyData( MRIS* mris, vtkPolyData* polydata,
                       vtkPolyData* polydata_verts = NULL,
                       vtkPolyData* polydata_wireframe = NULL );
  void UpdateVerticesAndNormals();
  void ComputeNormals();
  void NormalFace(int fac, int n, float *norm );
  float TriangleArea( int fac, int n );
  void Normalize( float v[3] );

  bool LoadVectors( const QString& filename );
  void LoadTargetSurface( const QString& filename );
  void UpdateVectors();
  void UpdateVertices();

  void SaveNormals ( MRIS* mris, int nSet );
  void RestoreNormals ( MRIS* mris, int nSet );
  void SaveVertices ( MRIS* mris, int nSet );
  void RestoreVertices( MRIS* mris, int nSet );

  bool ProjectVectorPoint2D( double* pt_in,
                             vtkPoints* contour_pts,
                             vtkCellArray* contour_lines,
                             double* pt_out );

  void PostEditProcess();

  MRIS*   m_MRIS;
  MRIS*   m_MRISTarget;

  double  m_SurfaceToRASMatrix[16];
  vtkSmartPointer<vtkTransform> m_SurfaceToRASTransform;

  // RAS bounds.
  bool    m_bBoundsCacheDirty;
  float   m_RASBounds[6];
  float   m_RASCenter[3];

  vtkSmartPointer<vtkPolyData> m_polydata;
  vtkSmartPointer<vtkPolyData> m_polydataVector;
  vtkSmartPointer<vtkPolyData> m_polydataVertices;
  vtkSmartPointer<vtkPolyData> m_polydataWireframes;
  vtkSmartPointer<vtkPolyData> m_polydataVector2D[3];
  vtkSmartPointer<vtkPolyData> m_polydataVertex2D[3];
  vtkSmartPointer<vtkPolyData> m_polydataTarget;

  // Hash table so we can look up vertices. Uses v->x,y,z.
  MRIS_HASH_TABLE* m_HashTable[5];

  bool m_bSurfaceLoaded[NUM_OF_VSETS];
  bool m_bCurvatureLoaded;

  int  m_nActiveSurface;

  FSVolume* m_volumeRef;

  struct VertexItem
  {
    float x;
    float y;
    float z;
  };

  bool SaveVertices( MRIS* mris, VertexItem*& pVertex );
  bool SaveVertices( MRI* mri,   VertexItem*& pVertex );
  bool ComputeVectors( MRIS* mris, VertexItem*& buffer );

  VertexItem*  m_fVertexSets[NUM_OF_VSETS];
  VertexItem*  m_fNormalSets[NUM_OF_VSETS];

  struct VertexVectorItem
  {
    QString  name;
    VertexItem*  data;

    VertexVectorItem()
    {
      data = NULL;
    }
  };

  std::vector<VertexVectorItem>  m_vertexVectors;
  int      m_nActiveVector;

  bool     m_bValidVolumeGeometry;
};

#endif


