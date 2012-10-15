/**
 * @file  kvlAtlasMeshSmoother.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:39 $
 *    $Revision: 1.3 $
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
 */
#ifndef __kvlAtlasMeshSmoother_h
#define __kvlAtlasMeshSmoother_h

#include "kvlAtlasMeshCollection.h"


namespace kvl
{


class AtlasMeshSmoother: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshSmoother  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshSmoother, itk::Object );

  //
  void SetMeshCollection( AtlasMeshCollection* meshCollection )
  {
    m_MeshCollection = meshCollection;
  }

  //
  const AtlasMeshCollection* GetMeshCollection() const
  {
    return m_MeshCollection;
  }

  //
  void  SetSigma( float sigma )
  {
    m_Sigma = sigma;
  }

  //
  float  GetSigma() const
  {
    return m_Sigma;
  }

  //
  AtlasMeshCollection::Pointer  GetSmoothedMeshCollection();


protected :
  // Constructor
  AtlasMeshSmoother();

  // Destructor
  virtual ~AtlasMeshSmoother();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;



private :
  AtlasMeshSmoother(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  // Data members
  AtlasMeshCollection::Pointer  m_MeshCollection;
  float  m_Sigma;

};



} // end namespace kvl


#endif
