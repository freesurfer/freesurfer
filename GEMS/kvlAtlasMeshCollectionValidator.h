/**
 * @file  kvlAtlasMeshCollectionValidator.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Koen Van Leemput
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2012/10/15 21:17:38 $
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
#ifndef __kvlAtlasMeshCollectionValidator_h
#define __kvlAtlasMeshCollectionValidator_h

#include "kvlAtlasMeshCollection.h"


namespace kvl
{


class AtlasMeshCollectionValidator: public itk::Object
{
public :

  /** Standard class typedefs */
  typedef AtlasMeshCollectionValidator  Self;
  typedef itk::Object  Superclass;
  typedef itk::SmartPointer< Self >  Pointer;
  typedef itk::SmartPointer< const Self >  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro( Self );

  /** Run-time type information (and related methods). */
  itkTypeMacro( AtlasMeshCollectionValidator, itk::Object );

  //
  bool Validate( const AtlasMeshCollection* meshCollection );

protected :
  // Constructor
  AtlasMeshCollectionValidator();

  // Destructor
  virtual ~AtlasMeshCollectionValidator();

  // Print
  void PrintSelf( std::ostream& os, itk::Indent indent ) const;



private :
  AtlasMeshCollectionValidator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

};



} // end namespace kvl


#endif
