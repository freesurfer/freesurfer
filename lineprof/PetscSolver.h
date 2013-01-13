
#ifndef _Solver_h
#define _Solver_h

#include <list>

#include <itkImage.h>
#include <itkImageRegionIterator.h>
#include <itkConstNeighborhoodIterator.h>
#include <itkNeighborhoodIterator.h>

#define USE_PETSC 1

#if USE_PETSC
#include "petscksp.h"
#endif

//--------------------------------------------------------

class PetscSolver
{
 public:
  PetscSolver();
  ~PetscSolver();
  int Update(double);

  static const unsigned int Dimension = 2;
  typedef unsigned char MaskPixelType;
  typedef itk::Image<MaskPixelType,Dimension> MaskImageType;
  typedef MaskImageType::Pointer MaskImagePointer;
  typedef double PixelType;
  typedef itk::Image<PixelType, Dimension> ImageType;
  typedef ImageType::Pointer ImagePointer;

  void SetInputData( ImagePointer inputData)
  { data = inputData; }
  void SetInputMask( MaskImagePointer inputMask,
		     MaskPixelType insideValue,
		     MaskPixelType zeroValue );
  
  inline ImagePointer GetOutputData() { return this->data;}
  inline MaskImagePointer GetOutputMask() { return this->mask; }

  inline MaskPixelType GetInsideValue() const { return labelInside; }
  inline MaskPixelType GetZeroValue() const { return labelZero; }


 private:
  ImagePointer data;
  MaskImagePointer mask;
  
  MaskPixelType labelInside;
  MaskPixelType labelZero;

  typedef itk::ImageRegionIterator<ImageType> Iterator;
  typedef itk::ImageRegionConstIterator<ImageType> ConstIterator;

  typedef itk::NeighborhoodIterator<ImageType> NeighborhoodIterator;
  typedef itk::ConstNeighborhoodIterator<ImageType> ConstNeighborhoodIterator;

  typedef itk::ImageRegionConstIterator<MaskImageType> MaskConstIterator;

  int SetupIndexCorrespondence();
  void SetupSystem();
  int DistributeSolution();

  Vec sol, rhs;
  Mat A;
  typedef unsigned int IndexPixelType;
  typedef itk::Image<IndexPixelType,Dimension> IndexImageType;
  typedef IndexImageType::Pointer IndexImagePointer;
  IndexImagePointer indexImage;

  typedef itk::ImageRegionConstIterator<MaskImageType>
    MaskConstIteratorType;
  typedef itk::ImageRegionIterator<IndexImageType>
    IndexIterator;
  typedef itk::ImageRegionConstIterator<IndexImageType>
    IndexConstIterator;

};

#endif
