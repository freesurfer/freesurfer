
// ITK includes
#include <itkImageRegionConstIterator.h>

#include "PetscSolver.h"

PetscSolver::PetscSolver()
{
  data = NULL;
  mask = NULL;

  sol = 0;
  rhs = 0;
  A = 0;
  indexImage = NULL;
}

PetscSolver::~PetscSolver()
{
  if ( A )   MatDestroy(A);
  if ( rhs ) VecDestroy( rhs );
  if ( sol ) VecDestroy( sol );
}

void
PetscSolver::SetInputMask(MaskImagePointer inputMask,
        MaskPixelType insideValue,
        MaskPixelType zeroValue)
{
  mask = inputMask;
  labelInside = insideValue;
  labelZero = zeroValue;
}


int
PetscSolver::Update(double convergence)
{
  PetscErrorCode ierr; 

  int pixelCount = this->SetupIndexCorrespondence();
  std::cout << " Linear System size = " << pixelCount << std::endl;

  ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
  ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE,
         pixelCount, pixelCount); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&rhs); CHKERRQ(ierr);
  ierr = VecSetSizes(rhs, PETSC_DECIDE, pixelCount); CHKERRQ(ierr);
  ierr = VecSetFromOptions(rhs); CHKERRQ(ierr);
  ierr = VecDuplicate(rhs,&sol); CHKERRQ(ierr);

  this->SetupSystem();
  // solve the system
  KSP ksp;
  ierr = KSPCreate( PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
  PC pc;
  ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
  ierr = PCSetType(pc, PCJACOBI); CHKERRQ(ierr);
  double tol = pow(10.0, -convergence);
  ierr = KSPSetTolerances(ksp, tol/pixelCount, 
        PETSC_DEFAULT,
        PETSC_DEFAULT,
        pixelCount/3); CHKERRQ(ierr);
  ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

  ierr = KSPSolve(ksp, rhs, sol); CHKERRQ(ierr);
  KSPConvergedReason reason;
  ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
  this->DistributeSolution();

  ierr = KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD); CHKERRQ(ierr);
  PetscInt its;
  ierr = KSPGetIterationNumber( ksp, &its); CHKERRQ(ierr);
  std::cout << " Solver iterations = " << its << std::endl;

  Vec vcheck;
  ierr = VecDuplicate( sol, &vcheck ); CHKERRQ(ierr);
  ierr = MatMult( A, sol, vcheck ); CHKERRQ(ierr);
  PetscScalar neg_one = -1.0;
  PetscReal norm;
  ierr = VecAXPY( vcheck, neg_one, rhs ); CHKERRQ(ierr);
  ierr = VecNorm(vcheck, NORM_2, &norm); CHKERRQ(ierr);
  ierr = VecDestroy(vcheck); CHKERRQ(ierr);


  ierr = KSPDestroy(ksp); CHKERRQ(ierr);

  return 0;
}

int
PetscSolver::SetupIndexCorrespondence()
{
  // allocate index image
  indexImage = IndexImageType::New();
  indexImage->SetRegions( data->GetRequestedRegion() );
  indexImage->Allocate();

  MaskConstIteratorType maskCiter( mask, mask->GetRequestedRegion() );
  IndexIterator indexIter( indexImage, indexImage->GetRequestedRegion() );
  
  unsigned int count = 0;
  for( maskCiter.GoToBegin(), indexIter.GoToBegin();
       !indexIter.IsAtEnd(); 
       ++indexIter, ++maskCiter )
  {
    if ( maskCiter.Get() == labelInside ) indexIter.Set(count++);
  } // next indexIter, maskCiter
  return count;
}

void
PetscSolver::SetupSystem()
{
  typedef itk::ConstNeighborhoodIterator<MaskImageType> MaskCNIterator;
  typedef itk::ConstNeighborhoodIterator<IndexImageType> IndexCNIterator;

  MaskCNIterator::RadiusType radius;
  radius.Fill(1);

  MaskCNIterator mask_cni( radius, mask, mask->GetRequestedRegion() );
  IndexCNIterator index_cni( radius, indexImage, indexImage->GetRequestedRegion() );
  ConstNeighborhoodIterator data_cni( radius, data, data->GetRequestedRegion() );
         
  NeighborhoodIterator::OffsetType offset_mx;
  NeighborhoodIterator::OffsetType offset_px;
  NeighborhoodIterator::OffsetType offset_my;
  NeighborhoodIterator::OffsetType offset_py;
  NeighborhoodIterator::OffsetType offset_zero;
  offset_mx[0] = -1; offset_mx[1] = 0;
  offset_px[0] = 1;  offset_px[1] = 0;
  offset_my[0] = 0;  offset_my[1] = -1;
  offset_py[0] = 0;  offset_py[1] = 1;
  offset_zero[0]=0;  offset_zero[1]=0;

  typedef std::vector<NeighborhoodIterator::OffsetType> OffsetVectorType;
  OffsetVectorType vecOffsets;

  vecOffsets.push_back(offset_mx);
  vecOffsets.push_back(offset_px);
  vecOffsets.push_back(offset_my);
  vecOffsets.push_back(offset_py);

  typedef std::map<int,double> MapType;
  MapType map_rhs;

  OffsetVectorType::const_iterator off_cit;
  int pj[5], counter=0;
  int row;
  PetscScalar value[5];
  //int dbgCount = 0;
  for( index_cni.GoToBegin(), mask_cni.GoToBegin(), data_cni.GoToBegin();
       !index_cni.IsAtEnd();
       ++index_cni, ++mask_cni, ++data_cni, counter=0 )
  { 
    if ( mask_cni.GetPixel(offset_zero) == labelInside )
    {
      row = index_cni.GetPixel(offset_zero);

      pj[counter]=row; 
      value[counter] = 4;
      ++counter;

      for( off_cit = vecOffsets.begin();
           off_cit != vecOffsets.end();
           ++off_cit )
      {
        if ( mask_cni.GetPixel( *off_cit ) == labelInside )
        {
          pj[counter] = index_cni.GetPixel( *off_cit );
          value[counter] = -1;
          ++counter;
        }
        else map_rhs[ row ] += data_cni.GetPixel( *off_cit );
      } // next off_cit

      MatSetValues( A, 1, &row, counter, pj, value, INSERT_VALUES );
  }
    } // next indexCiter, maskCiter

  // finish setup of the matrix
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  // collect values from RHS
  int *rhsIndices = new int[ map_rhs.size() ];
  double *rhsValues = new double[ map_rhs.size() ];
  counter = 0;
  for( MapType::const_iterator map_cit = map_rhs.begin();
       map_cit != map_rhs.end();
       ++map_cit, ++counter )
  {
    rhsIndices[counter] = map_cit->first;
    rhsValues[counter]  = map_cit->second;
  }
  VecSetValues(rhs, counter, rhsIndices, rhsValues, INSERT_VALUES);
  VecAssemblyBegin(rhs);
  VecAssemblyEnd(rhs);
  // cleanup
  delete[] rhsIndices;
  delete[] rhsValues;

}

int
PetscSolver::DistributeSolution()
{
  // quick & dirty - use input image
  Iterator iter(data, data->GetRequestedRegion() );
  MaskConstIterator mask_cit(mask, mask->GetRequestedRegion() );
  IndexConstIterator index_cit(indexImage, indexImage->GetRequestedRegion() );

  // borrow array from the rhs
  PetscScalar* array;
  PetscErrorCode ierr;
  ierr = VecGetArray(sol, &array); CHKERRQ(ierr);

  for( iter.GoToBegin(), mask_cit.GoToBegin(), index_cit.GoToBegin();
       !iter.IsAtEnd();
       ++iter, ++mask_cit, ++index_cit )
  {
    if ( mask_cit.Get() == labelInside )
    {
      iter.Set( array[ index_cit.Get() ] );
    }
  } // next iter, mask_cit, index_cit
  
  ierr = VecRestoreArray(sol, &array); CHKERRQ(ierr);
  return 0;
}
