function [ affineRegistrationMeshCollectionFileName, affineRegistrationTemplateFileName ] = createAtlasMeshForAffineRegistration( priors, transformMatrix, savePath, meshSize )
%  
%  Create a regular (low-res) tetrahedral mesh atlas from existing volumetric prior (uint16).
%  The result is a mesh with the following order of classes: background, WM, GM, CSF, [bone, soft tissue]
%  (last two are optional, depending on input volumetric prior that is given)
%

if ( nargin < 4 )
  meshSize = [ 30 30 30 ];
end

if ( nargin < 1 )
  % Use the priors distributed with SPM12. 
  
  % Read in the correct TPMs extracted manually from the 6-channel NIFTI atlas. 
  DIM = [ 121 145 121 ];
  numberOfClasses = 6;
  % numberOfClasses = 4;  
  priors = zeros( [ DIM numberOfClasses ] );
  
  priors( :, :, :, 1 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/background.nii' ) );
  priors( :, :, :, 2 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/wm.nii' ) );
  priors( :, :, :, 3 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/gm.nii' ) );
  priors( :, :, :, 4 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/csf.nii' ) );
  if ( numberOfClasses > 4 )
    priors( :, :, :, 5 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/bone.nii' ) );
    priors( :, :, :, 6 ) = kvlGetImageBuffer( kvlReadImage( '/data/testing/atlas/soft.nii' ) );
  end  
  
  % Also get the image-to-world transform of this template
  [ ~, transform ] = kvlReadImage( '/data/testing/atlas/wm.nii' );
  transformMatrix = kvlGetTransformMatrix( transform ); 

  
  % Add a rim of background
  if 0
    numberOfPaddingVoxels = 5;
    paddedDIM = DIM + 2 * numberOfPaddingVoxels;
    paddedPriors = zeros( [ paddedDIM numberOfClasses ] );
    paddedPriors( :, :, :, 1 ) = 1.0;
    paddedPriors( numberOfPaddingVoxels + [ 1 : DIM(1) ], ...
                  numberOfPaddingVoxels + [ 1 : DIM(2) ], ...
                  numberOfPaddingVoxels + [ 1 : DIM(3) ], ...
                  : ) = priors;
    priors = paddedPriors;
    DIM = paddedDIM;
    transformMatrix = [ 1 0 0 -numberOfPaddingVoxels; 0 1 0 -numberOfPaddingVoxels; 0 0 1 -numberOfPaddingVoxels; 0 0 0 1 ] * transformMatrix;
  end
  
  % Convert to uint16 format
  priors = uint16( priors * ( 2^16 - 1 ) + .5 );
  
  
end



%
meshCollection = kvlCreateMeshCollection( priors, meshSize );
outputFileNameBase = fullfile( savePath, 'atlasForAffineRegistration' );
affineRegistrationMeshCollectionFileName = [ outputFileNameBase '_meshCollection.txt' ];
kvlWriteMeshCollection( meshCollection, affineRegistrationMeshCollectionFileName );   

tmp = size( priors );
DIM = tmp( 1 : 3 );
numberOfClasses = tmp( 4 );
templateImageBuffer = reshape( reshape( double( priors ) / (2^16-1), [ prod( DIM ) numberOfClasses ] ) * [ 1 : numberOfClasses ]', ...
                               DIM );
template = kvlCreateImage( single( templateImageBuffer ) );                               
affineRegistrationTemplateFileName = [ outputFileNameBase '_template.nii' ];
transform = kvlCreateTransform( double( transformMatrix ) );
kvlWriteImage( template, affineRegistrationTemplateFileName, transform );
