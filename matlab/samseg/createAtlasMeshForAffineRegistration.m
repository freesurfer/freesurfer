%  
%  Create a regular (low-res) tetrahedral mesh atlas from existing mesh atlas and/or volumetric prior
%  The result is a mesh with the following order of classes: background, WM, GM, CSF, [bone, soft tissue]
%  (last two are optional, depending on input mesh/atlas that is given)
%

%
close all
clear all

% Options go here 
meshSize = [ 40 40 40 ];
meshType = 3;  % 1 = SPM12 with 6 classes; 2 = SPM12 with 4 classes; 3 = own with 4 classes 


%
addpath /home/koen/software/freesurferGit/freesurfer/GEMS2-Release/bin/


%
% Part I: get the existing prior in volumetric form (uint16)
%
if ( ( meshType == 1 ) || ( meshType == 2 ) )
  % Dealing with SPM12. 
  
  % Read in the correct TPMs extracted manually from the 6-channel NIFTI atlas. 
  DIM = [ 121 145 121 ];
  if ( meshType == 1 )
    numberOfClasses = 6;
  else
    numberOfClasses = 4;
  end
  
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
    transformMatrix = kvlGetTransformMatrix( transform ); 
    transformMatrix = [ 1 0 0 -numberOfPaddingVoxels; 0 1 0 -numberOfPaddingVoxels; 0 0 1 -numberOfPaddingVoxels; 0 0 0 1 ] * transformMatrix;
    transform = kvlCreateTransform( double( transformMatrix ) );
  end
  
  % Convert to uint16 format
  priors = uint16( priors * ( 2^16 - 1 ) + .5 );
  
  
else

  % Our own atlas
  templateFileName = '/data/testing/atlas/mni305_masked_autoCropped.mgz';
  meshCollectionFileName = '/data/testing/atlas/CurrentMeshCollection30New.txt.gz';
  compressionLookupTableFileName = '/data/testing/atlas/namedCompressionLookupTable.txt';

  meshCollection = kvlReadMeshCollection( meshCollectionFileName );
  mesh = kvlGetMesh( meshCollection, -1 );

  % Reduce number of classes to just 
  originalAlphas = kvlGetAlphasInMeshNodes( mesh );
  [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
  sameGaussianParameters = cell(0,0);
  sameGaussianParameters{1} = [ 0 85 ];  % Background is a separate class
  sameGaussianParameters{2} = [ 2 7 16 41 46 28 60 ]; % Force all white matter structures to have the same intensity model
  sameGaussianParameters{3} = [ 3 8 42 47 11 50 17 53 18 54 26 58 77 80 10 49 12 51 13 52 ]; % Same for several gray matter structures
  sameGaussianParameters{4} = [ 4 5 14 15 24 43 44 72 30 62 31 63 ]; % Same for CSF
  [ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
  size( reducedAlphas )
  max( abs( sum( reducedAlphas, 2 ) - 1 ) )  % Make sure these vectors really sum to 1
  kvlSetAlphasInMeshNodes( mesh, reducedAlphas )

  
  [ template, transform ] = kvlReadImage( templateFileName ); 
  templateImageBuffer = kvlGetImageBuffer( template );

  priors = kvlRasterizeAtlasMesh( mesh, size( templateImageBuffer ) );

end



%
% Part II: convert the volumetric priors into a mesh, and write out (incl. volumetric template)
%
meshCollection = kvlCreateMeshCollection( priors, meshSize );
meshTypeNames = str2mat( 'SPM12_6classes', 'SPM12_4classes', 'Own_4classes' );
outputFileNameBase = [ deblank( meshTypeNames( meshType, : ) ) '_' ...
                           num2str( meshSize( 1 ) ) 'x' num2str( meshSize( 2 ) ) 'x' num2str( meshSize( 3 ) ) ];
meshCollectionFileName = [ outputFileNameBase '_meshCollection.txt' ];
kvlWriteMeshCollection( meshCollection, meshCollectionFileName );   

tmp = size( priors );
DIM = tmp( 1 : 3 );
numberOfClasses = tmp( 4 );
templateImageBuffer = reshape( reshape( double( priors ) / (2^16-1), [ prod( DIM ) numberOfClasses ] ) * [ 1 : numberOfClasses ]', ...
                               DIM );
template = kvlCreateImage( single( templateImageBuffer ) );                               
templateFileName = [ outputFileNameBase '_template.nii' ];
kvlWriteImage( template, templateFileName, transform );




