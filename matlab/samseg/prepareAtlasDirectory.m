function prepareAtlasDirectory( directoryName, compressionLookupTableFileName, meshCollectionFileName, templateFileName, ...
                                FreeSurferLookupTableFileName, sharedGMMParameters, ...
                                missingStructuresMergeOptions, affineAtlasAdditionalMergeOptions, ...
                                smoothingSigmaForLevel1, meshCollectionFileNameForLevel1, ...
                                smoothingSigmaForAffine, meshCollectionFileNameForAffine )
%
%
% 

if ( nargin == 0 )
  % Test-drive ourselves

  % Point to various input and output files/directories
  directoryName = '/tmp/myAtlas';
  compressionLookupTableFileName = ...
    '/data/testing/atlas/Buckner39AtlasWithMoreClassesAndEyeballs/atlases/10SubjectAtlas3X/result/compressionLookupTable.txt';
  meshCollectionFileName = ... 
    '/data/testing/atlas/Buckner39AtlasWithMoreClassesAndEyeballs/atlases/10SubjectAtlas3X/result/CurrentMeshCollection30.gz';
  templateFileName = '/data/testing/atlas/Buckner39AtlasWithMoreClassesAndEyeballs/mni305_1mm_emptyBox.nii';
  FreeSurferLookupTableFileName = '/data/testing/atlas/Buckner39AtlasWithMoreClassesAndEyeballs/FreeSurferColorLUT.txt';


  % Specify which classes share the same intensity Gaussian mixture model (and the number of components within each model)
  sharedGMMParameters = struct;
  sharedGMMParameters( 1 ).mergedName = 'Unknown'; % Background and bone
  sharedGMMParameters( 1 ).searchStrings = { 'Unknown', 'Skull' };
  sharedGMMParameters( 1 ).numberOfComponents = 1;
  sharedGMMParameters( 2 ).mergedName = 'GlobalWM'; % WM
  sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
  sharedGMMParameters( 2 ).numberOfComponents = 2;
  sharedGMMParameters( 3 ).mergedName = 'GlobalGM'; % GM
  sharedGMMParameters( 3 ).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus', 'Amygdala', 'Accumbens', 'hypointensities' };
  sharedGMMParameters( 3 ).numberOfComponents = 3;
  sharedGMMParameters( 4 ).mergedName = 'GlobalCSF'; % CSF
  sharedGMMParameters( 4 ).searchStrings = { 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus', 'Fluid' };
  sharedGMMParameters( 4 ).numberOfComponents = 3;
  sharedGMMParameters( 5 ).mergedName = 'Thalamus'; % Thalamus
  sharedGMMParameters( 5 ).searchStrings = { 'Thalamus' };
  sharedGMMParameters( 5 ).numberOfComponents = 2;
  sharedGMMParameters( 6 ).mergedName = 'Pallidum'; % Pallidum
  sharedGMMParameters( 6 ).searchStrings = { 'Pallidum' };
  sharedGMMParameters( 6 ).numberOfComponents = 2;
  sharedGMMParameters( 7 ).mergedName = 'Putamen'; % Putamen
  sharedGMMParameters( 7 ).searchStrings = { 'Putamen' };
  sharedGMMParameters( 7 ).numberOfComponents = 2;
  sharedGMMParameters( 8 ).mergedName = 'Soft'; % Soft non-brain tissue
  sharedGMMParameters( 8 ).searchStrings = { 'Soft' };
  sharedGMMParameters( 8 ).numberOfComponents = 3;

  % If some structures are missing from your images (e.g., in ex vivo scans) you can specify those here.
  % They will then be removed from the atlases (by effectively merging them into the 'Unknown' class) and therefore
  % no attempt will be made to segment them.
  missingStructuresMergeOptions = struct;
  missingStructuresMergeOptions( 1 ).mergedName = 'Unknown';
  missingStructuresMergeOptions( 1 ).searchStrings = { 'Skull', 'Soft' };


  % For the atlas used in simple affine registration, a subset of the super-structures defined above (through sharedGMMParameters) 
  % is usually sufficient. Additional merge options for this purpose (e.g., moving some deep gray matter structures into global GM) 
  % can be specified here.
  affineAtlasAdditionalMergeOptions = struct;
  affineAtlasAdditionalMergeOptions( 1 ).mergedName = 'GlobalGM';
  affineAtlasAdditionalMergeOptions( 1 ).searchStrings = { 'Thalamus', 'Pallidum', 'Putamen' };


  % Specify how the smoother atlas used at the first multi-resolution level is obtained: by artifically smoothing the final
  % mesh (i.e., the one used at the second and final multi-resolution level) and/or by using a precomputed mesh with a simpler topology
  smoothingSigmaForLevel1 = 2.0; % In voxels
  if 0
    meshCollectionFileNameForLevel1 = [];
  else  
    meshCollectionFileNameForLevel1 = '/data/testing/atlas/Buckner39AtlasWithMoreClassesAndEyeballs/atlases/10SubjectAtlasMultiresolution/scratch/koenLogDir10SubjectAtlasMultires/CurrentMeshCollection3.gz'; % 1 is a factor 1.4, 3 is factor 2, 30 is factor 3  
  end

  % Specify how a simpler mesh for affine registration purposes is obtained (estimating 12 degrees of freedom really doesn't require
  % a supe-high-quality detailed tetrahedral mesh): either by re-meshing to a regular low-resolution mesh, or by using a precomputed 
  % mesh with a simpler topology
  smoothingSigmaForAffine = 2.0; % In voxels
  if 0
    meshCollectionFileNameForAffine = [];
  else  
    meshCollectionFileNameForAffine = '/data/testing/atlas/Buckner39AtlasWithMoreClassesAndEyeballs/atlases/10SubjectAtlasMultiresolution/scratch/koenLogDir10SubjectAtlasMultires/CurrentMeshCollection30.gz'; % 1 is a factor 1.4, 3 is factor 2, 30 is factor 3
  end
  
  
  % Let the beast go
  prepareAtlasDirectory( directoryName, compressionLookupTableFileName, meshCollectionFileName, templateFileName, ...
                         FreeSurferLookupTableFileName, sharedGMMParameters, ...
                         missingStructuresMergeOptions, affineAtlasAdditionalMergeOptions, ...
                         smoothingSigmaForLevel1, meshCollectionFileNameForLevel1, ...
                         smoothingSigmaForAffine, meshCollectionFileNameForAffine );
  
  return

end  % End test if test-drive ourselves




% Try to create directory
[ success, message ] = mkdir( directoryName );
if ~success
  error( [ 'Can''t create directory: ' directoryName ' (' message ')' ] )
end

% Copy the FreeSurferLookupTable
[ success, message ] = copyfile( FreeSurferLookupTableFileName, fullfile( directoryName, 'modifiedFreeSurferColorLUT.txt' ) );
if ~success
  error( [ 'Can''t copy file : ' FreeSurferLookupTableFileName ' to directory ' directoryName ' (' message ')' ] )
end

% Write sharedGMMParameters to file
kvlWriteSharedGMMParameters( sharedGMMParameters, fullfile( directoryName, 'sharedGMMParameters.txt' ) )


% Read mesh 
meshCollection = kvlReadMeshCollection( meshCollectionFileName );
mesh = kvlGetMesh( meshCollection, -1 );
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
alphas = kvlGetAlphasInMeshNodes( mesh );

% Only retain the reference (average across subjects) position -- the rest is nobody's bussiness and having to read
% dozens of warped mesh positions to the training subjects for no purpose at all is slow
referencePosition = kvlGetMeshNodePositions( mesh );
kvlSetMeshCollectionPositions( meshCollection, referencePosition, referencePosition );


% Remove non-existing structures
[ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, missingStructuresMergeOptions, FreeSurferLabels, colors );
kvlSetAlphasInMeshNodes( mesh, alphas );


% Write the ultimate atlas (at multi-resolution level 2)
kvlWriteMeshCollection( meshCollection, fullfile( directoryName, 'atlas_level2.txt' ) );


% Write the blurrier atlas (at multi-resolution level 1)
if ~isempty( meshCollectionFileNameForLevel1 )
  % Use a user-specified mesh with a simpler topology - follow the exact same pipeline 
  meshCollection = kvlReadMeshCollection( meshCollectionFileNameForLevel1 );
  downMesh = kvlGetMesh( meshCollection, -1 );

  [ downFreeSurferLabels, downNames, downColors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
  downAlphas = kvlGetAlphasInMeshNodes( downMesh );
  
  downReferencePosition = kvlGetMeshNodePositions( downMesh );
  kvlSetMeshCollectionPositions( meshCollection, downReferencePosition, downReferencePosition );
  
  [ downAlphas, downNames, downFreeSurferLabels, downColors ] = ...
       kvlMergeAlphas( downAlphas, downNames, missingStructuresMergeOptions, downFreeSurferLabels, ...
                       downColors );
  kvlSetAlphasInMeshNodes( downMesh, downAlphas );
end
% Use spatial smoothing
smoothedMeshCollection = kvlSmoothMeshCollection( meshCollection, smoothingSigmaForLevel1 );
kvlWriteMeshCollection( smoothedMeshCollection, fullfile( directoryName, 'atlas_level1.txt' ) );


  
% Write the corresponding compressionLookupTable
kvlWriteCompressionLookupTable( fullfile( directoryName, 'compressionLookupTable.txt' ), FreeSurferLabels, names, colors );


% For the affine atlas, simplify the mesh to only contain global tissue types specified by the sharedGMMParameters, and 
% additionally merge whatever the user thinks can be merged for affine registration purposes
[ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, sharedGMMParameters, FreeSurferLabels, colors );
[ alphas, names, FreeSurferLabels, colors ] = kvlMergeAlphas( alphas, names, affineAtlasAdditionalMergeOptions, FreeSurferLabels, colors );
kvlSetAlphasInMeshNodes( mesh, alphas );


% Create a custom template that color-codes the affine atlas label probabilities (just for visualization purposes). 
[ template, transform ] = kvlReadImage( templateFileName ); 
templateImageBuffer = kvlGetImageBuffer( template );
priors = kvlRasterizeAtlasMesh( mesh, size( templateImageBuffer ) );
tmp = size( priors );
DIM = tmp( 1 : 3 );
numberOfClasses = tmp( 4 );
templateImageBuffer = reshape( reshape( double( priors ) / (2^16-1), [ prod( DIM ) numberOfClasses ] ) * [ 1 : numberOfClasses ]', ...
                               DIM );
template = kvlCreateImage( single( templateImageBuffer ) );                               
kvlWriteImage( template, fullfile( directoryName, 'template.nii' ), transform );


if isempty( meshCollectionFileNameForAffine )
  % Smooth affine atlas
  kvlSmoothMesh( mesh, smoothingSigmaForAffine );
  priors = kvlRasterizeAtlasMesh( mesh, size( templateImageBuffer ) );

  % Re-mesh the affine atlas
  affineMeshCollection = kvlCreateMeshCollection( priors,  [ 30 30 30 ] );
else

  % Read mesh 
  affineMeshCollection = kvlReadMeshCollection( meshCollectionFileNameForAffine );
  affineMesh = kvlGetMesh( affineMeshCollection, -1 );
  [ affineFreeSurferLabels, affineNames, affineColors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
  affineAlphas = kvlGetAlphasInMeshNodes( affineMesh );

  % Only retain the reference (average across subjects) position -- the rest is nobody's bussiness and having to read
  % dozens of warped mesh positions to the training subjects for no purpose at all is slow
  affineReferencePosition = kvlGetMeshNodePositions( affineMesh );
  kvlSetMeshCollectionPositions( affineMeshCollection, affineReferencePosition, affineReferencePosition );

  % Remove non-existing structures
  [ affineAlphas, affineNames, affineFreeSurferLabels, affineColors ] = ...
        kvlMergeAlphas( affineAlphas, affineNames, missingStructuresMergeOptions, affineFreeSurferLabels, affineColors );

  % For the affine atlas, simplify the mesh to only contain global tissue types specified by the sharedGMMParameters, and 
  % additionally merge whatever the user thinks can be merged for affine registration purposes
  [ affineAlphas, affineNames, affineFreeSurferLabels, affineColors ] = ...
        kvlMergeAlphas( affineAlphas, affineNames, sharedGMMParameters, affineFreeSurferLabels, affineColors );
  [ affineAlphas, affineNames, affineFreeSurferLabels, affineColors ] = ...
        kvlMergeAlphas( affineAlphas, affineNames, affineAtlasAdditionalMergeOptions, affineFreeSurferLabels, affineColors );
  kvlSetAlphasInMeshNodes( affineMesh, affineAlphas );

  % Smooth affine atlas
  affineMeshCollection = kvlSmoothMeshCollection( affineMeshCollection, smoothingSigmaForAffine );
    
end
    
    
% Write out the affine atlas mesh
kvlWriteMeshCollection( affineMeshCollection, fullfile( directoryName, 'atlasForAffineRegistration.txt' ) );



