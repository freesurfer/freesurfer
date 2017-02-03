function transformedTemplateFileName = samseg_coregisterAtlas_SPM( imageFileName, templateFileName, meshCollectionFileName, compressionLookupTableFileName, initializeUsingCenterOfGravityAlignment, savePath )
%
% transformedTemplateFileName = coregisterAtlas_SPM( imageFileName, templateFileName, meshCollectionFileName, compressionLookupTableFileName )
%
% Affinely coregister a mesh-based whole-brain atlas to an MR image using SPM8's spm_maff function.
% This is equivalent to Oula's coregisterAtlas_SPM.m
%
% $Id: samseg_coregisterAtlas_SPM.m,v 1.1 2017/01/26 00:19:13 greve Exp $

showFigures = false;
useRegularization = true;

if nargin < 6
  savePath = '';
end

if nargin < 5
  initializeUsingCenterOfGravityAlignment = false;
end



% Read image
[ image, transform ] = kvlReadImage( imageFileName );
imageBuffer = kvlGetImageBuffer( image );
transformMatrix = kvlGetTransformMatrix( transform )
if showFigures
  figure
  showImage( imageBuffer )
end


% Write out again in NIFTI format
niftiImageFileName = [ tempname '.nii' ]; 
kvlWriteImage( image, niftiImageFileName )
%  [ image2, transform2 ] = kvlReadImage( niftiImageFileName );



% Compute the priors for WM, GM, and CSF
[ template, templateTransform ] = kvlReadImage( templateFileName );
templateTransformMatrix = kvlGetTransformMatrix( templateTransform )

% Read atlas mesh
meshCollection = kvlReadMeshCollection( meshCollectionFileName, kvlCreateTransform( eye(4) ), 0.1 );
mesh = kvlGetMesh( meshCollection, -1 );


% Compute the "reduced" alphas - those referring to the "super"-structures 
originalAlphas = kvlGetAlphasInMeshNodes( mesh );
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
sameGaussianParameters = cell(0,0);
sameGaussianParameters{1} = [ 0 ];  % Background is a separate class
sameGaussianParameters{2} = [ 2 7 16 41 46 28 60 85 ]; % Force all white matter structures to have the same intensity model
sameGaussianParameters{3} = [ 3 8 42 47 11 50 17 53 18 54 26 58 77 80 10 49 12 51 13 52 ]; % Same for several gray matter structures
sameGaussianParameters{4} = [ 4 5 14 15 24 43 44 72 30 62 31 63 ]; % Same for CSF
[ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
size( reducedAlphas )
max( abs( sum( reducedAlphas, 2 ) - 1 ) )  % Make sure these vectors really sum to 1


% Set the reduced alphas to be the alphas of the mesh
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )

if 0
  % Smooth mesh somewhat
  meshSmoothingSigma = 2;
  fprintf( 'Smoothing mesh with kernel size %f ...', meshSmoothingSigma )
  kvlSmoothMesh( mesh, meshSmoothingSigma )
  fprintf( 'done\n' )
end


priors = kvlRasterizeAtlasMesh( mesh, size( kvlGetImageBuffer( template ) ) );
if showFigures
  figure
  for reducedLabel = 1 : size( reducedAlphas, 2 )
    subplot( 2, 2, reducedLabel )
    showImage( priors( :, :, :, reducedLabel ) )
  end
end

% b0 should be 4x1 struct with double-precision priors for GM, WM, CSF, and other (<- not used)
b0 = cell( 4, 1 );
b0{1} = double( priors( :, :, :, 3 ) ) / ( 2^16 - 1 );
b0{2} = double( priors( :, :, :, 2 ) ) / ( 2^16 - 1 );
b0{3} = double( priors( :, :, :, 4 ) ) / ( 2^16 - 1 );
b0{4} = double( priors( :, :, :, 1 ) ) / ( 2^16 - 1 );
if 0
  % kvlWriteImage( kvlCreateImage( single( b0{2} ) ), 'debugWhite.nii', templateTransform );
  kvlWriteImage( kvlCreateImage( single( b0{1} + b0{2} + b0{3} ) ), 'debugBrain.nii', templateTransform );
end


% Prepare image and atlas according to how it's done in kvl_preproc.m
V = spm_vol( niftiImageFileName );
samp = 3;
d         = V(1).dim(1:3);
vx        = sqrt(sum(V(1).mat(1:3,1:3).^2));
sk        = max([1 1 1],round(samp*[1 1 1]./vx));
[x0,y0,o] = ndgrid(1:sk(1):d(1),1:sk(2):d(2),1);
z0        = 1:sk(3):d(3);

% Careful now - the image-to-world transform in SPM is slightly off the plain image-to-world transform.
% Two differences are happening: 
%   1) SPM uses 1-based indexing (e.g., first voxel has index (1,1,1) instead of (0,0,0)
%   2) The NIFTI reader in SPM somehow flips x and y axes
% Therefore, to compute the ITK image-to-world transform from the SPM one, we have that
% 
% M =  [ -1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1 ]  * M_spm * [eye(4,3) [ 1 1 1 1]' ] 
%
% To double-check that does is indeed the case for a non-trivial (2mm voxel size) dataset, do
%
%      b = spm_get_space( '/home/koen/software/spm8/tpm/grey.nii' )
%      [ tmp, tmpT ] = kvlReadImage( '/home/koen/software/spm8/tpm/grey.nii' );
%      a = kvlGetTransformMatrix( tmpT )
%      [ -1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1 ]  * b * [eye(4,3) [ 1 1 1 1]' ]
%


% Prepare as interface to spm_maff.m
P = V;
x = {x0,y0,z0};
MF = [ -1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1 ] ...
     * double( templateTransformMatrix ) ...
     * [eye(4,3) [ -1 -1 -1 1]' ]; % Image-to-world of priors (in SPM convention)
M = eye( 4 );  % Starting estimate
if ~useRegularization
  regtype = 'none';
  ff = 0;
else
  regtype = 'mni';
  ff = 5;
  ff     = max(1,ff^3/prod(sk)/abs(det(V.mat(1:3,1:3))));
end


if initializeUsingCenterOfGravityAlignment
  disp( 'Initializing by aligning centers of gravity in the images' )

  % Start by matching the center of gravity of the two images
  [ xtmp, ytmp, ztmp ] = ndgrid( 1 : size( imageBuffer, 1 ), 1 : size( imageBuffer, 2 ), 1 : size( imageBuffer, 3 ) );
  centerOfGravityImage = [ xtmp(:) ytmp(:) ztmp(:) ]' * imageBuffer(:) / sum( imageBuffer(:) ); % In image coordinates
  if 0
    figure
    tmp = imageBuffer;
    tmp( round( centerOfGravityImage(1) ), round( centerOfGravityImage(2) ), round( centerOfGravityImage(3) ) ) = ...
      3 * max( imageBuffer(:) );
    showImage( tmp, [ round( centerOfGravityImage(1) ) round( centerOfGravityImage(2) ) round( centerOfGravityImage(3) ) ] )
  end
  centerOfGravityImage = spm_get_space( niftiImageFileName ) * [ centerOfGravityImage; 1 ] % In world coordinates


  tmp = b0{1} + b0{2} + b0{3};
  [ xtmp, ytmp, ztmp ] = ndgrid( 1 : size( tmp, 1 ), 1 : size( tmp, 2 ), 1 : size( tmp, 3 ) );
  centerOfGravityAtlas = [ xtmp(:) ytmp(:) ztmp(:) ]' * tmp(:) / sum( tmp(:) ); % In image coordinates
  centerOfGravityAtlas = MF * [ centerOfGravityAtlas; 1 ] % In world coordinates

  % In order to map to centers of gravity initially, we have to initialize M so that
  %   centerOfGravityAtlas = M * centerOfGravityImage
  M = [ eye( 3 ) centerOfGravityAtlas( 1 : 3) - centerOfGravityImage( 1 : 3 ); [ 0 0 0 1 ] ]
end



% Do the registration
if 1
  disp( 'Starting SPM affine registration...' )
  if useRegularization
    M_spm = spm_maff( P, x, b0, MF, M, regtype, ff*100 );
    M_spm = spm_maff( P, x, b0, MF, M_spm, regtype, ff );
  else
    M_spm = spm_maff( P, x, b0, MF, M, regtype, ff );
  end
  disp( '...done' )
else
  M_spm = M;
end



% Get the result out
%  imageToTemplate_imageToImageTransform = MF \ M_spm * spm_get_space( niftiImageFileName );
%  templateToImage_imageToImageTransform = inv( imageToTemplate_imageToImageTransform );
%  kvlTransformMeshCollection( meshCollection, kvlCreateTransform( single( templateToImage_imageToImageTransform ) ) );
%  kvlWriteMeshCollection( meshCollection, 'affinelyCoregisteredMeshCollection.txt' );
%  kvlWriteImage( template, 'coregisteredTemplate.nii' );
% Set the image-to-world transform X so that 
%   X \ spm_get_space( niftiImageFileName ) = imageToTemplate_imageToImageTransform
% where
%   imageToTemplate_imageToImageTransform = MF \ M_spm * spm_get_space( niftiImageFileName );
% which means
%   inv(MF) * M_spm = inv( X )
% and therefore
%   X = M_spm \ MF
%
desiredImageToWorldTransform =  [ -1 0 0 0; 0 -1 0 0; 0 0 1 0; 0 0 0 1 ] ...
                                 * ( M_spm \ MF ) ...
                                 * [eye(4,3) [ 1 1 1 1]' ] 

[ pathstr, name, ext ] = fileparts( templateFileName );
transformedTemplateFileName = fullfile( savePath, [ name '_coregistered' ext ]);
kvlWriteImage( template, transformedTemplateFileName, kvlCreateTransform( desiredImageToWorldTransform ) );
%  [ dummy, dummyT ] = kvlReadImage( transformedTemplateFileName );
%  T = kvlGetTransformMatrix( dummyT )


% Clean up
delete( niftiImageFileName );
