function prepareExVivoImages( subjectPath, flipAngles, backgroundThreshold, rotationAngleAroundZAxis, rotationAngleAroundYAxis, rotationAngleAroundXAxis, templateFileName, voxelIndicesInImage, voxelIndicesInTemplate, useSimilarityTransform )
%
%  
%

if ( nargin < 1 )
  % Test drive
  addpath /home/koen/software/freesurferGit/freesurfer/GEMS2-Release/bin/
  
  subjectPath = '/data/testing/exvivo/I48_Bay3_20170217/mri/parameter_maps/1mm/';

  flipAngles = [ 20 5 ];
  
  backgroundThreshold = 30;

  rotationAngleAroundZAxis = 0 * pi/2;
  rotationAngleAroundYAxis = 0 * pi/2;
  rotationAngleAroundXAxis = 0 * pi/2;

  templateFileName = '/data/testing/exvivo/Own_4classes_30x30x30_exvivo_template.nii';
  voxelIndicesInImage = [ 54 151 48; 84 97 48; 73 116 97; 97 39 59 ];
  voxelIndicesInTemplate = [ 85 86 147; 85 75 80; 146 87 119; 91 103 12 ];
  useSimilarityTransform = true;
  
  prepareExVivoImages( subjectPath, backgroundThreshold, rotationAngleAroundZAxis, rotationAngleAroundYAxis, rotationAngleAroundXAxis, voxelIndicesInImage, voxelIndicesInTemplate, useSimilarityTransform );
  
  return
  
end


if ( nargin < 7 )
  templateFileName = '';
  voxelIndicesInImage = [];
  voxelIndicesInTemplate = [];
  useSimilarityTransform = true;
end


% 
numberOfFlipAngles = length( flipAngles );
fileNames = cell( numberOfFlipAngles, 1 );
for flipAngleNumber = 1 : numberOfFlipAngles
  flipAngle = flipAngles( flipAngleNumber );
  fileNames{ flipAngleNumber } = fullfile( subjectPath, [ 'flash' num2str( flipAngle ) '.mgz' ] );
end

%
for flipAngleNumber = 1 : numberOfFlipAngles
  fileName = fileNames{ flipAngleNumber };
  [ image, imageTransform ] = kvlReadImage( fileName );
  imageBuffer = kvlGetImageBuffer( image );
  mask = ( imageBuffer > 0 ) .* ( imageBuffer < 1e5 ); % Generate a reasonable range
  if ( flipAngleNumber > 1 )
    % Generate a background mask
    mask = ( imageBuffer > backgroundThreshold );
    
    % Largest connected component
    cc = bwconncomp( mask );
    numPixels = cellfun( @numel, cc.PixelIdxList );
    [ biggest, idx ] = max( numPixels );
    tmp = zeros( size( mask ) );
    tmp( cc.PixelIdxList{ idx } ) = 1;
    mask = tmp;
    
    % se = strel( 'sphere', 25 );
    % mask = imerode( mask, se );
    % mask = imdilate( mask, se );
    
    figure
    showImage( mask );
  end
    
  maskedImageBuffer = imageBuffer .* mask;
  figure
  showImage( maskedImageBuffer )
  maskedImage = kvlCreateImage( maskedImageBuffer );
  newFileName = [ fileName( 1 : end-4 ) '_masked.mgz' ];
  kvlWriteImage( maskedImage, newFileName, imageTransform );

  % Remember file name to pass on to next processing steps
  fileNames{ flipAngleNumber } = newFileName;
  
end


% 
for flipAngleNumber = 1 : numberOfFlipAngles
  fileName = fileNames{ flipAngleNumber };       
  [ image, transform ] = kvlReadImage( fileName );
  transformMatrix = double( kvlGetTransformMatrix( transform ) );
  rotationMatrix = [  cos( rotationAngleAroundZAxis ) -sin( rotationAngleAroundZAxis ) 0 0; ...
                      sin( rotationAngleAroundZAxis ) cos( rotationAngleAroundZAxis ) 0 0; ...
                      0 0 1 0; ...
                      0 0 0 1 ];
  transformMatrix = rotationMatrix * transformMatrix;
  rotationMatrix = [ cos( rotationAngleAroundYAxis ) 0 -sin( rotationAngleAroundYAxis ) 0; ...
                     0 1 0 0; ...
                     sin( rotationAngleAroundYAxis ) 0 cos( rotationAngleAroundYAxis ) 0; ...
                     0 0 0 1 ];
  transformMatrix = rotationMatrix * transformMatrix;
  rotationMatrix = [ 1 0 0 0; ...
                     0 cos( rotationAngleAroundXAxis ) -sin( rotationAngleAroundXAxis ) 0; ...
                     0 sin( rotationAngleAroundXAxis ) cos( rotationAngleAroundXAxis ) 0; ...
                     0 0 0 1 ];
  transformMatrix = rotationMatrix * transformMatrix;
  newFileName = [ fileName( 1 : end-4 ) '_rotated.mgz' ];
  kvlWriteImage( image, newFileName, kvlCreateTransform( transformMatrix ) );
  
  % Remember file name to pass on to next processing steps
  fileNames{ flipAngleNumber } = newFileName;

end


if ~isempty( voxelIndicesInImage )

  % Compute the image-to-image tranform (from image to template)
  if ~useSimilarityTransform
    % Fully affine (tricky - need many points to get a reasonable transformation)
    numberOfVoxelsSelected = size( voxelIndicesInTemplate, 1 );
    X = [ voxelIndicesInImage ones( numberOfVoxelsSelected, 1 ) ];
    Y = [ voxelIndicesInTemplate ones( numberOfVoxelsSelected, 1 ) ];
    imageToImageTransformMatrix = ( ( X' * X ) \ ( X' * Y ) )'
  else
    % Rigid or similarity (rigid + isotropic scaling) transformation
    numberOfVoxelsSelected = size( voxelIndicesInTemplate, 1 );
    
    [ template, templateTransform ] = kvlReadImage( templateFileName );
    templateTransformMatrix = double( kvlGetTransformMatrix( templateTransform ) );

    [ image, imageTransform ] = kvlReadImage( fileNames{ 1 } );
    imageTransformMatrix = double( kvlGetTransformMatrix( imageTransform ) );

    X = imageTransformMatrix * [ voxelIndicesInImage ones( numberOfVoxelsSelected, 1 ) ]'; X = X( 1 : 3, : )';
    Y = templateTransformMatrix * [ voxelIndicesInTemplate ones( numberOfVoxelsSelected, 1 ) ]'; Y = Y( 1 : 3, : )';
    
    meanX = sum( X ) / numberOfVoxelsSelected;
    meanY = sum( Y ) / numberOfVoxelsSelected;
    Xcentered = X - repmat( meanX, [ numberOfVoxelsSelected 1 ] );
    Ycentered = Y - repmat( meanY, [ numberOfVoxelsSelected 1 ] );
    [ U, S, V ] = svd( Xcentered' * Ycentered );
    rotationMatrix = V * diag( [ 1 1 det( V * U ) ] ) * U';
    useRigid = false;
    if useRigid
      isotropicScaling = 1;
    else
      tmp = rotationMatrix' * Ycentered';
      numerator = 0;
      denominator = 0;
      for selectedVoxelNumber = 1 : numberOfVoxelsSelected
        numerator = numerator + Xcentered( selectedVoxelNumber, : ) * tmp( :, selectedVoxelNumber );
        denominator = denominator + Xcentered( selectedVoxelNumber, : ) * Xcentered( selectedVoxelNumber, : )';
      end
      isotropicScaling = numerator / denominator;
    end  
    translationVector = meanY' - isotropicScaling * rotationMatrix * meanX';
    worldToWorldTransformMatrix = [ isotropicScaling * rotationMatrix translationVector; 0 0 0 1 ];

    % Double-check that distance in world space is indeed reduced
    Ypredicted = worldToWorldTransformMatrix * [ X ones( numberOfVoxelsSelected, 1 ) ]'; Ypredicted = Ypredicted( 1:3,:)'
    sum( ( Y(:) - X(:) ).^2 )
    sum( ( Y(:) - Ypredicted(:) ).^2 )

    imageToImageTransformMatrix = inv( templateTransformMatrix ) * worldToWorldTransformMatrix * imageTransformMatrix;
      
  end
    
    
  % 
  [ template, templateTransform ] = kvlReadImage( templateFileName );
  templateTransformMatrix = double( kvlGetTransformMatrix( templateTransform ) );
  desiredImageToWorldTransformMatrix = templateTransformMatrix * imageToImageTransformMatrix
  for flipAngleNumber = 1 : numberOfFlipAngles
    imageFileName = fileNames{ flipAngleNumber };
    [ image, imageTransform ] = kvlReadImage( imageFileName );
    currentImageToWorldTransformMatrix = kvlGetTransformMatrix( imageTransform );
    desiredWorlToWorldTransformation = desiredImageToWorldTransformMatrix * inv( currentImageToWorldTransformMatrix )
    kvlWriteImage( image, imageFileName, ...
                  kvlCreateTransform( desiredImageToWorldTransformMatrix ) );
  end


end
