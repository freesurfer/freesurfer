% This is a version of Oula's registerAtlas.m. The inputs (eg,
% imageFileName) are assumed to have been set prior to calling this
% function (eg, in (run_samseg.m) This is a frontend for
% samseg_coregisterAtlas_SPM.m. The only real purpose it serves is to
% show images if showFigs != 0
% 
% $Id: samseg_registerToAtlas.m,v 1.1 2017/01/26 00:19:20 greve Exp $
%

% 
%imageFileName = './exampleData/bert/001.mgz';
showFigs = 0;

% 
%templateFileName = 'mni305_masked_autoCropped.mgz';
%meshCollectionFileName = 'CurrentMeshCollection30New.txt.gz';
%compressionLookupTableFileName = 'namedCompressionLookupTable.txt';


kvlClear; % Clear all the wrapped C++ stuff
close all;

if(showFigs)
    % Show situation before doing any affine registration
    [ image, transform ] = kvlReadCroppedImage( imageFileName, templateFileName );
    imageBuffer = kvlGetImageBuffer( image );
    figure
    showImage( imageBuffer );
    K = 0.1;
    meshCollection = kvlReadMeshCollection( meshCollectionFileName, transform, K );
    mesh = kvlGetMesh( meshCollection, -1 );

    labelNumber = 0;
    backgroundPrior = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ), labelNumber );
    smoothingSigma = 2; % sqrt of the variance of a Gaussian blurring kernel
    smoothedBackgroundPrior = kvlSmoothImageBuffer( backgroundPrior, smoothingSigma );
    brainMask = ( 1 - single( smoothedBackgroundPrior ) / 65535 ) > .01;
    oldImageBuffer = imageBuffer;
    imageBuffer( find( ~brainMask ) ) = 0;
    figure
    subplot( 2, 2, 1 )

    showImage( oldImageBuffer );
    subplot( 2, 2, 2 )
    showImage( imageBuffer )
    subplot( 2, 2, 3 )
    showImage( oldImageBuffer .* single( ~brainMask ) );
end


% Do the job
[pathstr,name,ext] = fileparts( imageFileName );

% Add rician noise
imageHandle = kvlReadImage(imageFileName );
imageBuffer = kvlGetImageBuffer( imageHandle );
RicianNoise = ( ( randn( size( imageBuffer ) ) ).^2 + ( randn( size( imageBuffer ) ) ).^2 ).^(1/2);
perfectlyZeroMask = ( imageBuffer == 0 ) ;
imageBuffer( find( perfectlyZeroMask(:) ) ) = 3 * RicianNoise( find(perfectlyZeroMask(:) ) );
kvlSetImageBuffer( imageHandle, imageBuffer );
imageFileNoise = sprintf('%s/%s_noise.mgz',savePath,name);
kvlWriteImage( imageHandle, imageFileNoise );

transformedTemplateFileName = samseg_coregisterAtlas_SPM(imageFileNoise,templateFileName,meshCollectionFileName,compressionLookupTableFileName,initializeUsingCenterOfGravityAlignment,savePath);

 
% Show after registration
if(showFigs)
  [ image, transform ] = kvlReadCroppedImage( imageFileName, transformedTemplateFileName );
  imageBuffer = kvlGetImageBuffer( image );
  figure
  showImage( imageBuffer ); 
  K = 0.1;
  meshCollection = kvlReadMeshCollection( meshCollectionFileName, transform, K );
  mesh = kvlGetMesh( meshCollection, -1 );

  labelNumber = 0;
  backgroundPrior = kvlRasterizeAtlasMesh( mesh, size( imageBuffer ), labelNumber );
  figure
  subplot( 2, 2, 1 )
  showImage( backgroundPrior )
  subplot( 2, 2, 2 )
  showImage( mosaicImages( 2^16 - 1 - double( backgroundPrior ), double( imageBuffer ), 10 ) )
  smoothingSigma = 2; % sqrt of the variance of a Gaussian blurring kernel
  smoothedBackgroundPrior = kvlSmoothImageBuffer( backgroundPrior, smoothingSigma );
  subplot( 2, 2, 3 )
  showImage( smoothedBackgroundPrior )
  brainMask = ( 1 - single( smoothedBackgroundPrior ) / 65535 ) > .01;
  oldImageBuffer = imageBuffer;
  imageBuffer( find( ~brainMask ) ) = 0;
  subplot( 2, 2, 4 )
  showImage( imageBuffer )

  figure
  subplot( 2, 2, 1 )
  showImage( oldImageBuffer );
  subplot( 2, 2, 2 )
  showImage( imageBuffer )
  subplot( 2, 2, 3 )
  showImage( oldImageBuffer .* single( ~brainMask ) );
end