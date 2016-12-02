%
close all
clear all

%
addpath /home/koen/software/mergedCodeEugenioAndOula/richardGEMS-Release/bin

load /data/tmp.mat

%kvlSetMaximumNumberOfThreads( 1 );


for nima = 1:numberOfImages
    [ images(nima), transform ] = kvlReadCroppedImage( imageFileNames{nima}, transformedTemplateFileName ); %Get the pointers to image and the corresponding transform
    imageBuffers(:,:,:,nima) = kvlGetImageBuffer( images(nima) ); %Get the actual imageBuffer
end


meshCollection = kvlReadMeshCollection( meshCollectionFileName, transform, K );
mesh = kvlGetMesh( meshCollection, -1 );

for nima = 1 : numberOfImages
    imageBuffer = kvlGetImageBuffer( images(nima) );
    imageBuffer( find( ~brainMask ) ) = 0;
    kvlSetImageBuffer( images(nima), imageBuffer );
    imageBuffers(:,:,:,nima) = imageBuffer;
    clear imageBuffer;
end

%  kvlSetAlphasInMeshNodes( mesh, reducedAlphas );
kvlSetAlphasInMeshNodes(mesh, optimizerAlphas);

for nima = 1:numberOfImages
    buffer = kvlGetImageBuffer(images(nima));
    DIM = size(buffer);
    mask = zeros( DIM );
    mask( maskIndices ) = 1;
    buffer( maskIndices ) = 1000* log( buffer( maskIndices ) ); % The 1000 factor isn't really necessary; it's just
                                                                % easier to have Gaussian means of 100 350 than
                                                                % 0.1 and 0.35 for my human brain
    
    buffer = buffer .* mask;
    kvlSetImageBuffer( images(nima), buffer );
    imageBuffers(:,:,:,nima) = buffer;
    clear buffer
end

for nima = 1:numberOfImages
    buffer = kvlGetImageBuffer(images(nima));
    biasCorrectedImages(nima) = kvlCreateImage( buffer );
end



optimizer = kvlGetConjugateGradientOptimizer( mesh, biasCorrectedImages, transform );
kvlSetOptimizerProperties( optimizer, means, precisions );
[ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStep( optimizer );

optimizerGPU = kvlGetConjugateGradientOptimizerGPU( mesh, biasCorrectedImages, transform );
kvlSetOptimizerPropertiesGPU( optimizerGPU, means, precisions );
[ minLogLikelihoodTimesPriorGPU, maximalDeformationGPU ] = kvlDeformOneStepGPU( optimizerGPU );


if 0

  nodePositions = kvlGetMeshNodePositions( mesh );
  alphas = kvlGetAlphasInMeshNodes( mesh );

  kvlWriteMeshCollection( meshCollection, 'test.txt' )

  kvlWriteImage( biasCorrectedImages(1), 'test.nii' );
  
end


return



%
historyOfMinLogLikelihoodTimesPrior = [];
historyOfMinLogLikelihoodTimesPriorGPU = [];
totalRasterizationTime = 0;
totalRasterizationTimeGPU = 0;
figure
for i = 1 : 100

  tic
  [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStep( optimizer );
  elapsedTime = toc;
  disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
  totalRasterizationTime = totalRasterizationTime + elapsedTime;
  disp( [ 'Average rasterization time so far: ' num2str( totalRasterizationTime / i ) ] )
  
  tic
  [ minLogLikelihoodTimesPriorGPU, maximalDeformationGPU ] = kvlDeformOneStepGPU( optimizerGPU );
  elapsedTime = toc;
  disp( [ 'GPU Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
  totalRasterizationTimeGPU = totalRasterizationTimeGPU + elapsedTime;
  disp( [ 'GPU Average rasterization time so far: ' num2str( totalRasterizationTimeGPU / i ) ] )

  historyOfMinLogLikelihoodTimesPrior = [ historyOfMinLogLikelihoodTimesPrior; minLogLikelihoodTimesPrior ];
  historyOfMinLogLikelihoodTimesPriorGPU = [ historyOfMinLogLikelihoodTimesPriorGPU; minLogLikelihoodTimesPriorGPU ];

  hold off
  plot( historyOfMinLogLikelihoodTimesPrior )
  hold on
  plot( historyOfMinLogLikelihoodTimesPriorGPU, 'r' )
  grid
  drawnow
end


