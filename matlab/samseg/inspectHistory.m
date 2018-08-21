%
load history

% Show data we were trying to segment (inspect brain mask)
numberOfContrasts = size( history.imageBuffers, 4 );
for contrastNumber = 1 : numberOfContrasts
  figure
  showImage( history.imageBuffers( :, :, :, contrastNumber ) )
end
  
%
numberOfVoxels = sum( history.mask( : ) );
numberOfMultiResolutionLevels = length( history.historyWithinEachMultiResolutionLevel );


%
figure
%  desiredYLim = [ Inf -Inf ];
numberOfRows = ceil( sqrt( numberOfMultiResolutionLevels ) )
numberOfColumns = ceil( numberOfMultiResolutionLevels / numberOfRows )
for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
  numberOfVoxels = sum( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSampledMask( : ) )
  subplot( numberOfRows, numberOfColumns, multiResolutionLevel )
  %plot( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfCost )
  tmp = history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfCost / numberOfVoxels;
  plot( tmp )
  grid
  title( [ 'multiResolutionLevel: ' num2str( multiResolutionLevel ) ] )
  if 1
    hold on
    ylim = get( gca, 'ylim' );
    thresholds = [ 1e-4 1e-5 1e-6 ];
    thresholdColors = str2mat( 'k', 'b', 'g' );
    for thresholdNumber = 1 : length( thresholds )
      threshold = thresholds( thresholdNumber );
      thresholdColor = thresholdColors( thresholdNumber );
      x = find( ( tmp(1:end-1) - tmp(2:end) ) < threshold, 1 );
      if isempty( x )
        continue
      end
      line( x * [ 1 1 ], ylim, 'color', thresholdColor, 'linestyle', '-.' )
    end
    hold off
  end
  
  %  ylim = get( gca, 'ylim' );
  %  if ( ylim(1) < desiredYLim(1) )
  %    desiredYLim(1) = ylim(1);
  %  end
  %  if ( ylim(2) > desiredYLim(2) )
  %    desiredYLim(2) = ylim(2);
  %  end
  
end
%  for axesHandle = get( gcf, 'children' )
%    set( axesHandle, 'ylim', desiredYLim )
%  end



%
figure
for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
  subplot( 2, numberOfMultiResolutionLevels, multiResolutionLevel )
  timeTakenIntensity = history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenIntensityParameterUpdating;
  timeTakenDeformation = history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenDeformationUpdating;
  timeTakenIntensity = timeTakenIntensity / 60;
  timeTakenDeformation = timeTakenDeformation / 60;
  bar( timeTakenIntensity + timeTakenDeformation, 'r' )
  hold on
  bar( timeTakenDeformation, 'b' )
  grid
  ylabel( 'minutes' )
  l = legend( [ 'Total ' num2str( sum( timeTakenIntensity ) ) 'min' ], [ 'Total ' num2str( sum( timeTakenDeformation ) ) 'min' ] );
  subplot( 2, numberOfMultiResolutionLevels, numberOfMultiResolutionLevels + multiResolutionLevel )
  bar( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfMaximalDeformationApplied )
  grid
  ylabel( 'max. deformation' )
end
%  % Make sure the y-axes are comparable across multi-resolution levels
%  for barTypeNumber = 1 : 2
%    subplot( 2, 2, 1 + 2 * ( barTypeNumber - 1 ) )
%    ylim1 = get( gca, 'ylim' );
%    subplot( 2, 2, 2 + 2 * ( barTypeNumber - 1 ) )
%    ylim2 = get( gca, 'ylim' );
%    ylim = [ min( ylim1(1), ylim2(1) ) max( ylim1(2), ylim2(2) ) ];
%    subplot( 2, 2, 1 + 2 * ( barTypeNumber - 1 )  )
%    set( gca, 'ylim', ylim )
%    subplot( 2, 2, 2 + 2 * ( barTypeNumber - 1 ) )
%    set( gca, 'ylim', ylim )
%  end



%
multiResolutionLevel = 1;
figure
numberOfIterations = length( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration );
numberOfVoxels = sum( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).downSampledMask( : ) );
while true
  for iterationNumber = 1 : numberOfIterations
    subplot( 2, 2, 1 )
    tmp = history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration( iterationNumber ).historyOfDeformationCost / numberOfVoxels;
    plot( tmp )
    if 1
      hold on
      ylim = get( gca, 'ylim' );
      thresholds = [ 1e-4 1e-5 1e-6 ];
      thresholdColors = str2mat( 'k', 'b', 'g' );
      for thresholdNumber = 1 : length( thresholds )
        threshold = thresholds( thresholdNumber );
        thresholdColor = thresholdColors( thresholdNumber );
        x = find( ( tmp(1:end-1) - tmp(2:end) ) < threshold, 1 );
        if isempty( x )
          continue
        end
        line( x * [ 1 1 ], ylim, 'color', thresholdColor, 'linestyle', '-.' )
      end
      hold off
    end
    grid
    ylabel( 'deformation cost' )
    title( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
    subplot( 2, 2, 3 )
    plot( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration( iterationNumber ).historyOfMaximalDeformation )
    grid
    ylabel( 'max. deformation' )
    title( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
    subplot( 2, 2, 2 )
    plot( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration( iterationNumber ).historyOfEMCost / numberOfVoxels )
    grid
    ylabel( 'EM cost' )
    title( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
    drawnow
    disp( 'Press ENTER to continue' )
    pause
  end
end





% Show how prior has deformed
multiResolutionLevel = 1;
figure
numberOfClasses = length( history.input.modelSpecifications.sharedGMMParameters );
colors = 255 * [ hsv( numberOfClasses ) ones( numberOfClasses, 1 ) ];           
DIM = size( history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).downSampledImageBuffers( :, :, :, 1 ) );
downSampledMaskIndices = find( history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).downSampledMask );
dataImage = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).downSampledImageBuffers( :, :, :, 1 );
dataImage = exp( dataImage );
while true
  % for i = [ 1 length( history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).historyOfCost ) ]
  %  tmp = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).historyWithinEachIteration(i).priors;
  for i = 1 : 2
    if ( i == 1 )
      tmp = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).priorsAtStart;
    else
      tmp = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).priorsAtEnd;
    end
    priorImages = zeros( [ DIM numberOfClasses ] ); 
    for j = 1 : numberOfClasses
      priorImages( downSampledMaskIndices + (j-1)*prod( DIM ) ) = tmp( :, j );
    end  
    imageToShow = ( dataImage - min( dataImage(:) ) ) / ( max( dataImage(:) ) - min( dataImage(:) ) );
    alpha = .8;
    imageToShow = ( 1 - alpha ) * repmat( imageToShow, [ 1 1 1 3 ] ) + ...
                  alpha * kvlColorCodeProbabilityImages( priorImages, colors );
    showImage( imageToShow )
    drawnow
    pause( .1 )
  end
end




%
multiResolutionLevel = 1;
classNumber = 5;
mergedName = history.input.modelSpecifications.sharedGMMParameters( classNumber ).mergedName 
modelSpecifications = history.input.modelSpecifications;
numberOfGaussiansPerClass = [ modelSpecifications.sharedGMMParameters.numberOfComponents ];
numberOfClasses = length( numberOfGaussiansPerClass );
numberOfGaussians = sum( numberOfGaussiansPerClass );
downSampledImageBuffer = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).downSampledImageBuffers( :, :, :, 1 );
downSampledImageSize = size( downSampledImageBuffer );
downSampledMaskIndices = find( history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).downSampledMask );
if 1
  % Posterior
  probabilities = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).posteriorsAtEnd;
else
  % Prior
  probabilities = zeros( length( downSampledMaskIndices ), numberOfGaussians );
  priors = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).priorsAtEnd;
  mixtureWeights = ...
          history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).historyWithinEachIteration(end).mixtureWeights;
  for classNumber = 1 : numberOfClasses
    numberOfComponents = numberOfGaussiansPerClass( classNumber );
    for componentNumber = 1 : numberOfComponents
      gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
      probabilities( :, gaussianNumber ) = priors( :, classNumber ) * mixtureWeights( gaussianNumber );
    end
  end
end
numberOfComponents = numberOfGaussiansPerClass( classNumber );
colors = 255 * [ parula( numberOfComponents ) ones( numberOfComponents, 1 ) ];
probabilityImages = zeros( [ downSampledImageSize numberOfComponents ] );
figure
numberOfRows = ceil( sqrt( numberOfComponents ) );
numberOfColumns = ceil( numberOfComponents / numberOfRows );
for componentNumber = 1 : numberOfComponents
  gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
  probabilityImages( downSampledMaskIndices + ( componentNumber-1) * prod( downSampledImageSize ) ) = ...
        probabilities( :, gaussianNumber );
  
  subplot( numberOfRows, numberOfColumns, componentNumber )
  imageToShow = zeros( downSampledImageSize );
  imageToShow( downSampledMaskIndices ) = probabilities( :, gaussianNumber );
  showImage( imageToShow )
  
end
dataImage = ( downSampledImageBuffer - min( downSampledImageBuffer(:) ) ) / ...
              ( max( downSampledImageBuffer(:) ) - min( downSampledImageBuffer(:) ) );
for alpha = [ 1 0.5 ]
  figure
  imageToShow = ( 1 - alpha ) * repmat( dataImage, [ 1 1 1 3 ] ) + ...
                alpha * kvlColorCodeProbabilityImages( probabilityImages, colors );
  showImage( imageToShow )
end
  

  
  
  
%
%
multiResolutionLevel = 1;
modelSpecifications = history.input.modelSpecifications;
numberOfGaussiansPerClass = [ modelSpecifications.sharedGMMParameters.numberOfComponents ];
numberOfClasses = length( numberOfGaussiansPerClass );
numberOfGaussians = sum( numberOfGaussiansPerClass );
downSampledImageBuffer = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).downSampledImageBuffers( :, :, :, 1 );
downSampledImageSize = size( downSampledImageBuffer );
downSampledMaskIndices = find( history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).downSampledMask );
numberOfIterations = ...
   length( history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).historyWithinEachIteration );
binCenters = linspace( min( downSampledImageBuffer(:) ), max( downSampledImageBuffer(:) ), 100 );
numberOfRows = ceil( sqrt( numberOfClasses ) );
numberOfColumns = ceil( numberOfClasses / numberOfRows );
figure
for iterationNumber = 1 : numberOfIterations
  historyWithinIteration = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).historyWithinEachIteration( iterationNumber );
  mixtureWeights = historyWithinIteration.mixtureWeights;
  means = historyWithinIteration.means;
  variances = historyWithinIteration.variances;
  
  for classNumber = 1 : numberOfClasses
    numberOfComponents = numberOfGaussiansPerClass( classNumber );
    mergedName = history.input.modelSpecifications.sharedGMMParameters( classNumber ).mergedName; 
  
    subplot( numberOfRows, numberOfColumns, classNumber )
    model = zeros( size( binCenters ) );
    for componentNumber = 1 : numberOfComponents
      gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;

      mean = means( gaussianNumber );
      variance = variances( gaussianNumber );
      mixtureWeight = mixtureWeights( gaussianNumber );
      gauss = 1/sqrt( 2 * pi * variance ) * exp( -( binCenters - mean ).^2 / 2 / variance );
      plot( binCenters, gauss * mixtureWeight )
      hold on
      model = model + gauss * mixtureWeight;
    end
    plot( binCenters, model )
    hold off
    xlim = get( gca, 'xlim' );
    ylim = get( gca, 'ylim' );
    % t = text( xlim * [ 0.7 0.3 ]', ylim * [ 0.2 0.8 ]', ...
    %          [ 'iterationNumber: ' num2str( iterationNumber ) ] );
    title( [ mergedName ' (iter ' num2str( iterationNumber ) ')' ] )
    grid
  end    
  a = get( gcf, 'Children' );
  xlims = cell2mat( get( a, 'xlim' ) );
  ylims = cell2mat( get( a, 'ylim' ) );
  xlim = [ min( xlims(:,1) ) max( xlims(:,2) ) ];
  ylim = [ min( ylims(:,1) ) max( ylims(:,2) ) ];
  set( a, 'xlim', xlim, 'ylim', ylim )
  hold off
  drawnow
  pause( .1 )
  
end 







