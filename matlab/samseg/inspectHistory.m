%
load history

% Show data we were trying to segment (inspect brain mask)
dataImage = zeros( history.DIM );
dataImage( history.maskIndices ) = history.data;
figure
showImage( dataImage )


%
figure
for multiResolutionLevel = 1 : 2
  subplot( 1, 2, multiResolutionLevel )
  plot( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfCost )

end


%
figure
for multiResolutionLevel = 1 : 2
  subplot( 2, 2, 1 + ( multiResolutionLevel - 1 ) )
  timeTakenIntensity = history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenIntensityParameterUpdating;
  timeTakenDeformation = history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfTimeTakenDeformationUpdating;
  timeTakenIntensity = timeTakenIntensity / 60;
  timeTakenDeformation = timeTakenDeformation / 60;
  bar( timeTakenIntensity + timeTakenDeformation, 'r' )
  hold on
  bar( timeTakenDeformation, 'b' )
  grid
  ylabel( 'minutes' )
  subplot( 2, 2, 3 + ( multiResolutionLevel - 1 ) )
  bar( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyOfMaximalDeformationApplied )
  grid
  ylabel( 'max. deformation' )
end
% Make sure the y-axes are comparable across multi-resolution levels
for barTypeNumber = 1 : 2
  subplot( 2, 2, 1 + 2 * ( barTypeNumber - 1 ) )
  ylim1 = get( gca, 'ylim' );
  subplot( 2, 2, 2 + 2 * ( barTypeNumber - 1 ) )
  ylim2 = get( gca, 'ylim' );
  ylim = [ min( ylim1(1), ylim2(1) ) max( ylim1(2), ylim2(2) ) ];
  subplot( 2, 2, 1 + 2 * ( barTypeNumber - 1 )  )
  set( gca, 'ylim', ylim )
  subplot( 2, 2, 2 + 2 * ( barTypeNumber - 1 ) )
  set( gca, 'ylim', ylim )
end



%
figure
multiResolutionLevel = 1;
numberOfIterations = length( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration );
while true
  for iterationNumber = 1 : numberOfIterations
    subplot( 2, 2, 1 )
    plot( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration( iterationNumber ).historyOfDeformationCost )
    grid
    ylabel( 'deformation cost' )
    title( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
    subplot( 2, 2, 3 )
    plot( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration( iterationNumber ).historyOfMaximalDeformation )
    grid
    ylabel( 'max. deformation' )
    title( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
    subplot( 2, 2, 2 )
    plot( history.historyWithinEachMultiResolutionLevel( multiResolutionLevel ).historyWithinEachIteration( iterationNumber ).historyOfEMCost )
    grid
    ylabel( 'EM cost' )
    title( [ 'iterationNumber: ' num2str( iterationNumber ) ] )
    drawnow
    disp( 'Press ENTER to continue' )
    pause
  end
end





% Show how prior has deformed
figure
colors = [  128 0 128 0; ...
           0 255 0 255; ...
           0 0 255 255; ...
           255 0 0 255 
           128 255 128 255;
           128 128 255 255;
           255 128 128 255 ];
multiResolutionLevel = 2;
while true
  for i = [ 1 length( history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).historyOfCost ) ]
    tmp = history.historyWithinEachMultiResolutionLevel(multiResolutionLevel).historyWithinEachIteration(i).priors;
    priorImages = zeros( [ history.DIM 7 ] ); 
    for j = 1 : 7
      priorImages( history.maskIndices + (j-1)*prod( history.DIM ) ) = tmp( :, j );
    end  
    imageToShow = ( dataImage - min( dataImage(:) ) ) / ( max( dataImage(:) ) - min( dataImage(:) ) );
    alpha = .4;
    imageToShow = ( 1 - alpha ) * repmat( imageToShow, [ 1 1 1 3 ] ) + ...
                  alpha * kvlColorCodeProbabilityImages( priorImages, colors );
    showImage( imageToShow )
    drawnow
    %pause( .1 )
  end
end
 
 
