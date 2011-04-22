function [ f ] = makeVolumePlot( volumes, labels, diseaseStatuses, colors, compareStatistics )
%
%

numberOfLabels = size( labels, 1 );

if ( nargin < 4 )
  colors = zeros( numberOfLabels, 3 );
end
if ( nargin < 5 )
  compareStatistics = 0;
end


%
f = figure;
set( f, 'color', 'w' )
set( f, 'position', [ 0 0 900 600 ] )
legendHandles = [];
stuffToHide = [];
for labelNumber = 1 : numberOfLabels
  % Get the label
  label = deblank( labels( labelNumber, : ) );

  % Get the color
  color = colors( labelNumber, 1 : 3 ) / 255;

  % Decide which disease statuses to show
  diseaseStatusesToShow = [ 0 : 2 ];
  if ( compareStatistics == 1 )
    diseaseStatusesToShow = [ 0 1 ];
  elseif ( compareStatistics == 2 )
    diseaseStatusesToShow = [ 0 2 ];
  elseif ( compareStatistics == 3 )
    diseaseStatusesToShow = [ 1 2 ];
  end

  % Loop over all disease statuses
  volumesToCompare1 = [];
  volumesToCompare2 = [];
  maxYDrawing = 0;
  totalWidthOfThreePatches = 0.8;
  for diseaseStatus = 0 : 2
    % Calculate volume statistics to visualize
    thisVolumes = volumes( labelNumber, find( diseaseStatuses == diseaseStatus ) );
    averageVolume = sum( thisVolumes ) / length( thisVolumes );
    if 0
      minimumVolume = min( thisVolumes );
      maximumVolume = max( thisVolumes );
    else
      variance = sum( ( thisVolumes - averageVolume ).^2 ) / length( thisVolumes );
      minimumVolume = averageVolume - sqrt( variance );
      maximumVolume = averageVolume + sqrt( variance );
    end

    % Draw main patch
    xmin = labelNumber + ( 1 - totalWidthOfThreePatches ) / 2 + diseaseStatus * totalWidthOfThreePatches / 3;
    xmax = xmin + totalWidthOfThreePatches / 3;
    %p = patch( [ xmin xmax xmax xmin ], [ 0 0 averageVolume averageVolume ], color );
    startAlpha = 1;
    endAlpha = 0.5;
    alpha = ( endAlpha - startAlpha ) / 2.0 * diseaseStatus + startAlpha;
    thisColor = color * alpha + [ 1 1 1 ] * ( 1 - alpha );
    p = patch( [ xmin xmax xmax xmin ], [ 0 0 averageVolume averageVolume ], ...
                thisColor );
    set( p, 'EdgeColor', color );
    set( p, 'LineWidth', 3.5 );
    if ( diseaseStatus == 0 )
      legendHandles = [ legendHandles p ];
    end


    % Print disease state
    if ( diseaseStatus == 0 )
      diseaseStatusText = 'NC';
    elseif  ( diseaseStatus == 1 )
      diseaseStatusText = 'MCI';
    else
      diseaseStatusText = 'AD';
    end
    startAlpha = 0;
    endAlpha = 1;
    alpha = ( endAlpha - startAlpha ) / 2.0 * diseaseStatus + startAlpha;
    thisColor = color * alpha + [ 1 1 1 ] * ( 1 - alpha );
    goalPosition = [ ( xmin + xmax ) / 2; averageVolume / 2 ]';
    t = text( goalPosition( 1 ), goalPosition( 2 ), diseaseStatusText, ...
              'rotation', 90, 'fontsize', 10, 'fontweight', 'bold', ...
              'color', thisColor, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle' );
    extent = get( t, 'extent' );
    xShift = goalPosition(1) - ( extent(1) + ( extent(1)+extent(3) ) ) / 2;
    set( t, 'position', goalPosition + [ xShift 0 ] );


    % Draw bars
    xmiddle = ( xmin + xmax ) / 2;
    %l = line( [ xmiddle xmiddle ], [ minimumVolume averageVolume ], 'color', 'w' );
    l1 = line( [ xmiddle xmiddle ], [ minimumVolume averageVolume ], 'color', color );
    l2 = line( [ xmiddle xmiddle ], [ averageVolume maximumVolume ], 'color', color );
    delta = ( xmax - xmin ) / 4;
    %l = line( [ xmiddle-delta xmiddle+delta ], [ minimumVolume minimumVolume  ], 'color', 'w' );
    l3 = line( [ xmiddle-delta xmiddle+delta ], [ minimumVolume minimumVolume  ], 'color', color );
    l4 = line( [ xmiddle-delta xmiddle+delta ], [ maximumVolume maximumVolume  ], 'color', color );


    % Remember some stuff for later
    maxYDrawing = max( maximumVolume, maxYDrawing );
    if ( sum( diseaseStatusesToShow == diseaseStatus ) == 0 )
      stuffToHide = [ stuffToHide; [ p t l1 l2 l3 l4 ]' ];
    else
      if isempty( volumesToCompare1 )
        volumesToCompare1 = thisVolumes;
      else
        volumesToCompare2 = thisVolumes;
      end
    end

  end % End loop over disease status

  if ( compareStatistics ~= 0 )
    %  disp( [ 'label: ' label ] )
    %  length( volumesToCompare1 )
    %  length( volumesToCompare2 )

    %significance = permutationTest( oldThisVolumes', thisVolumes' );
    [ dummy, significance ] =  ttest2( volumesToCompare1', volumesToCompare2' );
    if ( significance < 0.05 )
        textToDisplay = [ 'p < ' sprintf( '%.4f', max( significance, 0.0001 ) ) ];
      t = text( labelNumber + 0.5, maxYDrawing + 0.09 * max( volumes( : ) ), ...
                textToDisplay, ...
                'color', color, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 12 );
      t = text( labelNumber + 0.5, maxYDrawing + 0.05 * max( volumes( : ) ), ...
                '*', ...
                'color', color, ...
                'HorizontalAlignment', 'center', ...
                'FontSize', 16 );
    end

  end

end % End loop over all labels




% Take care of some display details
%title( 'Volumes' )
set( gca, 'xticklabel', '' )
xlim = get( gca, 'xlim' );
set( gca, 'xlim', [ xlim(1)-1 xlim(end)+1 ] )
%  ylim = get( gca, 'ylim' );
%  set( gca, 'ylim', [ ylim(1) ylim(end)+0.1 ] )
%set( get( gca, 'ylabel' ), 'string', 'mm^3' )
grid

% Show labels
labelsToDisplay = [];
for i = 1 : size( labels, 1 )
  label = deblank( labels( i, : ) );
  label = strrep( label, '_', '\_' );
  labelsToDisplay = strvcat( labelsToDisplay, label );
end
l = legend( legendHandles, labelsToDisplay );


% Don't show stuff to hide
ylim = get( gca, 'ylim' );
set( stuffToHide, 'visible', 'off' )
set( gca, 'ylim', ylim )

