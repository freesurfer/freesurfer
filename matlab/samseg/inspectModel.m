function inspectModel( historyFilename )
%

if ( nargin == 0 )
  % Test
  addpath /home/koen/software/freesurferGit/freesurfer/matlab/samseg/
  addpath /home/koen/software/freesurferGit/freesurfer/GEMS-Release/bin/

  historyFilename = '/home/koen/data/testing/bigvent30/bert/history.mat'
  inspectModel( historyFilename );
  
  return
end


% Options here - global variables 
m_MultiResolutionLevel = 0;
m_ClassNumber = 0;
m_ComponentNumber = 0;
m_ShowPosterior = 0;
m_Location = [];
m_Alpha = 0.5;

% Cache here
m_NumberOfComponents = [];
m_CurrentLevelHistory = [];
m_DownSampledImageBuffer = [];
m_DownSampledImageSize = [];
m_DownSampledMaskIndices = [];
m_Probabilities = [];



%
load( historyFilename, 'history' )
modelSpecifications = history.input.modelSpecifications;
numberOfGaussiansPerClass = [ modelSpecifications.sharedGMMParameters.numberOfComponents ];
numberOfClasses = length( numberOfGaussiansPerClass );
numberOfGaussians = sum( numberOfGaussiansPerClass );


% Construct GUI elements
screenSize = get( 0, 'ScreenSize' );
figure( 'position', [ 100 100 screenSize(3) * .9 screenSize(4) * .9 ] )
%  figure

a1 = axes( 'position', [ .01 .26 .48 .73 ] );
%  location = showImage( randn(200,100,50) );


numberOfRows = ceil( sqrt( numberOfClasses ) );
numberOfColumns = ceil( numberOfClasses / numberOfRows );
targetPosition = [ .52 .05 .46 .93 ];
fillFactor = .95; % Leave a bit of empty space
as = [];
for classNumber = 1 : numberOfClasses
  columnNumber = floor( (classNumber-1)/numberOfRows ) + 1
  rowNumber = classNumber - ( columnNumber-1 ) * numberOfRows
  a = axes( 'position', ...
            [ targetPosition(1)+(columnNumber-1)*targetPosition(3)/numberOfColumns ...
              targetPosition(2)+(rowNumber-1)*targetPosition(4)/numberOfRows ...
              targetPosition(3)/numberOfColumns*fillFactor targetPosition(4)/numberOfRows*fillFactor ] );
  % get( a, 'position' )            
  plot( randn( 10, 1 ) )
  as = [ as a ];            
end

mergedNames = cell( 1, numberOfClasses );
for classNumber = 1 : numberOfClasses
  mergedName = history.input.modelSpecifications.sharedGMMParameters( classNumber ).mergedName;
  mergedNames{ classNumber } = mergedName;
end 
h1 = uicontrol( 'Style', 'popupmenu', ...
               'String', mergedNames, ...
               'Units', 'normalized', ...
               'Position', [ .01 .00 .3 .2 ], ...
               'Callback', @classNumberCallback );

h2 = uicontrol( 'Style', 'popupmenu', ...
               'String', {}, ...
               'Units', 'normalized', ...
               'Position', [ .01 -.05 .3 .2 ], ...
               'Callback', @componentNumberCallback );

h3 = uicontrol( 'Style', 'popupmenu', ...
               'String', { 'Prior', 'Posterior' }, ...
               'Units', 'normalized', ...
               'Position', [ .01 -.10 .3 .2 ], ...
               'Callback', @posteriorCallback );
               
numberOfMultiResolutionLevels = length( history.historyWithinEachMultiResolutionLevel );
multiResolutionNames = cell( 1, numberOfMultiResolutionLevels );
for level = 1 : numberOfMultiResolutionLevels
  multiResolutionNames{ level } = [ 'multiResolutionLevel: ' num2str( level ) ];
end
h4 = uicontrol( 'Style', 'popupmenu', ...
               'String', multiResolutionNames, ...
               'Units', 'normalized', ...
               'Position', [ .01 -.15 .3 .2 ], ...
               'Callback', @multiResolutionLevelCallback );

sld = uicontrol( 'Style', 'slider', 'Min', 0, 'Max', 1, ...
                 'Value', m_Alpha, ...
                 'String', 'Alpha', ...
                 'Units', 'normalized', ...
                 'Position', [ .01 .21 .3 .02 ], ...
                 'Callback', @alphaCallback );
               
ls = []; % Vertical line handles indicating intensity level in histogram












% Pretend someone's pressing
h1.Value = 1;
h3.Value = 2;
h4.Value = numberOfMultiResolutionLevels;
%  sld.Value = m_Alpha;
multiResolutionLevelCallback( h4 );



%
function alphaCallback( source, eventdata )
  % 
  m_Alpha = sld.Value;

  % 
  componentNumberCallback( h2 );

end




% Callback for multiResolutionLevel
function multiResolutionLevelCallback( source, eventdata ) 
  % 
  m_MultiResolutionLevel = source.Value;

  m_CurrentLevelHistory = history.historyWithinEachMultiResolutionLevel( m_MultiResolutionLevel );
  m_DownSampledImageBuffer =  m_CurrentLevelHistory.downSampledImageBuffers( :, :, :, 1 );
  m_DownSampledImageSize = size( m_DownSampledImageBuffer );
  m_DownSampledMaskIndices = find( m_CurrentLevelHistory.downSampledMask );
  m_Location = round( m_DownSampledImageSize/2 );

  %
  binCenters = linspace( min( m_DownSampledImageBuffer(:) ), max( m_DownSampledImageBuffer(:) ), 100 );
  means = m_CurrentLevelHistory.historyWithinEachIteration(end).means;
  variances = m_CurrentLevelHistory.historyWithinEachIteration(end).variances;
  mixtureWeights = m_CurrentLevelHistory.historyWithinEachIteration(end).mixtureWeights;
  for classNumber = 1 : numberOfClasses
    numberOfComponents = numberOfGaussiansPerClass( classNumber );
    %mergedName = history.input.modelSpecifications.sharedGMMParameters( classNumber ).mergedName; 

    axes( as( classNumber ) )
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
    %set( gca, 'yscale', 'log' )
    hold off
    xlim = get( gca, 'xlim' );
    ylim = get( gca, 'ylim' );
    % t = text( xlim * [ 0.7 0.3 ]', ylim * [ 0.2 0.8 ]', ...
    %          [ 'iterationNumber: ' num2str( iterationNumber ) ] );
    %title( mergedName )
    grid
  end    
  xlims = cell2mat( get( as, 'xlim' ) );
  ylims = cell2mat( get( as, 'ylim' ) );
  xlim = [ min( xlims(:,1) ) max( xlims(:,2) ) ];
  ylim = [ min( ylims(:,1) ) max( ylims(:,2) ) ];
  set( as, 'xlim', xlim, 'ylim', ylim )

  %
  for classNumber = 1 : numberOfClasses
    axes( as( classNumber ) )
    mergedName = mergedNames{ classNumber }
    t = text( xlim(1) + ( xlim(2) - xlim(1) ) / 10, ...
              ylim(2) - ( ylim(2) - ylim(1) ) / 10, ...
              mergedName, 'Interpreter', 'none' )
    hold off
  end

  
  %
  posteriorCallback( h3 );
  imageClickCallback( m_DownSampledImageBuffer, m_Location );
  
end



% Callback for showPosterior
function posteriorCallback( source, eventdata ) 
  % 
  if ( source.Value == 2 )
    m_ShowPosterior = true;
  else
    m_ShowPosterior = false;  
  end
    
  if m_ShowPosterior
    % Posterior
    m_Probabilities = m_CurrentLevelHistory.posteriorsAtEnd;
  else
    % Prior
    m_Probabilities = zeros( length( m_DownSampledMaskIndices ), numberOfGaussians );
    priors = m_CurrentLevelHistory.priorsAtEnd;
    mixtureWeights = m_CurrentLevelHistory.historyWithinEachIteration(end).mixtureWeights;
    for classNumber = 1 : numberOfClasses
      numberOfComponents = numberOfGaussiansPerClass( classNumber );
      for componentNumber = 1 : numberOfComponents
        gaussianNumber = sum( numberOfGaussiansPerClass( 1 : classNumber-1 ) ) + componentNumber;
        m_Probabilities( :, gaussianNumber ) = priors( :, classNumber ) * mixtureWeights( gaussianNumber );
      end
    end
  end

  %
  classNumberCallback( h1 );

end 



function classNumberCallback( source, eventdata ) 
  % 
  m_ClassNumber = source.Value
  
  % Populate the componentNumber GUI
  m_NumberOfComponents = numberOfGaussiansPerClass( m_ClassNumber )
  componentNames = cell( 1, m_NumberOfComponents+1 );
  for componentNumber = 1 : m_NumberOfComponents
    componentNames{ componentNumber } = [ 'Component ' num2str( componentNumber ) ];
  end 
  componentNames{ end } = 'Composite';
  h2.String = componentNames;
  
  % 
  h2.Value = m_NumberOfComponents+1;
  componentNumberCallback( h2 );

  %
  set( as, 'color', [1 1 1 ] )
  set( as( m_ClassNumber ), 'color', [ .9 1 .9 ] )

end


function componentNumberCallback( source, eventdata ) 
  % 
  m_ComponentNumber = source.Value
  m_NumberOfComponents

  %
  if ( m_ComponentNumber <= m_NumberOfComponents )
    % Show one component
    disp( '=========' )
    m_ComponentNumber
    m_NumberOfComponents
    disp( '=========' )    
    colors = [ 0 0 255 255 ];
    gaussianNumber = sum( numberOfGaussiansPerClass( 1 : m_ClassNumber-1 ) ) + m_ComponentNumber;
    probabilityImages = zeros( [ m_DownSampledImageSize 1 ] );
    probabilityImages( m_DownSampledMaskIndices ) = m_Probabilities( :, gaussianNumber );
  else
    % Show all components simultaniously
    colors = 255 * [ parula( m_NumberOfComponents ) ones( m_NumberOfComponents, 1 ) ];
    probabilityImages = zeros( [ m_DownSampledImageSize m_NumberOfComponents ] );
    for componentNumber = 1 : m_NumberOfComponents
      gaussianNumber = sum( numberOfGaussiansPerClass( 1 : m_ClassNumber-1 ) ) + componentNumber;
      probabilityImages( m_DownSampledMaskIndices + ( componentNumber-1) * prod( m_DownSampledImageSize ) ) = ...
            m_Probabilities( :, gaussianNumber );
    end
  end  
  % dataImage = m_DownSampledImageBuffer;
  dataImage = exp( m_DownSampledImageBuffer );
  dataImage = ( dataImage - min( dataImage(:) ) ) / ...
                ( max( dataImage(:) ) - min( dataImage(:) ) );
  colorCodedProbabilityImages = kvlColorCodeProbabilityImages( probabilityImages, colors );
  axes( a1 )
  imageToShow = ( 1 - m_Alpha ) * repmat( dataImage, [ 1 1 1 3 ] ) + ...
                  m_Alpha * colorCodedProbabilityImages;
  m_Location = showImage( imageToShow, m_Location, [], @imageClickCallback );
    
end




function imageClickCallback( data, location )
  %
  m_Location = location;
  
  %
  intensity = m_DownSampledImageBuffer( m_Location(1), m_Location(2), m_Location(3) )

  %
  delete( ls )
  ls = [];
  for classNumber = 1 : numberOfClasses
    axes( as( classNumber ) )
    hold on
    ylim = get( as( classNumber ), 'ylim' );
    l = line( intensity * [1 1 ], ylim, 'color', 'k', 'linestyle', '--' );
    ls = [ ls l ];
    hold off
  end
  
end





end


