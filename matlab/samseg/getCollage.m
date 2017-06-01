function collage = getCollage( image, numberOfCrossSectionsToShow );
%

%
borderSize = 10; % Measured in pixels

% Get the ingredients: a collage picture for each of the directions
crossSections = cell( 0, 0 );
numberOfRows = 0;
numberOfColumns = 2 * borderSize;
for directionNumber = 1 : 3
  result = getCrossSections( image, directionNumber, numberOfCrossSectionsToShow );
  numberOfRows = max( numberOfRows, size( result, 1 ) );
  numberOfColumns = numberOfColumns + size( result, 2 );
  crossSections{ directionNumber } = result;
end

%  save everything


collage = zeros( numberOfRows, numberOfColumns, 3 ) + .5;
columnStartIndex = 1;
for directionNumber = 1 : 3
  result = crossSections{ directionNumber };
  if ( size( result, 3 ) == 1 )
    % Make an RGB color image if we don't already have one
    result = repmat( result, [ 1 1 3 ] );
  end
  collage( [ 1 : size( result, 1 ) ], columnStartIndex + [ 0 : size( result, 2 )-1 ], : ) = result;
  columnStartIndex = columnStartIndex + size( result, 2 ) + borderSize;
end




function crossSections = getCrossSections( image, dimensionNumber, numberOfCrossSectionsToShow )
%

DIM = size( image );
sliceNumbers = round( [ 1 : numberOfCrossSectionsToShow ] * ...
                      DIM( dimensionNumber ) / ( numberOfCrossSectionsToShow + 1 ) );

% Shuffle dimensions around so that the one you're subsampling from is the first                      
permutationIndices = [ 1 : 4 ];
permutationIndices( 1 ) = dimensionNumber;
permutationIndices( dimensionNumber ) = 1;
permutatedImage = permute( image, permutationIndices );

% Get the slices you need
slices = permutatedImage( sliceNumbers, :, :, : );

% Reshuffle so that subsequent slices are vertically (i.e., row-wise) stacked 
tmp = permute( slices, [ 2 1 3 4 ] );
crossSections = reshape( tmp, [ size( tmp, 1 ) * size( tmp, 2 ) size( tmp, 3 ) size( tmp, 4 ) ] );
