function mosaicedImage = mosaicImages( image1, image2, numberOfBlocks, rangeImage1, rangeImage2 )
%
%

if nargin < 5
  rangeImage2 = [ min( image2(:) ) max( image2(:) ) ];
end

if nargin < 4
  rangeImage1 = [ min( image1(:) ) max( image1(:) ) ];
end

if nargin < 3
  numberOfBlocks = 3;
end

% Normalize images
image1 = ( image1 - rangeImage1( 1 ) ) / ( rangeImage1( 2 ) - rangeImage1( 1 ) );
image2 = ( image2 - rangeImage2( 1 ) ) / ( rangeImage2( 2 ) - rangeImage2( 1 ) );

% Create 3-D repetitive pattern using seperable basis functions
DIM = size( image1 );
patternX = sin( [ 0 : DIM( 1 )-1 ] / DIM( 1 ) * 2 * pi * ( numberOfBlocks / 2 ) );
patternY = sin( [ 0 : DIM( 2 )-1 ] / DIM( 2 ) * 2 * pi * ( numberOfBlocks / 2 ) );
patternZ = sin( [ 0 : DIM( 3 )-1 ] / DIM( 3 ) * 2 * pi * ( numberOfBlocks / 2 ) );

pattern = zeros( DIM );
tmp = patternX' * patternY; 
for sliceNumber = 1 : DIM( 3 )
  pattern( :, :, sliceNumber ) = tmp * patternZ( sliceNumber );
end

mosaicedImage = image1;
image2Indices = find( pattern < 0 ); 
mosaicedImage( image2Indices ) = image2( image2Indices );




