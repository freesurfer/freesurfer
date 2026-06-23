function rgbBuffer = kvlColorCodeProbabilityImages( probabilityImages, colors )
%

DIM = size( probabilityImages );
if ( length( DIM ) < 4 )
  numberOfClasses = 1;
else
  numberOfClasses = DIM( 4 );
end  
DIM = DIM( 1 : 3 );

if ( nargin < 2 )
  % Auto-generate some color lookup table
  colors = [ colormap( lines( numberOfClasses ) ) ones( numberOfClasses, 1 ) ] * 255;
  colors( 1, end ) = 128; % Make first class translucent
end


redChannel = zeros( DIM );
greenChannel = zeros( DIM );
blueChannel = zeros( DIM );
for classNumber = 1 : numberOfClasses
  red = colors( classNumber, 1 ) / 255;
  green = colors( classNumber, 2 ) / 255;
  blue = colors( classNumber, 3 ) / 255;
  alpha = colors( classNumber, 4 ) / 255;

  probabilityImage = probabilityImages( :, :, :, classNumber );
  if isa( probabilityImage, 'uint16' )
    probabilityImage = single( probabilityImage ) / 65535;
  end
  redChannel(:) = redChannel(:) + ( red * alpha ) * probabilityImage(:);
  greenChannel(:) = greenChannel(:) + ( green * alpha ) * probabilityImage(:);
  blueChannel(:) = blueChannel(:) + ( blue * alpha ) * probabilityImage(:);

end
rgbBuffer = zeros( [ DIM 3 ] );
rgbBuffer(:,:,:,1) = redChannel;
rgbBuffer(:,:,:,2) = greenChannel;
rgbBuffer(:,:,:,3) = blueChannel;
