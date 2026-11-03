function location = showImage( data, location, range, callback )
%
%
if ( nargin == 0 )
  % Test run
  data = randn( 200, 100, 50 ) + [ 1 : 200 ]';
  % showImage( data, [], [], @() disp('hello') );
  myCallback = @(data,location) disp( [ 'data(' num2str(location(1)) ',' num2str(location(2)) ',' num2str(location(3)) '): ' num2str( data( location(1), location(2), location(3) ) ) ] );
  showImage( data, [], [], myCallback );
  
  
  
  return
end


if ( nargin < 4 )
  callback = [];
end

if ( isa( data, 'int64' ) && isscalar( data ) )
  data = kvlGetImageBuffer( data );
end

Nx = size( data, 1 );
Ny = size( data, 2 );
Nz = size( data, 3 );

x = round( Nx / 2 );
y = round( Ny / 2 );
z = round( Nz / 2 );

if ( ( nargin >= 2 ) & ~isempty( location ) )
  x = location( 1 );
  y = location( 2 );
  z = location( 3 );
end

if ( ( nargin < 3 ) | isempty( range ) )
  range = double( [ min( data(: ) ) max( data(:) ) ] );
end

%
[ slicesHandle, lineHandles ] = showPosition( x, y, z );





function [ slicesHandle, lineHandles ] = showPosition( x, y, z )
%

% Show image data
xySlice = squeeze( data( :, :, z, : ) );
xzSlice = squeeze( data( :, y, :, : ) );
yzSlice = squeeze( data( x, :, :, : ) );

if ( ndims( data ) < 4 )
  patchedSlices = [ xySlice xzSlice; ...
                    yzSlice' zeros( Nz, Nz ) + range( 1 ) ];
  patchedSlices = ( single( patchedSlices ) - range( 1 ) ) / ( range( 2 ) - range( 1 ) );
  slicesHandle = imshow( patchedSlices',  'InitialMagnification', 'fit' );
else
  % RGB data
  patchedSlices = [ xySlice xzSlice; ...
                    permute( yzSlice, [ 2 1 3 ] ) zeros( Nz, Nz, 3 ) ];  
  slicesHandle = imshow( permute( patchedSlices, [ 2 1 3 ] ), 'InitialMagnification', 'fit' );

end


% Also show cross-hairs. Remember: size layout is a follows:
%
%   [ Ny x Nx ] [ Ny x Nz ]
%   [ Nz x Nx ] [ Nz x Nz ]
%
l1 = line( [ 1 ( Nx + Nz ) ], y * [ 1 1 ] );
l2 = line( x * [ 1 1 ], [ 1 ( Ny + Nz ) ] );
l3 = line( [ 1 Nx ], ( Ny + z ) * [ 1 1 ] );
l4 = line( ( Nx + z ) * [ 1 1 ], [ 1 Ny ] );
lineHandles = [ l1; l2; l3; l4 ];
set( lineHandles, 'color', 'r' )

% Add callback
set( slicesHandle, 'ButtonDownFcn', @myImageClickCallback )

location = [ x y z ];


end



function myImageClickCallback( objectHandle , eventData )
%

% Retrieve coordinate clicked
axesHandle  = get( objectHandle, 'Parent' );
location = get( axesHandle, 'CurrentPoint' ); location = round( location( 1, [ 1 2 ] ) );
if ( location( 1 ) <= Nx )
  x = location( 1 );
  if ( location (2 ) <= Ny )
    y = location( 2 );
  else
    z = location( 2 ) - Ny;
  end

else
  if ( location (2 ) <= Ny )
    z = location( 1 ) - Nx;
    y = location( 2 );
  end
end

delete( lineHandles );
[ slicesHandle, lineHandles ] = showPosition( x, y, z );

% 
if ~isempty( callback )
  disp( 'doing callback' )
  % callback
  feval( callback, data, location )
end


end


end






