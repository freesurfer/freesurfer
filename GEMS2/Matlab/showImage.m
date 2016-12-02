function showImage( data, location, range )
%
%

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

xySlice = squeeze( data( :, :, z, : ) );
xzSlice = squeeze( data( :, y, :, : ) );
yzSlice = squeeze( data( x, :, :, : ) );

if ( nargin < 3 )
  range = double( [ min( data(: ) ) max( data(:) ) ] );
end

if ( ndims( data ) < 4 )
  patchedSlices = [ xySlice xzSlice; ...
                    yzSlice' zeros( Nz, Nz ) + range( 1 ) ];
  patchedSlices = ( single( patchedSlices ) - range( 1 ) ) / ( range( 2 ) - range( 1 ) );
  imshow( patchedSlices',  'InitialMagnification', 'fit' )
else
  % RGB data
  patchedSlices = [ xySlice xzSlice; ...
                    permute( yzSlice, [ 2 1 3 ] ) zeros( Nz, Nz, 3 ) ];  
  imshow( permute( patchedSlices, [ 2 1 3 ] ), 'InitialMagnification', 'fit' )

end


