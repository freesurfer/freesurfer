function [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( fileName )
%  
% function [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( fileName )
%
% The format is: FreeSurferLabel compressedLabel name R G B A
%
% Things are sorted according compressedLabels
% 

disp( [ 'Reading contexts of file ' fileName ] )
FreeSurferLabels = [];
compressedLabels = [];
names = [];
colors = [];

fid = fopen( fileName );
if ( fid < 0 )
  error( [ 'Couldn''t read from file ' fileName ] )
end
while 1
  % Get the next line
  textLine = fgetl( fid );
  if ~ischar( textLine ) 
    break
  end

  % Extract the contents
  [ result, count, errorMessage ] = sscanf( textLine, '%d %d %s %d %d %d %d' );
  if ( ~isempty( errorMessage ) )
    error( errorMessage )
  end

  FreeSurferLabel = result( 1 );
  compressedLabel = result( 2 );
  name = ( char( result( 3 : end-4 ) ) )';
  R = result( end-3 );
  G = result( end-2 );
  B = result( end-1 );
  A = result( end );

  % Add contents to output matrices
  FreeSurferLabels = [ FreeSurferLabels; FreeSurferLabel ]; 
  compressedLabels = [ compressedLabels; compressedLabel ];
  names = strvcat( names, name );
  colors = [ colors; [ R G B A ] ];

end
fclose( fid );


% Sort output according to compressedLabel
[ dummy, sortedIndices ] = sort( compressedLabels );
FreeSurferLabels = FreeSurferLabels( sortedIndices, : );
names = names( sortedIndices, : );
colors = colors( sortedIndices, : );



