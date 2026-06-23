function kvlWriteCompressionLookupTable( fileName, FreeSurferLabels, names, colors )
%
% function kvlWriteCompressionLookupTable( fileName, FreeSurferLabels, names, colors )
%
% The format is: FreeSurferLabel compressedLabel name R G B A
%
% Lines are sorted according to the FreeSurferLabel
%  

if ( nargin == 0 )
  % Test ourselves
  fileName = '/tmp/compressionLookupTable.txt';
  
  [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( ...
    '/data/testing/atlas/Buckner39AtlasWithMoreClassesAndEyeballs/atlases/10SubjectAtlas3X/result/compressionLookupTable.txt' );
  
  kvlWriteCompressionLookupTable( fileName, FreeSurferLabels, names, colors );
  
  return
end


% Sort according to FreeSurferLabels
[ dummy, sortedIndices ] = sort( FreeSurferLabels );
FreeSurferLabels = FreeSurferLabels( sortedIndices, : );
names = names( sortedIndices, : );
colors = colors( sortedIndices, : );
compressedLabels = sortedIndices - 1;

%
disp( [ 'Writing file ' fileName ] )
fid = fopen( fileName, 'w' );
if ( fid < 0 )
  error( [ 'Couldn''t write to file ' fileName ] )
end
for labelNumber = 1 : size( FreeSurferLabels, 1 )
  FreeSurferLabel = FreeSurferLabels( labelNumber );
  compressedLabel = compressedLabels( labelNumber );
  name = names( labelNumber, : );
  color = colors( labelNumber, : ); 
  
  textLine = [ num2str( FreeSurferLabel ) ' ' num2str( compressedLabel ) ' ' name ' ' num2str( color ) ];             

  fprintf( fid, [ textLine '\n' ] );
end
fclose( fid );

return

