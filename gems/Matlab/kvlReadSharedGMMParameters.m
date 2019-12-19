function sharedGMMParameters = kvlReadSharedGMMParameters( fileName )
%
% function sharedGMMParameters = kvlReadSharedGMMParameters( fileName )
%
% Read text file where each line is formatted as:
%  
%  mergedName numberOfComponents searchString(s)
%

if ( nargin == 0 )
  % Test ourselves
  fileName = '/data/testing/atlas/koenAtlases/10SubjectsSmoothing/sharedGMMParameters.txt';
  sharedGMMParameters = kvlReadSharedGMMParameters( fileName );
  return
end



disp( [ 'Reading contexts of file ' fileName ] )
sharedGMMParameters = struct( [] );

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
  
  % Remove leading and trailing blanks
  textLine = strtrim( textLine ); 

  if ( ( length( textLine ) == 0 ) | ( textLine(1) == '#' ) )
    % Skip comment and empty lines 
    continue;
  end

  % Extract the contents
  whiteSpaceIndices = findstr( textLine, ' ' );
  if ( length( whiteSpaceIndices ) < 2 )
    error( [ 'Can''t parse line: ' textLine ] )
  end
  mergedName = textLine( 1 : whiteSpaceIndices(1) );
  numberOfComponents = str2num( textLine( whiteSpaceIndices(1)+1 : whiteSpaceIndices(2) ) );
  textLine = textLine( whiteSpaceIndices(2)+1 : end );
  searchStrings = ( textscan( textLine, '%s' ) );                           
  searchStrings = searchStrings{ : }';
  
  % Add to sharedGMMParameters structure array
  index = length( sharedGMMParameters ) + 1;
  sharedGMMParameters( index ).mergedName = mergedName;
  sharedGMMParameters( index ).numberOfComponents = numberOfComponents;
  sharedGMMParameters( index ).searchStrings = searchStrings;
  
end % End loop over all lines

fclose( fid );


return

