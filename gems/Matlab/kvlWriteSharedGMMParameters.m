function kvlWriteSharedGMMParameters( sharedGMMParameters, fileName )
%
% function kvlWriteSharedGMMParameters( sharedGMMParameters, fileName )
%
% Write text file where each line is formatted as:
%  
%  mergedName numberOfComponents searchString(s)
%
% from input of the form
%  
%    sharedGMMParameters = struct;
%    sharedGMMParameters( 1 ).mergedName = 'Unknown';
%    sharedGMMParameters( 1 ).searchStrings = { 'Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
%    sharedGMMParameters( 1 ).numberOfComponents = 1;
%    sharedGMMParameters( 2 ).mergedName = 'GlobalWM';
%    sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
%    sharedGMMParameters( 2 ).numberOfComponents = 1;
%

if ( nargin == 0 )
  % Test ourselves
  fileName = '/tmp/sharedGMMParameters.txt';
  
  sharedGMMParameters = struct;
  sharedGMMParameters( 1 ).mergedName = 'Unknown';
  sharedGMMParameters( 1 ).searchStrings = { 'Unknown', 'Ventricle', 'Inf-Lat-Vent', 'CSF', 'vessel', 'choroid-plexus' };
  sharedGMMParameters( 1 ).numberOfComponents = 1;
  sharedGMMParameters( 2 ).mergedName = 'GlobalWM';
  sharedGMMParameters( 2 ).searchStrings = { 'White', 'Brain-Stem', 'VentralDC', 'Optic-Chiasm' };
  sharedGMMParameters( 2 ).numberOfComponents = 1;
  
  kvlWriteSharedGMMParameters( sharedGMMParameters, fileName );
  
  return
end



disp( [ 'Writing file ' fileName ] )

fid = fopen( fileName, 'w' );
if ( fid < 0 )
  error( [ 'Couldn''t write to file ' fileName ] )
end
fprintf( fid, '# The format is: mergedName numberOfComponents searchStrings\n\n' );
numberOfStructures = length( sharedGMMParameters );
for structureNumber = 1 : numberOfStructures
  textLine = [ sharedGMMParameters( structureNumber ).mergedName ' ' ...
               num2str( sharedGMMParameters( structureNumber ).numberOfComponents ) ];
               
  searchStrings = sharedGMMParameters( structureNumber ).searchStrings;
  for searchStringNumber = 1 : length( searchStrings )
    searchString = searchStrings{ searchStringNumber };
    textLine = [ textLine ' ' searchString ];
  end

  fprintf( fid, [ textLine '\n' ] );
end
fclose( fid );


return

