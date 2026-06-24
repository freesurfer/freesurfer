function [ mergedAlphas, mergedNames, mergedFreeSurferLabels, mergedColors, translationTable ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
%
%  [ mergedAlphas, mergedNames, [ mergedFreeSurferLabels, mergedColors ] ] = ...
%          kvlMergeAlphas( alphas, names, mergeOptions, [ FreeSurferLabels, colors ] );
%
% where mergeOptions is of the form
%  
%    mergeOptions = struct;
%    mergeOptions(1).mergedName = 'Global WM';
%    mergeOptions(1).searchStrings = { 'White', 'Brain-Stem' };     
%    mergeOptions(2).mergedName = 'Global GM';
%    mergeOptions(2).searchStrings = { 'Cortex', 'Caudate', 'Hippocampus' };     
%    mergeOptions(3).mergedName = 'Unknown';
%    mergeOptions(3).searchStrings = { 'CSF', '3', '4' };
%
% creates a 'mergedAlphas' matrix where one or more columns of the 'alphas' matrix have been "merged" (i.e., added together).
%  


%  if 0
%    save xxx.mat alphas names mergeOptions FreeSurferLabels colors
%  
%    exit;
%  else  
%    close all
%    clear all
%    addpath /home/koen/software/freesurferGit/freesurfer/GEMS2/Matlab/
%    load xxx.mat
%  end
%  


%
if ( nargin < 4 )
  FreeSurferLabels = nan( size( alphas, 2 ), 1 );
  colors = nan( size( alphas, 2 ), 4 );
end


% Make sure we're dealing with a valid mesh
if ( max( abs( sum( alphas, 2 ) - 1 ) ) > 1e-5 )
  error( 'alphas invalid: class probabilities in the mesh should sum to one in all nodes' )
end




%
numberOfClasses = length( mergeOptions );  % Super-structures aka mixture models
numberOfStructures = size( alphas, 2 ); % Number of freesurfer labels being segmented/outputed
translationTable = zeros( numberOfClasses, numberOfStructures );
mergedNames = [];
for classNumber = 1 : numberOfClasses
  % 
  className = mergeOptions( classNumber ).mergedName;
  mergedNames = strvcat( mergedNames, className );
  
  %
  searchStrings = mergeOptions( classNumber ).searchStrings;
  for searchStringNumber = 1 : length( searchStrings )
    searchString = searchStrings{ searchStringNumber };
    structureNumbers = find( ~cellfun( 'isempty', strfind( cellstr( names ), searchString ) ) );
    translationTable( classNumber, structureNumbers ) = 1;
  end

end
if sum( sum(translationTable) == 0 )
  error( 'Some structures are not associated with any super-structures' )
end
translationTable = translationTable ./ sum( translationTable );



%
mergedAlphas = alphas * translationTable';


%
mergedFreeSurferLabels = -[ 1 : numberOfClasses ]';
mergedColors = 255 * [ hsv( numberOfClasses ) ones( numberOfClasses, 1 ) ];


% Print out 
for classNumber = 1 : numberOfClasses
  % 
  disp( mergedNames( classNumber, : ) )
  for structureNumber = 1 : numberOfStructures
    percentage = translationTable( classNumber, structureNumber );
    if ( percentage > 0 )
      disp( [ '    ' names( structureNumber, : ) ' (' num2str( percentage * 100 )  '%)' ] )
    end    
  end
     
end 


