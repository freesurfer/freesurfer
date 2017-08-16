function [ mergedAlphas, mergedNames, mergedFreeSurferLabels, mergedColors, mergingLookupTable ] = kvlMergeAlphas( alphas, names, mergeOptions, FreeSurferLabels, colors );
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
% If the 'mergedName' field contains an existing name (i.e., one that already occurs in variable 'names'), the corresponding FreeSurfer 
% label and color is preserved, and columns in alphas matching the 'searchStrings' are merged into the already existing column. 
% Conversely, if the 'mergedName' field contains a non-existing name, a new column is created that contains the sum of the 'alphas' 
% columns matching the 'searchStrings' - in that case a bogus (negative) "FreeSurfer" label and color is invented.
% 
% The merging is done in first-order-first-serve basis: mergeOptions(1) is executed first, and the result is then fed into mergeOptions(2)
% etc. 
%
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


% Start by copying everything
mergedAlphas = alphas;
mergedNames = names;
mergedFreeSurferLabels = FreeSurferLabels;
mergedColors = colors;

% Also keep track of what class ends up being merged into what
mergingTable = num2cell( [ 1 : size( mergedAlphas, 2 ) ]' );


% 
numberOfMergingTasks = length( mergeOptions );
for mergingTaskNumber = 1 : numberOfMergingTasks

  %
  mergedName = mergeOptions( mergingTaskNumber ).mergedName;
  searchStrings = mergeOptions( mergingTaskNumber ).searchStrings;

  % Retrieve the class number corresponding to the merged label. If it doesn't exist yet,
  % create a new class
  classNumberOfMerged = find( strcmp( cellstr( mergedNames ), mergedName ) );
  if ( length( classNumberOfMerged ) ~= 1 )
    classNumberOfMerged = size( mergedAlphas, 2 ) + 1;
    
    mergedAlphas = [ mergedAlphas zeros( size( mergedAlphas, 1 ), 1 ) ];
    mergedNames = strvcat( mergedNames, mergedName );
    mergedFreeSurferLabels = [ mergedFreeSurferLabels; NaN ]; % We'll take care of these when all is done
    mergedColors = [ mergedColors; NaN NaN NaN NaN ]; % We'll take care of these when all is done
    
    mergingTable = [ mergingTable; cell(1) ];

  end


  % Get class numbers of all structures to be merged into classNumberOfMerged
  classIsNotMerging = ones( size( mergedAlphas, 2 ), 1 );
  classIsNotMerging( classNumberOfMerged ) = 0;
  for searchStringNumber = 1 : length( searchStrings )
    searchString = searchStrings{ searchStringNumber };
    classIsNotMerging = classIsNotMerging .* cellfun( 'isempty', strfind( cellstr( mergedNames ), searchString ) );
  end
  mergingClassNumbers = find( ~classIsNotMerging );
  

  % Now merge (i.e., add the probabilities) of those class numbers to that of classNumberOfMerged
  mergedAlphas( :, classNumberOfMerged ) = sum( mergedAlphas( :, mergingClassNumbers ), 2 );
  mergingTable{ classNumberOfMerged } = [ mergingTable{ mergingClassNumbers' } ];
  
  % Remove the corresponding rows/columns in the resulting variables
  classDisappears = ( ~classIsNotMerging ); classDisappears( classNumberOfMerged ) = 0; 
  disappearingClassNumbers = find( classDisappears );

  disp( [ 'Moved the following structures into ''' mergedName ''':' ] )
  disp( '-------------------------------------------------' )
  disp( mergedNames( disappearingClassNumbers, : ) )
  disp( ' ' )
  
  mergedAlphas( :, disappearingClassNumbers ) = [];                                       
  mergedFreeSurferLabels( disappearingClassNumbers ) = [];
  mergedNames( disappearingClassNumbers, : ) = [];
  mergedColors( disappearingClassNumbers, : ) = [];
  mergingTable( disappearingClassNumbers ) = [];
  
end % End loop over all merging tasks
  
  
% Make sure we still have a valid mesh
if ( max( abs( sum( mergedAlphas, 2 ) - 1 ) ) > 1e-5 )
  error( 'mergedAlphas invalid: class probabilities in the mesh should sum to one in all nodes' )
end


% Take care of NaN's we temporily put in the mergedFreeSurferLabels and mergedColors
classNumberOfNewlyInventedClasses = find( isnan( mergedFreeSurferLabels ) );
numberOfNewlyInvitedClasses = length( classNumberOfNewlyInventedClasses );
mergedFreeSurferLabels( classNumberOfNewlyInventedClasses ) = -[ 1 : numberOfNewlyInvitedClasses ];
mergedColors( classNumberOfNewlyInventedClasses, : ) = ...
      255 * [ hsv( numberOfNewlyInvitedClasses ) ones( numberOfNewlyInvitedClasses, 1 ) ];

      
% Also compute lookup table indicating for each original class number (column number in alphas) the final 
% class number (column number) in mergedAlphas         
numberOfClasses = size( alphas, 2 );
mergingLookupTable = zeros( numberOfClasses, 1 );
for mergedClassNumber = 1 : length( mergingTable )
  mergingLookupTable( mergingTable{ mergedClassNumber } ) = mergedClassNumber;
end
