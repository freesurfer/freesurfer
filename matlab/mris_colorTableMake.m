function [aS_ct] = mris_colorTableMake( a_labelVar, varargin)
%
% NAME
% function [aS_ct] = mris_colorTableMake( a_labelVar, varargin)
%
%
%
% ARGUMENTS
%       
%       INPUT
%       a_labelVar      var             multi-type:
%                                       + int: number of labels
%                                       + cstr_labelName: cell list of label
%                                                         names
%
%       OPTIONAL
%       ab_labelUnknown bool            if true, add an extra entry at
%                                       [0 0 0 0 0] called 'unknown'.
%
%       OUTPUT
%       aS_ct           struct          color table structure
%
% DESCRIPTION
%
%       'mris_colorTableMake' constructs a color table structure.
%       
%       In the simplest case, its input argument is an integer denoting
%       the size of the table, in which case the script names labels
%       'region-1', 'region-2', ... 'region-N'.
%       
%       Alternatively, a cell array of strings can be passed, in which
%       case these are used as label names and also define the size
%       of the table.
%              
% PRECONDITIONS
%       o FreeSurfer environment.
%
% POSTCONDITIONS
%       o FreeSurfer colortable struct is returned.
%         
% SEE ALSO
%       o read_annotation / write_annotation for a description of the
%         color table format.
%
% HISTORY
% 07 June 2011
% o Initial design and coding.
%

% ---------------------------------------------------------

sys_printf('mris_colorTableMake: START\n');
 

% Parse optional arguments
b_labelUnknown  = 0;
if length(varargin),    b_labelUnknown  = varargin{1}; end
if b_labelUnknown,      b_labelUnknown  = 1; end

% Determine the table size
b_labelNames    = 0;
if (isfloat(a_labelVar))
    tableSize   = int32(a_labelVar);
end
if (iscell(a_labelVar))
    tableSize   = numel(a_labelVar);
    b_labelNames = 1;
end

% Build the color table using permutations in the [256 256 256] space
dimension       = 0;
rows            = 0;
while rows < tableSize
    dimension   = dimension + 1;
    [I, D]      = permutations_find(3, dimension, [127 127 127]);
    [rows cols] = size(I{dimension});
end

M_RGBfull       = int32(V_normalize(I{dimension})*255);
% Skip the first _RGBfull entry which is [ 0 0 0 ]
M_RGB           = M_RGBfull(2:tableSize+1, :);

v_index         = [tableSize,1];
% This anonymous function calculates the color index based on
% R + G*2^8 + B*2^16
fRGB            = @(X) X(:,1) + X(:,2)*2^8 + X(:,3) * 2^16;
v_index         = fRGB(M_RGB);
M_ct            = double([M_RGB zeros(tableSize, 1) v_index]);

% Label names
if (~b_labelNames)
    Cstr_labelName = cell(tableSize, 1);
    for i=1:tableSize
        Cstr_labelName{i}       = sprintf('region-%d', i);
    end
else
    Cstr_labelName              = a_labelVar;
end

if b_labelUnknown && ~b_labelNames
    Cstr_labelName      = [ 'unknown' Cstr_labelName' ]';
    v_unknown           = [0 0 0 0 0];
    M_ct                = [v_unknown' M_ct']';
end

aS_ct                   = struct;
aS_ct.numEntries        = tableSize + b_labelUnknown;
aS_ct.orig_tab          = 'none';
aS_ct.struct_names      = Cstr_labelName;
aS_ct.table             = M_ct;

sys_printf('mris_colorTableMake: END\n');

end
% ---------------------------------------------------------


