function [hf, hp, av_curv, av_filtered] = mris_display( astr_mris, ...
                                                        astr_curv, ...
                                                        varargin)
%
% NAME
%       [hf, hp, av_curv, av_filtered]   =              ...
%       mris_display(     	astr_mris,              ... 
%                               astr_curv               ...
%                               <, astr_title,          ...
%                               a_az,                   ...
%                               a_el,                   ...
%                               a_bandFilter,           ...
%                               ab_invMap,              ...
%                               astr_colorMap
%                         )
%
%
% ARGUMENTS
%       
%       INPUT
%       astr_mris       string          filename of surface file to load
%       astr_curv       string          filename of curvature file to load
%
%       OPTIONAL
%       astr_title      string          title of plot. If empty string,
%                                       title will be constructed from
%                                       surface and curvature file names
%       a_az            float           azimuth of viewpoint
%       a_el            float           elevation of viewpoint
%       a_bandFilter    variable        apply a bandpass filter over data:
%                                           if vector, defines range
%                                           if single float, defines
%                                              symmetrical pulse [-f, +f]
%                                           if string == 'gz'
%                                              define range as [>0, max] 
%                                       between vector range.
%       ab_invMap       bool            if true, invert the sign of the
%                                       curvature data -- useful if 'neg'
%                                       values are on gyri and 'pos' values
%                                       in sulci.
%       astr_colorMap   string          colormap override to use.
%
%       OUTPUT
%       hf              handle          handle to generated figure
%       hp              handle          handle to patch
%       av_curv         vector          curvature vector
%       av_filtered     vector          number of vertices in original curv
%                                       that have been filtered out by the
%                                       <af_bandFilter>.
%
% DESCRIPTION
%
%       'mris_display' reads a FreeSurfer surface structure and
%       displays it as seen from an optional <a_az, a_el>.
%
%       The <av_curv> can be further band pass filtered by passing
%	an <af_bandFilter> parameter, in which case the curvature
%	vector is hardlimited to [-<af_bandFilter>, ... <af_bandFilter>].
%
% PRECONDITIONS
%       o <astr_mris> and <astr_curv> should be valid filenamesi.
%       o FreeSurfer environment.
%
% POSTCONDITIONS
%       o Figure is generated.
%       o handle to figure is returned.
%
% HISTORY
% 26 August 2009
% o Initial design and coding.
% 
% 02 June 2011
% o Fixed colormap handling for cases where mapping toolbox is not
%   available.
%

% ---------------------------------------------------------

%%%%%%%%%%%%%% 
%%% Nested functions :START
%%%%%%%%%%%%%% 
	function error_exit(	str_action, str_msg, str_ret)
		fprintf(1, '\tFATAL:\n');
		fprintf(1, '\tSorry, some error has occurred.\n');
		fprintf(1, '\tWhile %s,\n', str_action);
		fprintf(1, '\t%s\n', str_msg);
		error(str_ret);
	end

%%%%%%%%%%%%%% 
%%% Nested functions :END
%%%%%%%%%%%%%% 


sys_print('mris_display: START\n');
 
az              = 270;
el              = 0;

b_bandFilter    = 0;
av_bandFilter   = [-1 1];
b_invCurv       = 0;
b_colorMap      = 0;
str_colorMap    = 'Jet';

% Parse optional arguments
if length(varargin) >= 1, str_title = varargin{1};      end
if length(varargin) >= 2, az = varargin{2};             end
if length(varargin) >= 3, el = varargin{3};             end
if length(varargin) >= 4
    b_bandFilter        = 1;
    a_bandFilter        = varargin{4};
end
if length(varargin) >= 5, b_invCurv = varargin{5};      end
if length(varargin) >= 6
    b_colorMap          = 1;      
    str_colorMap        = varargin{6};
end

% Read curvature file
cprintsn('Reading curvature file', astr_curv);
[v_curv, fnum] = read_curv(astr_curv);
cprintdn('Number of curv elements', numel(v_curv));
if b_invCurv
    cprintdn('Invert curvature data sign', b_invCurv);
    v_curv = v_curv * -1;
end

% Read surface
cprintsn('Reading mris file', astr_mris);
[v_vertices, v_faces] = read_surf(astr_mris);
v_vertSize      = size(v_vertices);
v_faceSize      = size(v_faces);
str_vertSize    = sprintf('%d x %d', v_vertSize(1), v_vertSize(2));
str_faceSize    = sprintf('%d x %d', v_faceSize(1), v_faceSize(2));
cprintsn('Size of vert struct', str_vertSize);
cprintsn('Size of face struct', str_faceSize);

if numel(v_curv) ~= v_vertSize(1)
    error_exit( 'reading inputs',        ...
                'mismatch between curvature size and surf vertices', ...
                '1');
end

v_vertices      = single(v_vertices);
v_faces         = int32(v_faces+1);  % for matlab compliance

av_filtered     = [0 0];
if b_bandFilter
    v_bandFilter        = [min(v_curv) max(v_curv)];
    if isfloat(a_bandFilter)
        v_bandFilter    = [-a_bandFilter a_bandFilter];
    end
    if isvector(a_bandFilter)
        v_bandFilter    = a_bandFilter;
    end
    if ischar(a_bandFilter)
        if strcmp(a_bandFilter, 'gz')
            v_bandFilter = [min(v_curv(find(v_curv>0))) max(v_curv)];
        end
    end
    lowerCount  = numel(find(v_curv<v_bandFilter(1)));
    upperCount  = numel(find(v_curv>v_bandFilter(2)));
    av_filtered = [lowerCount upperCount];
    v_curv      = filter_bandPass(v_curv,                       ...
                                 v_bandFilter(1), v_bandFilter(2), 1);
end

% Display:
hf              = figure;
hp              = patch('vertices',     v_vertices,             ...
                        'faces',        v_faces(:,[1 3 2]),     ...
                        'facevertexcdata',      v_curv,         ...
                        'edgecolor',    'none',                 ...
                        'facecolor',    'interp'); 
axis equal; 
grid;
try
    demcmap(v_curv);
catch ME
    colormap(str_colorMap);
end

if b_colorMap
    colormap(str_colorMap);
end

if length(str_title)
    title(str_title);
else
    title(sprintf('%s: %s', astr_mris, astr_curv));
end
colorbar;
view(az, el);

sys_print('mris_display: END\n');

av_curv = v_curv;
end
% ---------------------------------------------------------


