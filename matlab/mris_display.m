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
%                               a_mapOperation,         ...
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
%                                       surface and curvature file names.
%                                       If title string starts with "save:"
%                                       then the graph is also saved to
%                                       filesystem using the title string
%                                       as a filestem.
%       a_az            float           azimuth of viewpoint
%       a_el            float           elevation of viewpoint
%       a_bandFilter    variable        apply a bandpass filter over data:
%                                           if vector, defines range
%                                           if single float, defines
%                                              symmetrical pulse [-f, +f]
%                                           if string == 'gz'
%                                              define range as [>0, max] 
%                                       between vector range.
%       a_mapOperation  string          apply an operation on the curvature
%                                       data post filtering.
%                                           'none' | '':        no operation 
%                                           'inv'      :        invert sign
%                                           'signed'   :        shift values
%                                                               symmetrically
%                                                               about 0.
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
b_mapOperation  = 0;
b_colorMap      = 0;
str_colorMap    = 'Jet';
str_title       = '';

% Parse optional arguments
if length(varargin) >= 1, str_title = varargin{1};      end
if length(varargin) >= 2, az = varargin{2};             end
if length(varargin) >= 3, el = varargin{3};             end
if length(varargin) >= 4
    b_bandFilter        = 1;
    a_bandFilter        = varargin{4};
end
if length(varargin) >= 5,
    b_mapOperation      = 1;
    a_mapOperation      = varargin{5};      
end
if length(varargin) >= 6
    b_colorMap          = 1;      
    str_colorMap        = varargin{6};
end

% Parse title string for 'save:' directive
b_save			= 0;
c_title			= regexp(str_title, 'save:', 'split');
if numel(c_title) == 2
    b_save		= 1;
    str_title		= c_title{2};
end

% Read curvature file
colprintf('40;40', 'Reading curvature file', '[ %s ]\n', astr_curv);
[v_curv, fnum] = read_curv(astr_curv);
colprintf('40;40', 'Number of curv elements', '[ %d ]\n', numel(v_curv));

% Read surface
colprintf('40;40', 'Reading mris file', '[ %s ]\n', astr_mris);
[v_vertices, v_faces] = read_surf(astr_mris);
v_vertSize      = size(v_vertices);
v_faceSize      = size(v_faces);
str_vertSize    = sprintf('%d x %d', v_vertSize(1), v_vertSize(2));
str_faceSize    = sprintf('%d x %d', v_faceSize(1), v_faceSize(2));
colprintf('40;40', 'Size of vert struct', '[ %s ]\n', str_vertSize);
colprintf('40;40', 'Size of face struct', '[ %s ]\n', str_faceSize);

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
    if isvector(a_bandFilter) & numel(a_bandFilter) == 2
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

if b_mapOperation
    if strcmp(a_mapOperation, 'inv')    
        colprintf('40;40', 'Inverting curvature data sign', '[ ok ]\n');
        v_curv = v_curv * -1;
    end
    if strcmp(a_mapOperation, 'signed')
        colprintf('40;40', 'Shifting curvs evenly about zero', '[ ok ]\n');
        f_range = max(v_curv) - min(v_curv);
        v_curv = v_curv - min(v_curv) - f_range/2;
    end
end 


% Display:
hf              = figure;
hp              = patch('vertices',             v_vertices,             ...
                        'faces',                v_faces(:,[1 3 2]),     ...
                        'facevertexcdata',      v_curv,                 ...
                        'edgecolor',            'none',                 ...
                        'facecolor',            'interp'); 
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

if b_save
    str_epsFile	= [ str_title '.eps' ];
    str_pngFile	= [ str_title '.png' ];
    colprintf('40;40', 'Saving eps snapshot', '');
    print('-depsc2', str_epsFile);
    colprintf('40;40', '', '[ %s ]\n', str_epsFile);
    colprintf('40;40', 'Saving png snapshot', '');
    print('-dpng', str_pngFile);
    colprintf('40;40', '', '[ %s ]\n', str_pngFile);
end

sys_print('mris_display: END\n');

av_curv = v_curv;
end
% ---------------------------------------------------------


