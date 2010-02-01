function [hf, hp, av_filtered] = mris_display(astr_mris, astr_curv, varargin)
%
% NAME
%       [hf, hp, av_filtered]   =                       ...
%       mris_display(     	astr_mris,              ... 
%                               astr_curv               ...
%                               <, a_az,                ... 
%                               a_el,                   ...
%                               af_bandFilter,          ...
%                               ab_invMap>
%                         )
% $Id: mris_display.m,v 1.3 2010/02/01 21:16:13 rudolph Exp $
%
%
% ARGUMENTS
%       
%       INPUT
%       astr_mris       string          filename of surface file to load
%       astr_curv       string          filename of curvature file to load
%
%       OPTIONAL
%       a_az            float           azimuth of viewpoint
%       a_el            float           elevation of viewpoint
%       af_bandFilter   float           apply a symmetrical band pass filt
%                                       between +- <af_bandFilter>
%       ab_invMap       bool            if true, invert the sign of the
%                                       curvature data -- useful if 'neg'
%                                       values are on gyri and 'pos' values
%                                       in sulci.
%
%       OUTPUT
%       hf              handle          handle to generated figure
%       hp              handle          handle to patch
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

	function vprintf(level, str_msg)
	    if verbosity >= level
		fprintf(1, str_msg);
	    end
	end

%%%%%%%%%%%%%% 
%%% Nested functions :END
%%%%%%%%%%%%%% 


sys_printf('mris_display: START\n');
 
az              = 270;
el              = 0;

b_bandFilter    = 0;
af_bandFilter   = 1;
b_invCurv       = 0;

% Parse optional arguments
if length(varargin) >= 1, az = varargin{1};             end
if length(varargin) >= 2, el = varargin{2};             end
if length(varargin) >= 3
    b_bandFilter        = 1;
    af_bandFilter       = varargin{3};
end
if length(varargin) >= 4, b_invCurv = varargin{4};      end

% Read curvature file
colprintf('40;40', 'Reading curvature file', astr_curv);
[v_curv, fnum] = read_curv(astr_curv);
colprintf('40;40', 'Number of curv elements', sprintf('%d', numel(v_curv)));
if b_invCurv
    colprintf('40;40', 'Invert curvature data sign', sprintf('%d', b_invCurv));
    v_curv = v_curv * -1;
end

% Read surface
colprintf('40;40', 'Reading mris file', astr_mris);
[v_vertices, v_faces] = read_surf(astr_mris);
v_vertSize      = size(v_vertices);
v_faceSize      = size(v_faces);
str_vertSize    = sprintf('%d x %d', v_vertSize(1), v_vertSize(2));
str_faceSize    = sprintf('%d x %d', v_faceSize(1), v_faceSize(2));
colprintf('40;40', 'Size of vert struct', str_vertSize);
colprintf('40;40', 'Size of face struct', str_faceSize);

if numel(v_curv) ~= v_vertSize(1)
    error_exit( 'reading inputs',        ...
                'mismatch between curvature size and surf vertices', ...
                '1');
end

v_vertices      = single(v_vertices);
v_faces         = int32(v_faces+1);  % for matlab compliance

av_filtered     = [0 0];
if b_bandFilter
    lowerCount  = numel(find(v_curv<-af_bandFilter));
    upperCount  = numel(find(v_curv>af_bandFilter));
    av_filtered = [lowerCount upperCount];
    v_curv      = filter_bandPass(v_curv,                       ...
                                -af_bandFilter, af_bandFilter, 1);
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
%colormap('winter');
demcmap(v_curv);
title(sprintf('%s: %s', astr_mris, astr_curv));
colorbar;
view(az, el);

sys_printf('mris_display: END\n');

end
% ---------------------------------------------------------


