function [av_curvDiff, v_verticesA, v_verticesB] = mris_XYZdiff(astr_mrisA, astr_mrisB, astr_curvDiff, ...
                                 varargin)
%
% NAME
%
%       [av_curvDiff] = mris_XYZdiff(
%                               astr_mrisA,             ...
%                               astr_mrisB,             ...
%                               astr_curvDiff
%                               <ab_displayDistorion>
%                         )
%
%
% ARGUMENTS
%       
%       INPUT
%       astr_mrisA      string          filename of input surface A
%       astr_mrisA      string          filename of input surface B
%       astr_curvDiff   string          filename of output curv file containing
%                                               difference
%
%       OPTIONAL
%       ab_displayDistorion bool        if true, display the distortion, as
%                                               projected onto <astr_mrisA>
%
%       OUTPUT
%       av_curvDiff     vector          the per-vertex difference
%
% DESCRIPTION
%
%       'mris_XYZdiff' is used to quantify the metric distortion between two
%       representations of a surface structure. Typically <astr_mrisA> is the
%       original surface, and <astr_mrisB> is the original surface subjected
%       to some distorition.
%       
%       Simply stated, the Euclidean distance between each corresponding vertex
%       in each surface is computed and stored as a FreeSurfer 'curvature' type
%       file in <astr_curvDiff>.
%       
% PRECONDITIONS
%       o <astr_mrisA> and <astr_mrisB> must have the same number of vertices
%         and faces.
%       o FreeSurfer environment.
%
% POSTCONDITIONS
%       o The metric distortion is stored in <astr_curvDiff>.
%       o The metric distortion is also (optionally) displayed.
%
% HISTORY
% 11 September 2009
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

        function [v_vertices, v_faces] = mris_read(astr_surf)
            cprintsn('Reading mris file', astr_surf);
            [v_vertices, v_faces] = read_surf(astr_surf);
            v_vertSize      = size(v_vertices);
            v_faceSize      = size(v_faces);
            str_vertSize    = sprintf('%d x %d', v_vertSize(1), v_vertSize(2));
            str_faceSize    = sprintf('%d x %d', v_faceSize(1), v_faceSize(2));
            cprintsn('Size of vert struct', str_vertSize);
            cprintsn('Size of face struct', str_faceSize);
        end

        function [v_curv] = curv_write(astr_curvDiff, av_curvDiff, afnum)
            cprintsn('Writing curvature difference file', astr_curvDiff);
            v_vertSize      = size(av_curvDiff);
            str_vertSize    = sprintf('%d x %d', v_vertSize(1), v_vertSize(2));
            cprintsn('Size of curvDiff struct', str_vertSize);
            [v_curv]        = write_curv(astr_curvDiff, av_curvDiff, afnum);
        end


%%%%%%%%%%%%%% 
%%% Nested functions :END
%%%%%%%%%%%%%% 


sys_print('mris_XYZdiff: START\n');
 
az              = 270;
el              = 0;

b_display       = 0;

% Parse optional arguments
if length(varargin) >= 1, b_display = varargin{1};      end

% Read surfaces
[v_verticesA, v_facesA] = mris_read(astr_mrisA);
[v_verticesB, v_facesB] = mris_read(astr_mrisB);


if numel(v_verticesA) ~= numel(v_verticesB)
    error_exit( 'reading inputs',        ...
                'mismatch between surf vertices', ...
                '1');
end
if numel(v_facesA) ~= numel(v_facesB)
    error_exit( 'reading inputs',        ...
                'mismatch between surf faces', ...
                '1');
end

v_fnum          = size(v_facesA);
v_vertices      = size(v_verticesA);
fnum            = v_fnum(1);
v_curvDiff      = zeros(v_vertices(1), 1);
for v=1:v_vertices(1)
    v_curvDiff(v) = vector_distance(v_verticesA(v,:), v_verticesB(v,:));
end

curv_write(astr_curvDiff, v_curvDiff, fnum);
av_curvDiff     = v_curvDiff;
if b_display, mris_display(astr_mrisA, astr_curvDiff); end


sys_print('mris_XYZdiff: END\n');

end
% ---------------------------------------------------------


