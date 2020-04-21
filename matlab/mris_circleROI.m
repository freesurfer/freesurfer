function [aCv_indicesROI, aCv_indicesNROI, aCv_tracker, aS_ct] =        ...
                                  mris_circleROI(                       ...
                                        astr_hemi,                      ...
                                        astr_mrisFile,                  ...
                                        av_center,                      ...
                                        af_radius,                      ...
                                        varargin)
%
% NAME
%  function [aCv_indicesROI, aCv_indicesNROI, aCv_tracker, aS_ct] =     ...
%                                    mris_circleROI(                    ...
%                                          astr_hemi,                   ...  
%                                          astr_mrisFile,               ...
%                                          av_center,                   ...
%                                          af_radius                    ...
%                                          [, astr_annotFile            ...
%                                          [, astr_labelStem,           ...
%                                          [, ab_progress]]])
%
%
%
% ARGUMENTS
%       
%       INPUT
%       astr_hemi       string          hemisphere to process
%       astr_mrisFile   string          filename of surface file to load
%                                       + (w/o hemisphere string)
%       av_center       vector          list of vertex indices about which
%                                       to determine ROIs.
%       af_radius       float           radius length
%
%       OPTIONAL
%       astr_annotFile  string          annotation file name
%                                       + (w/o hemisphere string)
%       astr_labelStem  string          stem to use for the label files
%       ab_progress     bool            if true, show running progress
%
%       OUTPUT
%       aCv_indicesROI  cell            indices within the ROI 
%       aCv_indicesNROI cell            indices outside the ROI
%       aCv_tracker     cell            indices across whole vertex space:
%                                       + 1 indicates corresponding vertex
%                                       +   index is in the ROI;
%                                       + 0 indicates corresponding vertex
%                                       +   index in in the NROI       
%       aS_ct           struct          color table structure that was
%                                       written to file. If no annotation
%                                       was created, this is empty.
%
% DESCRIPTION
%
%       'mris_circleROI' accepts a list of vertex points and
%       generates circular ROIs about each point. The vertex indices
%       are returned as cell arrays with ROIs and non-ROIs.
%       
%       If an optional annotation output filename is provided,
%       the ROIs are added to a FreeSurfer annotation file suitable
%       for uploading onto surfaces. If a color table structure is
%       passed, this is used for the annotation, otherwise a
%       color table is created (or attempted).
%       
% PRECONDITIONS
%       o <astr_mris> and <aS_colorTable> should be valid.
%       o FreeSurfer environment.
%
% POSTCONDITIONS
%       o Vertex indices are returned and optional annotation file
%         is created.
%         
% SEE ALSO
%       o read_annotation / write_annotation for a description of the
%         color table format.
%
% HISTORY
% 06 June 2011
% o Initial design and coding.
%   WARNING! Passing external colorTable is not fully supported yet!
%

% ---------------------------------------------------------

sys_printf('mris_circleROI: START\n');
 
aS_ct           = struct;
b_annotate      = 0;
b_colorTable    = 0;
b_progress      = 0;
str_labelStem   = '';
% Parse optional arguments
if length(varargin) >=1 && length(varargin{1})
    b_annotate          = 1;
    str_annotationFile  = [astr_hemi '.' varargin{1}];
end
if length(varargin) >=2 && length(varargin{2})
    str_labelStem       = varargin{2};
end
if length(varargin) >=3 && length(varargin{3})
    b_progress          = varargin{3};
end

astr_mrisFile   = [astr_hemi '.' astr_mrisFile];

% Read surface
colprintf('40;40', 'Reading mris file', '[ %s ]\n', astr_mrisFile);
[v_vertices, v_faces] = read_surf(astr_mrisFile);
v_vertSize      = size(v_vertices);
v_faceSize      = size(v_faces);
str_vertSize    = sprintf('%d x %d', v_vertSize(1), v_vertSize(2));
str_faceSize    = sprintf('%d x %d', v_faceSize(1), v_faceSize(2));
colprintf('40;40', 'Size of vert struct', '[ %s ]\n', str_vertSize);
colprintf('40;40', 'Size of face struct', '[ %s ]\n', str_faceSize);

numROIcenters   = numel(av_center);

if numROIcenters
    aCv_indicesROI      = cell(1, numROIcenters);
    aCv_indicesNROI     = cell(1, numROIcenters);
    aCv_tracker         = cell(1, numROIcenters);
    for ROI=1:numROIcenters
        aCv_tracker{ROI} = zeros(1, v_vertSize(1));
    end
    if ~b_progress
        colprintf('40;40', 'Processing vertices...', '');
    end
    for vi = 1:length(v_vertices)
        if b_progress
            colprintf('40;40', 'Vertex process count', '[ %6d (%3.2f%s) ]\r', ... 
                    vi, vi/v_vertSize(1)*100, '%');
        end
        for ROI=1:numROIcenters
            v_ROIcenter = v_vertices(av_center(ROI), :);
            v_vertices(vi, :);
            f_dist      = norm(v_vertices(vi, :) - v_ROIcenter);
            if (f_dist < af_radius)
                aCv_tracker{ROI}(vi)    = 1;
                aCv_indicesROI{ROI}     = [ aCv_indicesROI{ROI} vi];
            else
                aCv_tracker{ROI}(vi)    = 0;
                aCv_indicesNROI{ROI}    = [ aCv_indicesNROI{ROI} vi];
            end
        end
    end
    if ~b_progress
        colprintf('40;40', '', '[ done ]');
    end
end

fprintf('\n');
if b_annotate
    b_useUnknownLabel   = 0;
    if ~b_colorTable
        b_useUnknownLabel = 1;
        S_ct    = mris_colorTableMake(numROIcenters, b_useUnknownLabel);
    end
    v_annotV    = [0:v_vertSize(1)-1];
    v_elementID = zeros(1, v_vertSize(1));
    for ROI=1:numROIcenters
        label   = S_ct.table(ROI+b_useUnknownLabel, 5);
        colprintf('40;40', 'Annotating ROI:lookup     ', '');
        v_elementID(aCv_indicesROI{ROI}) = label;
        colprintf('40;40', '', '[ %d:%d ]\n', ROI, label);
        str_hemi        = basename(strtok(str_annotationFile, '-'));
        str_labelFile   = sprintf('%s/%s_%s%s.label',                   ...
                                dirname(str_annotationFile),            ...
                                astr_hemi,                              ...
                                str_labelStem,                          ...
                                S_ct.struct_names{ROI+b_useUnknownLabel});
        colprintf('40;40', 'Writing label file', '');
        v_verticesLabel = v_vertices(aCv_indicesROI{ROI},:);
        v_labelVals     = zeros(length(aCv_indicesROI{ROI}),1);
        write_label(aCv_indicesROI{ROI}'-1, v_verticesLabel, v_labelVals,   ...
                    str_labelFile);
        colprintf('40;40', '', '[ %s ]\n', str_labelFile);
        colprintf('40;40', 'Label elements', '[ %d ]\n', numel(v_labelVals));
    end
    colprintf('40;40', 'Writing annotation file', '');
    write_annotation(str_annotationFile, v_annotV, v_elementID, S_ct);
    colprintf('40;40', '', '[ %s ]\n', str_annotationFile);
    aS_ct       = S_ct;
end

sys_printf('mris_circleROI: END\n');

end
% ---------------------------------------------------------


