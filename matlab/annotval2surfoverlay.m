function surfoverlay = annotval2surfoverlay(annotvals,annotnames,annotfile)
% surfoverlay = annotval2surfoverlay(annotvals,annotnames,annotfile)
%
% Example:
% I have three annotations
%   annotnames = strvcat('superiortemporal','insula','postcentral');
% with matchine values  
%   annotvals = [1.1 2.2 3.7];
% I want to create a surface overlay with all the vertices in a
% given annotation having its corresponding value (ie, all vertices
% in superior temporal gyrus being a value of 1.1, etc)
%
% Run these commands in matlab
% annotnames = strvcat('superiortemporal','insula','postcentral');
% annotvals = [1.1 2.2 3.7];
% annotfile = '~/subjects/fsaverage/label/lh.aparc.annot';
% surfoverlay = annotval2surfoverlay(annotvals,annotnames,annotfile);
% clear mri
% mri.vol = surfoverlay;
% MRIwrite(mri,'vals.mgh');
%
% Run this command from a linux shell:
%   tksurferfv fsaverage lh inflated -aparc -ov vals.mgh -fminmax .01 1
%

surfoverlay = [];

% Check that the number of values equals the number of annotation names
nannots = length(annotvals);
if(nannots ~= size(annotnames,1))
  fprintf('ERROR: size mismatch between annotvals and annotnames\n');
  return;
end

% Read in the annotation file
[vertices, labelcode, colortable] = read_annotation(annotfile);
if(isempty(vertices))
  fprintf('ERROR: reading %s\n',annotfile);
  reutrn;
end
nvertices = length(vertices);

surfoverlay = zeros(1,nvertices);
for nthannot = 1:nannots
  annotname = deblank(annotnames(nthannot,:));
  % Get the index into the color table
  indctab = strmatch(annotname,char(colortable.struct_names),'exact');
  if(isempty(indctab))
    fprintf('ERROR: count not find annotname %s in annot\n',annotname);
    surfoverlay = [];
    reutrn;
  end
  % Get the label code
  annotcode = colortable.table(indctab,5);
  % Find matching indices
  ind = find(labelcode == annotcode);
  % Set all indices to that value
  surfoverlay(ind) = annotvals(nthannot);
end

return

