function err = fast_svbslice(y,stem,sliceno,bext,bhdrstr)
% err = fast_svbslice(y,stem,sliceno,<bext>)
%
% size(y) = [rows cols frames]
% sliceno is zero-based
% if bext is not specificed or is null, it is set to bfloat
% bhdrstr is printed to stem.bhdr

err = 1;

if(nargin < 2 & nargin > 4)
  fprintf('err = fast_svbslice(y,stem,sliceno,<bext>,<bhdrstr>)\n');
  return;
end

if(exist('bext') ~= 1) bext = ''; end
if(isempty(bext)) bext = 'bfloat'; end

if(exist('bhdrstr') ~= 1) bhdrstr = ''; end

if(strcmp(bext,'bfloat') == 0 & strcmp(bext,'bshort') == 0)
  fprintf('ERROR: bext = %s, must be bfloat or bshort\n',bext);
  return;
end
  
if(sliceno < 0)
  fprintf('ERROR: cannot have negative slice number\n');
  return;
end
  

fname = sprintf('%s_%03d.%s',stem,sliceno,bext);
fmri_svbfile(y,fname);

if(~isempty(bhdrstr))
  % Should make sure nframes is correct %
  bhdrfile = sprintf('%s.bhdr',stem);
  fid = fopen(bhdrfile,'w');
  fprintf(fid,'%s',bhdrstr);
  fclose(fid);
end


err = 0;

return;







