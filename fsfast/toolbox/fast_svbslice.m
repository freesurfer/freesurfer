function err = fast_svbslice(y,stem,sliceno,bext)
% err = fast_svbslice(y,stem,sliceno,<bext>)
%
% size(y) = [rows cols frames]
% sliceno is zero-based
% if bext is not specificed, it is set to bfloat

err = 1;

if(nargin ~= 2 & nargin ~= 3)
  fprintf('err = fast_svbslice(y,stem,sliceno,<bext>)\n');
  return;
end

if(exist('bext') ~= 1) bext = 'bfloat'; end

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

err = 0;

return;







