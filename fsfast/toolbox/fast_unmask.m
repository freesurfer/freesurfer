function vol = fast_unmask(data,imask,volsize)
% vol = fast_unmask(data,imask,volsize)
%
% data is nf-by-nmask
% imask has nmask elements
% Creates vol which is ns-nr-nc-nf with zeros
%  in places without mask. 
% If imask is empty, converts data into a volume.
%

vol = [];

if(nargin ~= 3)
  fprintf('vol = fast_unmask(data,imask,volsize)\n');
  return;
end

if(length(volsize) < 3)
  fprintf('ERROR: volsize must have at least 3 dimensions \n');
  return;
end
volsize = volsize(1:3);
nv = prod(volsize);
[nf ndata] = size(data);

if(isempty(imask))
  if(ndata ~= nv)
    fprintf('ERROR: ndata ~= nv\n');
    return;
  end
  vol = reshape(data', [volsize nf]); %'
  return;
end

if(max(imask) > nv)
  fprintf('ERROR: max(imask) > nv \n');
  return;
end

if(size(data,2)==1) data = data'; end %'

nmask = length(imask);
if(ndata ~= nmask)
  fprintf('ERROR: ndata ~= nmask \n');
  return;
end

vol = zeros([nf nv]);
vol(:,imask) = data;
vol = reshape(vol', [volsize nf]); %'

return;
