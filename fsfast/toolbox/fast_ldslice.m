function slice = fast_ldslice(volid,sliceno)
% slice = fast_ldslice(volid,sliceno)

slice = [];

if(nargin ~= 2)
  msg = 'USAGE: slice = fast_ldslice(volid,sliceno)'
  qoe(msg) ; error(msg);
end

fmt = fast_getvolformat(voldid);
if(isempty(fmt))
  msg = sprintf('Could not determine format of %s',volid);
  qoe(msg) ; error(msg);
end

switch(fmt)

end



return;
