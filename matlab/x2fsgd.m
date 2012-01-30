function x2fsgd(X,fp)
% x2fsgd(X,fp)
%
% Converts a design matrix to an FSGD file assuming that each
% column codes for a separate class
%
% Example:
% fp = fopen('my.fsgd','w');
% X = zeros(3,2); X(1,1) = 1;X(2:3,2) = 1;
% x2fsgd(X,fp);
% fclose(fp);
%
% $Id: x2fsgd.m,v 1.1 2012/01/30 17:08:50 greve Exp $

if(nargin ~= 2)
  fprintf('x2fsgd(X,fp)\n');
end

nclasses = size(X,2);

fprintf(fp,'GroupDescriptorFile 1\n');
for nthclass = 1:nclasses
  fprintf(fp,'Class%02d\n',nthclass);
end

for nthclass = 1:nclasses
  ind = find(X(:,nthclass));
  for nthsubject = 1:length(ind)
    fprintf(fp,'Input Subject%04d Class%02d\n',ind(nthsubject),nthclass);
  end
end

return;

