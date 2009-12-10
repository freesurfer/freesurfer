function [oxy deoxy] = fast_nirs2oxy(rawnirs)
% [oxy deoxy] = fast_nirs2oxy(rawnirs)
%
% rawnirs is ntp by nSDs*2 arranged with matching wavelengths in the
% nth and nth+nSD locations. Assumes 690 and 830. Based on code
% from D Boas.
%
% $Id: fast_nirs2oxy.m,v 1.1 2009/12/10 22:32:14 greve Exp $

if(nargin ~= 1)
  fprintf('[oxy deoxy] = fast_nirs2oxy(rawnirs)\n');
  return;
end

[ntp nSD2] = size(rawnirs);
nSD = nSD2/2;

n1 = 1:nSD;
n2 = nSD+1:2*nSD;

w1 = rawnirs(:,n1);
w2 = rawnirs(:,n2);

% HB0,690 HBR,690
% HB0,830 HBR,830

e = [ [    0.635    4.723];
    [    2.243    1.596] ];

%einv = inv( e'*e )*e';
einv = inv(e);

oxy = zeros(ntp,nSD);
deoxy = zeros(ntp,nSD);
for nthSD = 1:nSD
  w1nth = w1(:,nthSD);
  w1nth = w1nth/mean(w1nth);
  w1nth = log(w1nth);
  w2nth = w2(:,nthSD);
  w2nth = w2nth/mean(w2nth);
  w2nth = log(w2nth);
  w = [w1nth w2nth];
  od = w*einv';
  oxy(:,nthSD)   = od(:,1);
  deoxy(:,nthSD) = od(:,2);
end



return;

% d = rawnirs;
% mean1 = mean(d,1);
% [M,N]=size(d);
% dm1 = d./(ones(M,1)*mean1);  %set dm1 to the normalized data list
% dd = log( dm1 );            %dod1 is the log of the normalized data
% %ADD YOUR FAVORITE BANDPASS FILTER HERE

% lst = find( ml(:,4)==1 );
% for idx=1:length(lst)
%     idx1 = lst(idx);
%     idx2 = find( ml(:,4)>1 & ml(:,1)==ml(idx1,1) & ml(:,2)==ml(idx1,2) );
%     concs(:,:,idx) = ( einv * dd(:,[idx1 idx2])' )';
% end

