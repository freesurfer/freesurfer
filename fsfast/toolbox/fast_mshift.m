function m2 = fast_mshift(m1, shift, wrap)
% m2 = fast_mshift(m1, shift, wrap)
% 
% Shifts m1 by shift. Uses wrap-around if wrap=1. The number of 
% elements in shift must equal to the dimension of m1. This is similar
% the matlab native wshift() function except that wshift can only handle
% two dimensions and forces wrap-around.
%
% Positives shift down or to the right.
%
% Examples:
%
%  x = [1 2 3 4 5];
%  x2 = fast_mshift(x,[0 +2])   --> x2 = [0 0 1 2 3]
%  x2 = fast_mshift(x,[0 +2],1) --> x2 = [4 5 1 2 3]
%  x2 = fast_mshift(x,[0 -2])   --> x2 = [3 4 5 0 0]
%  x2 = fast_mshift(x,[0 -2],1) --> x2 = [3 4 5 1 2]
% 
%  x = [1 2 3; 
%       4 5 6];
%  x2 = fast_mshift(x,[0 +1])
%     0     1     2
%     0     4     5
%  x2 = fast_mshift(x,[0 +1],1)
%     3     1     2
%     6     4     5
%  x2 = fast_mshift(x,[-1 +1],1)
%     6     4     5
%     3     1     2


%
% fast_mshift.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

m2 = [];

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: m2 = fast_mshift(m1, shift, <wrap>)';
  qoe(msg);error(msg);
end

if(nargin ~= 3) wrap = 0; end
     
szm1 = size(m1);
m1dim = length(szm1);
lenshift = length(shift);

if(m1dim ~= lenshift)
  msg = sprintf('m1 dim (%d) != length of shift vector (%d)',m1dim,lenshift);
  qoe(msg);error(msg);
end

if(~wrap)
  % If not going to wrap, then return zeros if any
  % shift value is greater than the length of 
  % its dimension.
  ind = find(shift >= szm1);
  if(~isempty(ind)) 
    m2 = zeros(size(m1));
    return;
  end
end

% Make all shifts betwee -szm1  to szm1 %
for n = 1:m1dim
  shift(n) = sign(shift(n))*mod(abs(shift(n)),szm1(n));
end

m2 = m1;

for d = 1:m1dim
  szm2 = size(m2);
  nd = szm2(1);
  ndshift = shift(d);

  if(ndshift ~= 0 | nd == 1)
    if(ndshift < 0) ndshift = ndshift + nd; end

    m2 = reshape2d(m2);
    tmp = m2;

    if(shift(d) > 0)
      m2(ndshift+1:nd,:) =  tmp(1:nd-ndshift,:);
      if(wrap) 
        m2(1:ndshift,:) = tmp(nd-ndshift+1:nd,:); 
      else
        m2(1:ndshift,:) = zeros(size(m2(1:ndshift,:)));
      end
    end

    if(shift(d) < 0)
      m2(1:ndshift,:) = tmp(nd-ndshift+1:nd,:); 
      if(wrap) 
        m2(ndshift+1:nd,:) =  tmp(1:nd-ndshift,:); 
      else
        m2(ndshift+1:nd,:) =  zeros(size(m2(ndshift+1:nd,:)));
      end
    end

    m2 = reshape(m2, szm2);
  end

  m2 = shiftdim(m2,1);
end


return;


