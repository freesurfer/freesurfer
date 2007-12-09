function rte = pred2path(P,s,t)
%PRED2PATH Convert predecessor indices to shortest paths from 's' to 't'.
%   rte = pred2path(P,s,t)
%     P = |s| x n matrix of predecessor indices (from DIJK)
%     s = FROM node indices
%       = [] (default), paths from all nodes
%     t = TO node indices
%       = [] (default), paths to all nodes
%   rte = |s| x |t| cell array of paths (or routes) from 's' to 't', where
%         rte{i,j} = path from s(i) to t(j)
%                  = [], if no path exists from s(i) to t(j)
%
% (Used with output of DIJK)

% Copyright (c) 1994-2006 by Michael G. Kay
% Matlog Version 9 13-Jan-2006 (http://www.ie.ncsu.edu/kay/matlog)

% Input Error Checking ****************************************************
error(nargchk(1,3,nargin));

[rP,n] = size(P);

if nargin < 2 || isempty(s), s = (1:n)'; else s = s(:); end
if nargin < 3 || isempty(t), t = (1:n)'; else t = t(:); end

if any(P < 0 | P > n)
   error(['Elements of P must be integers between 1 and ',num2str(n)]);
elseif any(s < 1 | s > n)
   error(['"s" must be an integer between 1 and ',num2str(n)]);
elseif any(t < 1 | t > n)
   error(['"t" must be an integer between 1 and ',num2str(n)]);
end
% End (Input Error Checking) **********************************************

rte = cell(length(s),length(t));

[ans,idxs] = find(P==0);

for i = 1:length(s)
%    if rP == 1
%       si = 1;
%    else
%       si = s(i);
%       if si < 1 | si > rP
%          error('Invalid P matrix.')
%       end
%    end
   si = find(idxs == s(i));
   for j = 1:length(t)
      tj = t(j);
      if tj == s(i)
         r = tj;
      elseif P(si,tj) == 0
         r = [];
      else
         r = tj;
         while tj ~= 0
            if tj < 1 || tj > n
               error('Invalid element of P matrix found.')
            end
            r = [P(si,tj) r];
            tj = P(si,tj);
         end
         r(1) = [];
      end
      rte{i,j} = r;
   end
end

if length(s) == 1 && length(t) == 1
   rte = rte{:};
end

%rte = t;
while 0%t ~= s
   if t < 1 || t > n || round(t) ~= t
      error('Invalid "pred" element found prior to reaching "s"');
   end
   rte = [P(t) rte];
   t = P(t);
end

