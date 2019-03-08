function [D,P] = dijk(A,s,t)
%DIJK Shortest paths from nodes 's' to nodes 't' using Dijkstra algorithm.
% [D,P] = dijk(A,s,t)
%     A = n x n node-node weighted adjacency matrix of arc lengths
%         (Note: A(i,j) = 0   => Arc (i,j) does not exist;
%                A(i,j) = NaN => Arc (i,j) exists with 0 weight)
%     s = FROM node indices
%       = [] (default), paths from all nodes
%     t = TO node indices
%       = [] (default), paths to all nodes
%     D = |s| x |t| matrix of shortest path distances from 's' to 't'
%       = [D(i,j)], where D(i,j) = distance from node 'i' to node 'j'
%                                = Inf, if no path from 'i' to 'j'
%     P = |s| x n matrix of predecessor indices, where P(i,j) is the
%         index of the predecessor to node 'j' on the path from 's(i)' to
%         'j',where P(i,i) = 0 and P(i,j) = NaN is 'j' not on path to 
%         's(i)' (use PRED2PATH to convert P to paths)
%       = path from 's' to 't', if |s| = |t| = 1
%
% Example:
% A = [0 1 3 0
%      0 0 0 2
%      0 0 0 4
%      0 0 0 0]
% [d,p] = dijk(A,1,4)   % (Single path) d =  3
%                       %               p =  1   2   4
%
% [D,P] = dijk(A)       % (All paths)   D =  0   1   3   3
%                       %                  Inf   0 Inf   2
%                       %                  Inf Inf   0   4
%                       %                  Inf Inf Inf   0
%                       %               P =  0   1   1   2
%                       %                  NaN   0 NaN   2
%                       %                  NaN NaN   0   3
%                       %                  NaN NaN NaN   0
% p = pred2path(P,1,4)  %               p =  1   2   4
%
% (If A is a triangular matrix, then computationally intensive node
%  selection step not needed since graph is acyclic (triangularity is a 
%  sufficient, but not a necessary, condition for a graph to be acyclic)
%  and A can have non-negative elements)
%
% (If |s| >> |t|, then DIJK is faster if DIJK(A',t,s) used, where D is now
%  transposed and P now represents successor indices)
%
% (Based on Fig. 4.6 in Ahuja, Magnanti, and Orlin, Network Flows,
%  Prentice-Hall, 1993, p. 109.)

% Copyright (c) 1994-2006 by Michael G. Kay
% Matlog Version 9 13-Jan-2006 (http://www.ie.ncsu.edu/kay/matlog)

% Input Error Checking ****************************************************
error(nargchk(1,3,nargin))

[n,cA] = size(A);

if nargin < 2 || isempty(s), s = (1:n)'; else s = s(:); end
if nargin < 3 || isempty(t), t = (1:n)'; else t = t(:); end

if ~any(any(tril(A) ~= 0))       % A is upper triangular
   isAcyclic = 1;
elseif ~any(any(triu(A) ~= 0))   % A is lower triangular
   isAcyclic = 2;
else                             % Graph may not be acyclic
   isAcyclic = 0;
end

if n ~= cA
   error('A must be a square matrix');
elseif ~isAcyclic && any(any(A < 0))
   error('A must be non-negative');
elseif any(s < 1 | s > n)
   error(['"s" must be an integer between 1 and ',num2str(n)]);
elseif any(t < 1 | t > n)
   error(['"t" must be an integer between 1 and ',num2str(n)]);
end
% End (Input Error Checking) **********************************************

A = A';    % Use transpose to speed-up FIND for sparse A

D = zeros(length(s),length(t));
if nargout > 1, P = NaN*ones(length(s),n); end

for i = 1:length(s)
   j = s(i);
   
   Di = Inf*ones(n,1); Di(j) = 0;
   
   isLab = false(length(t),1);
   if isAcyclic ==  1
      nLab = j - 1;
   elseif isAcyclic == 2
      nLab = n - j;
   else
      nLab = 0;
      UnLab = 1:n;
      isUnLab = true(n,1);
   end
   
   if nargout > 1, P(i,s(i)) = 0; end  %Change from NaN to indicate no pred
   
   while nLab < n && ~all(isLab)
      if isAcyclic
         Dj = Di(j);
      else	% Node selection
         [Dj,jj] = min(Di(isUnLab));
         j = UnLab(jj);
         UnLab(jj) = [];
         isUnLab(j) = 0;
      end
      
      nLab = nLab + 1;
      if length(t) < n, isLab = isLab | (j == t); end
      
      [jA,kA,Aj] = find(A(:,j));
      Aj(isnan(Aj)) = 0;
            
      if isempty(Aj), Dk = Inf; else Dk = Dj + Aj; end
      
      if nargout > 1, P(i,jA(Dk < Di(jA))) = j; end
      Di(jA) = min(Di(jA),Dk);
      
      if isAcyclic == 1       % Increment node index for upper triangular A
         j = j + 1;
      elseif isAcyclic == 2   % Decrement node index for lower triangular A
         j = j - 1;
      end
   end
   D(i,:) = Di(t)';
end

if nargout > 1 && length(s) == 1 && length(t) == 1
   P = pred2path(P,s,t);
end
