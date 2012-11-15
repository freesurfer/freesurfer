function [costs,paths] = dijkstra(AorV,xyCorE,SID,FID,iswaitbar)
%DIJKSTRA Calculate Minimum Costs and Paths using Dijkstra's Algorithm
%
%   Inputs:
%     [AorV] Either A or V where
%         A   is a NxN adjacency matrix, where A(I,J) is nonzero (=1)
%               if and only if an edge connects point I to point J
%               NOTE: Works for both symmetric and asymmetric A
%         V   is a Nx2 (or Nx3) matrix of x,y,(z) coordinates
%     [xyCorE] Either xy or C or E (or E3) where
%         xy  is a Nx2 (or Nx3) matrix of x,y,(z) coordinates (equivalent to V)
%               NOTE: only valid with A as the first input
%         C   is a NxN cost (perhaps distance) matrix, where C(I,J) contains
%               the value of the cost to move from point I to point J
%               NOTE: only valid with A as the first input
%         E   is a Px2 matrix containing a list of edge connections
%               NOTE: only valid with V as the first input
%         E3  is a Px3 matrix containing a list of edge connections in the
%               first two columns and edge weights in the third column
%               NOTE: only valid with V as the first input
%     [SID] (optional) 1xL vector of starting points
%         if unspecified, the algorithm will calculate the minimal path from
%         all N points to the finish point(s) (automatically sets SID = 1:N)
%     [FID] (optional) 1xM vector of finish points
%         if unspecified, the algorithm will calculate the minimal path from
%         the starting point(s) to all N points (automatically sets FID = 1:N)
%     [iswaitbar] (optional) a scalar logical that initializes a waitbar if nonzero 
%
%   Outputs:
%     [costs] is an LxM matrix of minimum cost values for the minimal paths
%     [paths] is an LxM cell array containing the shortest path arrays
%
%   Revision Notes:
%     (4/29/09) Previously, this code ignored edges that have a cost of zero,
%     potentially producing an incorrect result when such a condition exists.
%     I have solved this issue by using NaNs in the table rather than a
%     sparse matrix of zeros. However, storing all of the NaNs requires more
%     memory than a sparse matrix. This may be an issue for massive data
%     sets, but only if there are one or more 0-cost edges, because a sparse
%     matrix is still used if all of the costs are positive.
%
%   Note:
%     If the inputs are [A,xy] or [V,E], the cost is assumed to be (and is
%       calculated as) the point-to-point Euclidean distance
%     If the inputs are [A,C] or [V,E3], the cost is obtained from either
%       the C matrix or from the edge weights in the 3rd column of E3
%
%   Example:
%       % Calculate the (all pairs) shortest distances and paths using [A,xy] inputs
%       n = 7; A = zeros(n); xy = 10*rand(n,2)
%       tri = delaunay(xy(:,1),xy(:,2));
%       I = tri(:); J = tri(:,[2 3 1]); J = J(:);
%       IJ = I + n*(J-1); A(IJ) = 1
%       [costs,paths] = dijkstra(A,xy)
%
%   Example:
%       % Calculate the (all pairs) shortest distances and paths using [A,C] inputs
%       n = 7; A = zeros(n); xy = 10*rand(n,2)
%       tri = delaunay(xy(:,1),xy(:,2));
%       I = tri(:); J = tri(:,[2 3 1]); J = J(:);
%       IJ = I + n*(J-1); A(IJ) = 1
%       a = (1:n); b = a(ones(n,1),:);
%       C = round(reshape(sqrt(sum((xy(b,:) - xy(b',:)).^2,2)),n,n))
%       [costs,paths] = dijkstra(A,C)
%
%   Example:
%       % Calculate the (all pairs) shortest distances and paths using [V,E] inputs
%       n = 7; V = 10*rand(n,2)
%       I = delaunay(V(:,1),V(:,2));
%       J = I(:,[2 3 1]); E = [I(:) J(:)]
%       [costs,paths] = dijkstra(V,E)
%
%   Example:
%       % Calculate the (all pairs) shortest distances and paths using [V,E3] inputs
%       n = 7; V = 10*rand(n,2)
%       I = delaunay(V(:,1),V(:,2));
%       J = I(:,[2 3 1]);
%       D = sqrt(sum((V(I(:),:) - V(J(:),:)).^2,2));
%       E3 = [I(:) J(:) D]
%       [costs,paths] = dijkstra(V,E3)
%
%   Example:
%       % Calculate the shortest distances and paths from the 3rd point to all the rest
%       n = 7; V = 10*rand(n,2)
%       I = delaunay(V(:,1),V(:,2));
%       J = I(:,[2 3 1]); E = [I(:) J(:)]
%       [costs,paths] = dijkstra(V,E,3)
%
%   Example:
%       % Calculate the shortest distances and paths from all points to the 2nd
%       n = 7; A = zeros(n); xy = 10*rand(n,2)
%       tri = delaunay(xy(:,1),xy(:,2));
%       I = tri(:); J = tri(:,[2 3 1]); J = J(:);
%       IJ = I + n*(J-1); A(IJ) = 1
%       [costs,paths] = dijkstra(A,xy,1:n,2)
%
%   Example:
%       % Calculate the shortest distance and path from points [1 3 4] to [2 3 5 7]
%       n = 7; V = 10*rand(n,2)
%       I = delaunay(V(:,1),V(:,2));
%       J = I(:,[2 3 1]); E = [I(:) J(:)]
%       [costs,paths] = dijkstra(V,E,[1 3 4],[2 3 5 7])
%
%   Example:
%       % Calculate the shortest distance and path between two points
%       n = 1000; A = zeros(n); xy = 10*rand(n,2);
%       tri = delaunay(xy(:,1),xy(:,2));
%       I = tri(:); J = tri(:,[2 3 1]); J = J(:);
%       D = sqrt(sum((xy(I,:)-xy(J,:)).^2,2));
%       I(D > 0.75,:) = []; J(D > 0.75,:) = [];
%       IJ = I + n*(J-1); A(IJ) = 1;
%       [cost,path] = dijkstra(A,xy,1,n)
%       gplot(A,xy,'k.:'); hold on;
%       plot(xy(path,1),xy(path,2),'ro-','LineWidth',2); hold off
%       title(sprintf('Distance from 1 to 1000 = %1.3f',cost))
%
% Web Resources:
%   <a href="http://en.wikipedia.org/wiki/Dijkstra%27s_algorithm">Dijkstra's Algorithm</a>
%   <a href="http://en.wikipedia.org/wiki/Graph_%28mathematics%29">Graphs</a>
%   <a href="http://en.wikipedia.org/wiki/Adjacency_matrix">Adjacency Matrix</a>
%
% See also: gplot, gplotd, gplotdc, distmat, ve2axy, axy2ve
%
% Author: Joseph Kirk
% Email: jdkirk630@gmail.com
% Release: 1.1
% Date: 4/29/09

% Process Inputs
error(nargchk(2,5,nargin));
all_positive = 1;
[n,nc] = size(AorV);
[m,mc] = size(xyCorE);
[E,cost] = processInputs(AorV,xyCorE);
if nargin < 5
    iswaitbar = 0;
end
if nargin < 4
    FID = (1:n);
end
if nargin < 3
    SID = (1:n);
end
if max(SID) > n || min(SID) < 1
    eval(['help ' mfilename]);
    error('Invalid [SID] input. See help notes above.');
end
if max(FID) > n || min(FID) < 1
    eval(['help ' mfilename]);
    error('Invalid [FID] input. See help notes above.');
end

isreversed = 0;
if length(FID) < length(SID)
    E = E(:,[2 1]);
    cost = cost';
    tmp = SID;
    SID = FID;
    FID = tmp;
    isreversed = 1;
end

L = length(SID);
M = length(FID);
costs = zeros(L,M);
paths = num2cell(nan(L,M));

% Find the Minimum Costs and Paths using Dijkstra's Algorithm
if iswaitbar, wbh = waitbar(0,'Please Wait ... '); end
for k = 1:L
    % Initializations
    if all_positive, TBL = sparse(1,n); else TBL = NaN(1,n); end
    min_cost = Inf(1,n);
    settled = zeros(1,n);
    path = num2cell(nan(1,n));
    I = SID(k);
    min_cost(I) = 0;
    TBL(I) = 0;
    settled(I) = 1;
    path(I) = {I};

    while any(~settled(FID))
        % Update the Table
        TAB = TBL;
        if all_positive, TBL(I) = 0; else TBL(I) = NaN; end
        nids = find(E(:,1) == I);
        % Calculate the Costs to the Neighbor Points and Record Paths
        for kk = 1:length(nids)
            J = E(nids(kk),2);
            if ~settled(J)
                c = cost(I,J);
                if all_positive, empty = ~TAB(J); else empty = isnan(TAB(J)); end
                if empty || (TAB(J) > (TAB(I) + c))
                    TBL(J) = TAB(I) + c;
                    if isreversed
                        path{J} = [J path{I}];
                    else
                        path{J} = [path{I} J];
                    end
                else
                    TBL(J) = TAB(J);
                end
            end
        end
        
        if all_positive, K = find(TBL); else K = find(~isnan(TBL)); end
        % Find the Minimum Value in the Table
        N = find(TBL(K) == min(TBL(K)));
        if isempty(N)
            break
        else
            % Settle the Minimum Value
            I = K(N(1));
            min_cost(I) = TBL(I);
            settled(I) = 1;
        end
    end
    % Store Costs and Paths
    costs(k,:) = min_cost(FID);
    paths(k,:) = path(FID);
    if iswaitbar, waitbar(k/L,wbh); end
end
if iswaitbar, close(wbh); end

if isreversed
    costs = costs';
    paths = paths';
end

if L == 1 && M == 1
    paths = paths{1};
end

% -------------------------------------------------------------------
    function [E,C] = processInputs(AorV,xyCorE)
        C = sparse(n,n);
        if n == nc
            if m == n
                if m == mc % Inputs: A,cost
                    A = AorV;
                    A = A - diag(diag(A));
                    C = xyCorE;
                    all_positive = all(C(logical(A)) > 0);
                    E = a2e(A);
                else % Inputs: A,xy
                    A = AorV;
                    A = A - diag(diag(A));
                    xy = xyCorE;
                    E = a2e(A);
                    D = ve2d(xy,E);
                    all_positive = all(D > 0);
                    for row = 1:length(D)
                        C(E(row,1),E(row,2)) = D(row);
                    end
                end
            else
                eval(['help ' mfilename]);
                error('Invalid [A,xy] or [A,cost] inputs. See help notes above.');
            end
        else
            if mc == 2 % Inputs: V,E
                V = AorV;
                E = xyCorE;
                D = ve2d(V,E);
                all_positive = all(D > 0);
                for row = 1:m
                    C(E(row,1),E(row,2)) = D(row);
                end
            elseif mc == 3 % Inputs: V,E3
                E3 = xyCorE;
                all_positive = all(E3 > 0);
                E = E3(:,1:2);
                for row = 1:m
                    C(E3(row,1),E3(row,2)) = E3(row,3);
                end
            else
                eval(['help ' mfilename]);
                error('Invalid [V,E] inputs. See help notes above.');
            end
        end
    end

    % Convert Adjacency Matrix to Edge List
    function E = a2e(A)
        [I,J] = find(A);
        E = [I J];
    end

    % Compute Euclidean Distance for Edges
    function D = ve2d(V,E)
        VI = V(E(:,1),:);
        VJ = V(E(:,2),:);
        D = sqrt(sum((VI - VJ).^2,2));
    end
end
