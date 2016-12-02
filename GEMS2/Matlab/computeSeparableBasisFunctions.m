function [B, B2, indexMapB2] = computeSeparableBasisFunctions(numBFuncs, domain, distance, dataSize, smoothingType,  step, start)

if nargin < 5
    smoothingType = 'spm';
end
if nargin < 6
    step = [1, 1, 1];
end
if nargin < 7
    start = [0, 0, 0];
end

numDims = length(dataSize);
B = cell(numDims, 1);
B2 = cell(numDims, 1);

smoothingType = lower(smoothingType);

switch smoothingType
    case 'bspline'
        % now we prepare the knot vectors
        knots = zeros(max(numBFuncs) + 4, numDims);
        
        % knot positions along each dimension in mm. Note how the domain is padded on both sides with this approach
        for i = 1:numDims
            knots(1:(numBFuncs(i)+4),i) = (0.5 * (domain(i, 1) + domain(i, 2) - distance * (numBFuncs(i) + 3))) + distance * (0:(numBFuncs(i)+3));
        end
        
        % now we evaluate all 1D b-splines along the voxel centers for each dimension
        for i = 1:numDims
            % one cell covers all voxels evaluated for each of the B-Splines
            B{i} = zeros(dataSize(i), numBFuncs(i));
            
            x = start(i) + step(i)*(0:(dataSize(i)-1));
            
            % doing it the N3 way
            normalizer =  1 / (distance^3);
            for k = 1:length(x)
                nearest = ceil((x(k) - knots(4,i)) ./ distance);
                temp = normalizer * (knots(nearest+4,i) - x(k))^3;
                B{i}(k,nearest) = temp;
                B{i}(k,nearest+1) = normalizer * (knots(nearest+5,i) - x(k))^3 -4*temp;
                temp = normalizer * (x(k) - knots(nearest+3,i))^3;
                B{i}(k,nearest+2) = normalizer * (x(k) - knots(nearest+2,i))^3 - 4*temp;
                B{i}(k,nearest+3) = temp;
            end
            %             % precompute factorial values for the B splines
            %             fac = factorial(4) ./ (factorial(0:4).*factorial(4:-1:0));
            %             % irving-hall pdf (N3 doesn't do this)
            %             normalizer =  1 / (2 * fac(3) * (dist^3));
            %
            %             % the 2 stems from doing mu = -1, if val < 0, instead of mu = 0, if val < 0
            %             normalizer =  1 / (2 * (dist^3));
            %             for j = 1:numBFuncs(i)
            %                 for s = 0:4
            %                     val = x - knots(j+4-s, i);
            %
            %                     % mus are not defined properly in article - but seems to be properly implemented in code
            %                     % fixed here
            %                     mu = zeros(size(x));
            %                     mu(val > 0) = 1;
            %                     mu(val < 0) = -1;
            %
            %                     B{i}(indices, j) = B{i}(indices, j) + ((-1)^s) * normalizer * fac(s+1) * (val.^3) .* mu;
            %                 end
            %             end
            
            % doing stuff in closed form + normalized - for verification
            %             for j = 1:numBFuncs(i)
            %                 val = (x - knots(j,i)) / dist;
            %                 nonzero = ceil(val) == 1;
            %                 B{i}(nonzero, j) = val(nonzero).^3;
            %
            %                 val = (x - knots(j+1,i)) / dist;
            %                 nonzero = ceil(val) == 1;
            %                 B{i}(nonzero, j) = -3*val(nonzero).^3 + 3*val(nonzero).^2 + 3*val(nonzero) + 1;
            %
            %                 val = (x - knots(j+2,i)) / dist;
            %                 nonzero = ceil(val) == 1;
            %                 B{i}(nonzero, j) = 3*val(nonzero).^3 -6*val(nonzero).^2 + 4;
            %
            %                 val = (x - knots(j+3,i)) / dist;
            %                 nonzero = ceil(val) == 1;
            %                 B{i}(nonzero, j) = -val(nonzero).^3 + 3*val(nonzero).^2 - 3*val(nonzero) + 1;
            %             end
        end
    case 'spm'
        % Construct DCT basis functions Let's use the definition of Matlab:
        % http://www.mathworks.se/help/signal/ref/dct.html
        
        % according to my current implementation, based on the way N3 does things, I always assume that the basis
        % functions are evaluated at the center of the voxel. This means that the cosines actually extend a bit to both
        % sides of the voxels at the beginning and at the end - which makes sense intuitively
        % This is different from the old cosine code - but more elegant and more flexible, since we depend on the step
        % size
        
        % Note: First column weight 1 / sqrt(Nx) == sqrt(2 / Nx) * oneOver2sqrt
        oneOverSqrt2 = 1 / sqrt(2);
        for i = 1:numDims
            
            % compute the start + end of span to evaluate basis functions in in
            indices = 1:dataSize(i);
            x = start(i) + indices*step(i);
            realStart = domain(i,1);
            realEnd = (domain(i, 2) - domain(i, 1));
            
            % using step size implicitly, but same as from webpage
            tmp = (2 * x' - realStart) * pi * 0.5 / realEnd;
            
            %             % straight from the webpage
            %             tmp2 = (2 * (1:N(i))' - 1) * pi * 0.5 / N(i);
            %
            %             % old code for expanding back - seems like a bug...
            %             tmp3 = (2 * (((1:(N(i)*4))' - 1) * 0.25 + 1) - 1) * 0.5 * pi / N(i);
            
            B{i} = cos(tmp * (0:numBFuncs(i)-1)) * sqrt(2 / realEnd);
            B{i}(:, 1) = B{i}(:, 1) * oneOverSqrt2;
        end
    otherwise
        error('Unsupported basis function smoothing type');
end



% Exploit the separability of the basis functions along each dimension for building AWA
for d = 1:numDims
    B2{d} = zeros(dataSize(d), numBFuncs(d));
    for i = 1:numBFuncs(d)
        shapeI = B{d}(:, i);
        for j = 1:numBFuncs(d)
            shapeJ = B{d}(:, j);
            B2{d}(:, i + (j-1) * numBFuncs(d)) = shapeI .* shapeJ;
        end
    end
end

Nb = numBFuncs;

% get the index map between the 3D volume and the 2D matrix
indexMapB2 = zeros(prod(Nb)^2,1);

count = 1;
% map to columns Z
for i = 0:Nb(3)-1
    % map to rows Z
    for j = 0:Nb(3)-1
        % map to columns Y
        for k = 0:Nb(2)-1
            % map to rows Y
            for l = 0:Nb(2)-1
                % map to columns X
                for m = 0:Nb(1)-1
                    % map to rows X
                    for n = 0:Nb(1)-1
                        % we're essentially mapping a 3D volume, where each dimension contains all combinations of
                        % basis functions for that particular dimension, to the corresponding 2D matrix where each
                        % dimension contains the possible combinations between the basis functions in three
                        % dimensions.
                        
                        % now, lets assume that both matrices are mapped to their corresponding 1D array. This means
                        % that the 3D volume will have the combinations follow the structure of this 6D loop, where
                        % each pair of loops correspond the the combined basis funcs along one dimension.
                        % the 2D volume then follows indices that combine given every second loop, corresponding to
                        % rows/columns respectively.
                        
                        % indexMap(count) = i*X*Y*Z*X*Y + k*X*Y*Z*X + m*X*Y*Z + j*X*Y + l*X + n;
                        indexMapB2(count) = Nb(1)*(Nb(2)*(Nb(3)*(Nb(1)*(Nb(2)*i + k) + m) + j) + l) + n;
                        count = count + 1;
                    end
                end
            end
        end
    end
end
% conform to Matlab indexing
indexMapB2 = indexMapB2 + 1;

% end of computeSeparableBasisFunctions
end

