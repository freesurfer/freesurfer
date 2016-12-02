function [numBFuncs, domain] = computeNumberOfSeparableBasisFunctions(dataSize, step, smoothingType, distance)

numDims = length(dataSize);
smoothingType = lower(smoothingType);

% Compute the actual domain (in mm) that the data spans. The beginning and ending of the domain needs to be defined for the spline knots
% we're essentially shifting the domain, since we assume that a voxel has its origin at its center
% so assuming  a 'normal' 1 mm^3 volume, we're starting in -0.5 and running until 199.5 for each dimension
domain = zeros(numDims, 2);
for i = 1:numDims
    if(step(i) > 0)
        domain(i,1) = -0.5 * step(i);
        domain(i,2) = (dataSize(i) - 0.5) * step(i);
    else
        domain(i,2) = -0.5 * step(i);
        domain(i,1) = (dataSize(i) - 0.5) * step(i);
    end
end

switch smoothingType
    case 'bspline'
        % this is actually cheating a bit by the original authors, since they add 3 additional bsplines (more flexibility)
        % for this example: 200 / 50 + 3 = 7
        % well, they do it to make sure that there's full flexibility over the domain, i.e. the splines extend to both
        % sides of the domain
        % I have no clue why they do the distance * (1 + eps)... shouldn't make any practical difference, since they ceil anyway
        numBFuncs = ceil((domain(:, 2) - domain(:, 1)) ./ distance) + 3;
    case 'spm'
        %numBFuncs = round((2 * (dataSize .* stepSize) / dist + 1) * 0.5);
        numBFuncs = ceil( 2 * (domain(:, 2) - domain(:, 1)) ./ distance);
    otherwise
        error('Unsupported basis function smoothing type');
end

end