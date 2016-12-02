function J = computeRegularizers(smoothingType, numBFuncs, domain, voxelSize)

smoothingType = lower(smoothingType);

switch smoothingType
    case 'bspline'
        J = regularizerBSpline(numBFuncs);
    case 'spm'             
        J = regularizerSPM(numBFuncs, domain, voxelSize);
    otherwise
        error('Unsupported basis function smoothing type');
end

end

% Penalizes curvature by implicit regularization of the second order derivatives only in the diagonal
% The second order derivatives are never actually used,
% But exp is used instead which elegantly ensures that we penalize more on higher order cosines.
% Then, taking the kronecker product of that ensures that the combined basis functions are further penalized, given
% the product of basis functions
function J = regularizerSPM(numBFuncs, domain, voxelSize)

sd = 2 * (domain(:, 2) - domain(:, 1));

% sd is the standard deviation of the cosines, d3 contains the number of basis functions in each dimension and krn_?
% represents the second-order derivatives we want to penalize
krn_x = exp(-(0:(numBFuncs(1)-1)).^2/sd(1).^2) / sqrt(voxelSize(1));
krn_y = exp(-(0:(numBFuncs(2)-1)).^2/sd(2).^2) / sqrt(voxelSize(2));
krn_z = exp(-(0:(numBFuncs(3)-1)).^2/sd(3).^2) / sqrt(voxelSize(3));

% This is how we form the regularizer with the kronecker product
J = diag(kron(krn_z, kron(krn_y, krn_x)).^(-2));
%J = sparse(1:length(J), 1:length(J), J, length(J), length(J));
end

% thin-plate spline regularization (N3 way), which penalizes curvature given the second order derivatives, + adds boundary
% conditions for the remaining partial derivatives
% TODO: Verify that this implementation is properly generalized
function J = regularizerBSpline(numBFuncs)
nDimensions = length(numBFuncs);

%  deal with special case first
if(nDimensions == 1)
    J = bendingEnergyCubic(numBFuncs, 2);
    return;
end

% compute 1-D bending energy matrices
bm = cell(nDimensions,3);
order = 2;

for i = 1:nDimensions
    for j = 0:order
        bm{i,j+1} = bendingEnergyCubic(numBFuncs(i), j);
    end
end

% form kronecker products
% eg. in 3D   x"yz + xy"z + xyz" + 2x'y'z + 2xy'z' + 2x'yz'
J = kron(bm{3,3}, kron(bm{2,1},bm{1,1})) + kron(bm{3,1}, kron(bm{2,3},bm{1,1})) + kron(bm{3,1}, kron(bm{2,1},bm{1,3})) ...
    + 2*kron(bm{3,2}, kron(bm{2,2},bm{1,1})) + 2*kron(bm{3,2}, kron(bm{2,1},bm{1,2})) + 2*kron(bm{3,1}, kron(bm{2,2},bm{1,2}));
J = J';
end

function energy = bendingEnergyCubic(s, order)

if(s < 4 || order < 0 || order > 2)
    error('bending energy not defined for s < 4');
end

% standardized cubic B-spline defined on [-4 0] with knots at -4 -3 -2 -1 0
% each unit segment is written as a cubic defined on [0 1]  eg. -x^3 + 3x^2 -3x + 1
B = [-1,  3, -3, 1; ...
    3, -6,  0, 4; ...
    -3,  3,  3, 1; ...
    1,  0,  0, 0];

% take order'th derivative of B
for i = 1:order
    for j = 0:3
        for k = 1:3
            B(j+1, 4-k+1) = k*B(j+1,4-k);
        end
        B(j+1,1) = 0;
    end
end

% compute product integral for each pair of segments actually we only compute upper triangle since D will be symmetric
C = zeros(7,1);
sizeC = length(C);
D = zeros(4);

for i = 0:3
    for j = i:3
        % convolve each pair of polynomials
        for k = 0:(sizeC-1)
            C(k+1) = 0;
            for l = 0:(min(k+1,sizeC-k)-1)
                C(k+1) = C(k+1) + B(i+1, max(0,k-4+1)+l+1) * B(j+1, 4-1-l-max(0,4-1-k)+1);
            end
        end
        
        for k = 0:(sizeC-1)
            D(i+1,j+1) = D(i+1,j+1) + (C(k+1) / (sizeC-k));
        end
    end
end


% define 6 regions:  [-1,0] [-2,0] [-3,0] [-4,0] [-2,-1] [-3,-1] splines can be shifted with respect to each by 0 1 2 or 3
% Compute product interval on each interval for each offset
integ = zeros(6,4);

% integrals for region one to four
for offset = 0:3
    for region = 0:(4-offset-1)
        for interval = 0:region
            integ(region+1,offset+1) = integ(region+1,offset+1) + D(interval+1,interval+offset+1);
        end
    end
    for region = (4-offset):3
        integ(region+1,offset+1) = integ(region, offset+1);
    end
end

% compute integrals for regions five and six separately
integ(5,1) = D(2,2);
integ(5,2) = D(2,3);
integ(5,3) = D(2,4);
integ(5,4) = 0;

integ(6,1) = D(2,2) + D(3,3);
integ(6,2) = D(2,3) + D(3,4);
integ(6,3) = D(2,4);
integ(6,4) = 0;

% form integrals into bending energy matrix
energy = zeros(s);

% set elements common to matrices of all sizes first
for i = 0:2
    
    energy(1,i+1) = integ(1,i+1);
    energy(i+1,1) = integ(1,i+1);
    energy(s,s-i) = integ(1,i+1);
    energy(s-i,s) = integ(1,i+1);
    
end

if(s == 4)
    % deal with special cases
    energy(2,2) = integ(5,1);
    energy(3,3) = integ(5,1);
    energy(2,3) = integ(5,2);
    energy(3,2) = integ(5,2);
    energy(4,1) = integ(4,4);
    energy(1,4) = integ(4,4);
    
elseif(s == 5)
    % deal with special cases
    energy(2,2) = integ(2,1);
    energy(4,4) = integ(2,1);
    energy(3,3) = integ(6,1);
    energy(3,2) = integ(2,2);
    energy(2,3) = integ(2,2);
    energy(3,4) = integ(2,2);
    energy(4,3) = integ(2,2);
    energy(2,4) = integ(4,3);
    energy(4,2) = integ(4,3);
    
    for i = 0:1
        energy(i+1,i+4) = integ(4,4);
        energy(i+4,i+1) = integ(4,4);
    end
else
    % n > 5  (general case)
    % fill in bulk region
    for j = 0:3
        for i = (3-j):(s-4)
            
            energy(i+1,i+j+1) = integ(4,j+1);
            energy(i+j+1,i+1) = integ(4,j+1);
        end
    end
    
    % fill in boundary region
    energy(3,3) = integ(3,1);
    energy(s-2,s-2) = integ(3,1);
    energy(2,3) = integ(2,2);
    energy(3,2) = integ(2,2);
    energy(s-1,s-2) = integ(2,2);
    energy(s-2,s-1) = integ(2,2);
    energy(2,2) = integ(2,1);
    energy(s-1,s-1) = integ(2,1);
end

end
