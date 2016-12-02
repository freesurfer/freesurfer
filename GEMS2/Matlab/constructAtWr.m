function AtWr = constructAtWr(data, weights, B)

% Modified by: Christian Thode Larsen (cthla), christian.thode.larsen@gmail.com

[Nx, NumBFuncsX] = size(B{1});
[Ny, NumBFuncsY] = size(B{2});
[Nz, NumBFuncsZ] = size(B{3});

% Calculate coefficients using separability of basis functions to speed things up
% x direction
tmp = permute(reshape(B{1}' * reshape(sum(data .* weights,4), [Nx, Ny * Nz]), [NumBFuncsX, Ny, Nz]), [2, 1, 3]);

% y direction
tmp = permute(reshape(B{2}' * reshape(tmp, [Ny, NumBFuncsX*Nz]), [NumBFuncsY, NumBFuncsX, Nz]), [2, 1, 3]);

% z direction
AtWr = reshape(reshape(tmp, [NumBFuncsX * NumBFuncsY, Nz]) * B{3}, [NumBFuncsX, NumBFuncsY, NumBFuncsZ]);

AtWr = AtWr(:);
 % cthla: I'm abusing the fact that the smoothing operation is really the same, but for the weights, the hessian needs to be recomputed.
 % No need to do that here though, so hessian is provided as R   
% if(isempty(Rtranspose))
%     C(:) = R \ C(:);
% else
%     C(:) = R \ (Rtranspose \ C(:));
% end



