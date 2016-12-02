function smoothField = expandSeparableBasisFunctions(C, B)
%
%
%

[Nx, Nbx] = size(B{1});
[Ny, Nby] = size(B{2});
[Nz, Nbz] = size(B{3});

% Expand the basis into the image domain
% x direction
tmp = permute(reshape(B{1} * reshape(C, [Nbx, Nby * Nbz]), [Nx, Nby, Nbz]), [2, 1, 3]);

% y direction
tmp = permute(reshape(B{2} * reshape(tmp, [Nby, Nx * Nbz]), [Ny, Nx, Nbz]), [2, 1, 3]);

% z direction
smoothField = reshape(reshape(tmp, [Nx * Ny, Nbz]) * B{3}', [Nx, Ny, Nz]);