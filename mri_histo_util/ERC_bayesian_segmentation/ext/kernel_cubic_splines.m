%%
% Spline-encoded field:
%     F(x,y,z) = \sum_{i,j,k} c_{i,j,k} f(x-x_i) f(y-y_j) f(z-z_j)
% where f is the cubic B-spline and c_{i,j,k} are the control coefficients 
% centered at (x_i, y_j, z_k).
%
% Membrane energy of a spline-encoded field:
%     L = Lx + Ly + Lz, where
%     Lx = 0.5 \int (\sum_{i,j,k} c_{i,j,k} g(x-x_i) f(y-y_j) f(z-z_j))^2 dx dy dz
%     Ly = 0.5 \int (\sum_{i,j,k} c_{i,j,k} f(x-x_i) g(y-y_j) f(z-z_j))^2 dx dy dz
%     Lz = 0.5 \int (\sum_{i,j,k} c_{i,j,k} f(x-x_i) f(y-y_j) g(z-z_j))^2 dx dy dz
% where g is the derivative of f.
%
% This "clearly" can be written as the quadratic c(:)' * L * c(:), where L
% is Toeplitz (i.e., a convolution kernel), assuming circulant boundaries.
% To find the kernel, all we have to do is differentiate L wrt c(:).
%
% dLx/dc_{i,j,k} 
%   = \int g(x-x_i) f(y-y_j) f(z-z_k) (\sum_{i',j',k'} c_{i',j',k'} g(x-x_i') f(y-y_j') f(z-z_k')) dx dy dz
% [get the sum out of the integral]
%   = \sum_{i',j',k'} c_{i',j',k'} \int g(x-x_i) f(y-y_j) f(z-z_k) g(x-x_i') f(y-y_j') f(z-z_k') dx dy dz
% [separable in x/y/z]
%   = \sum_{i',j',k'} c_{i',j',k'} (\int g(x-x_i) g(x-x_i') dx) (\int f(y-y_j) f(y-y_j') dy) (\int f(z-z_k) f(z-z_k') dz)
% [clearly a convolution]
% 
% We just need to compute the terms:
%   * \int g(x-x_i) g(x-x_i') dx
%   * \int f(y-y_j) f(y-y_j') dy
%   * \int f(z-z_k) f(z-z_k') dz
% and put the kernel together. Since the B-splines have a small support,
% the kernel is small (7x7x7).
% Note that the convolution kernels corresponding to Ly and Lz can be 
% computed in a similar way.
%
% If we want the bending energy (second derivatives) instead of the 
% membrane energy (first derivatives), we can put second derivatives in the
% initial integral and process similarly.


%%
syms x 'real'

%% functions
f3_lo(x)    = (x.*x.*(x-2.).*3 + 4.) ./ 6.;
f3_up(x)    = (2. - x).^3 ./ 6.;
g3_lo(x)    = x.*(3.*x - 4.) ./ 2.;
g3_up(x)    = -(2. - x).^2 ./ 2.;
h3_lo(x)    = (3.*x - 2.);
h3_up(x)    = (2. - x);
f3_pos(x)   = piecewise(x < 2, piecewise(x < 1, f3_lo(x), f3_up(x)), 0);
g3_pos(x)   = piecewise(x < 2, piecewise(x < 1, g3_lo(x), g3_up(x)), 0);
h3_pos(x)   = piecewise(x < 2, piecewise(x < 1, h3_lo(x), h3_up(x)), 0);
f3(x)       = f3_pos(abs(x));
g3(x)       = g3_pos(abs(x)) * sign(x);
h3(x)       = h3_pos(abs(x));

%% check derivatives

assert(simplify(g3_lo(x) - diff(f3_lo(x), x), 1000) == 0);
assert(simplify(g3_up(x) - diff(f3_up(x), x), 1000) == 0);

assert(simplify(h3_lo(x) - diff(g3_lo(x), x), 1000) == 0);
assert(simplify(h3_up(x) - diff(g3_up(x), x), 1000) == 0);

%% integrate
lim = 2;

F0 = int(f3(x) * f3(x+0), x, -lim, lim);
F1 = int(f3(x) * f3(x+1), x, -lim, lim);
F2 = int(f3(x) * f3(x+2), x, -lim, lim);
F3 = int(f3(x) * f3(x+3), x, -lim, lim);
F4 = int(f3(x) * f3(x+4), x, -lim, lim);
assert(F4 == 0)

G0 = int(g3(x) * g3(x+0), x, -lim, lim);
G1 = int(g3(x) * g3(x+1), x, -lim, lim);
G2 = int(g3(x) * g3(x+2), x, -lim, lim);
G3 = int(g3(x) * g3(x+3), x, -lim, lim);
G4 = int(g3(x) * g3(x+4), x, -lim, lim);
assert(G4 == 0)

H0 = int(h3(x) * h3(x+0), x, -lim, lim);
H1 = int(h3(x) * h3(x+1), x, -lim, lim);
H2 = int(h3(x) * h3(x+2), x, -lim, lim);
H3 = int(h3(x) * h3(x+3), x, -lim, lim);
H4 = int(h3(x) * h3(x+4), x, -lim, lim);
assert(H4 == 0)

fprintf("F0 = %s\n", rats(F0));
fprintf("F1 = %s\n", rats(F1));
fprintf("F2 = %s\n", rats(F2));
fprintf("F3 = %s\n", rats(F3));
fprintf("G0 = %s\n", rats(G0));
fprintf("G1 = %s\n", rats(G1));
fprintf("G2 = %s\n", rats(G2));
fprintf("G3 = %s\n", rats(G3));
fprintf("H0 = %s\n", rats(H0));
fprintf("H1 = %s\n", rats(H1));
fprintf("H2 = %s\n", rats(H2));
fprintf("H3 = %s\n", rats(H3));

F = double([F3 F2 F1 F0 F1 F2 F3]);
G = double([G3 G2 G1 G0 G1 G2 G3]);
H = double([H3 H2 H1 H0 H1 H2 H3]);


% putting it together (Lx + Ly + Lz)
kernel_membrane_2d = G' .* F + F' .* G;
kernel_membrane_3d = G' .* F .* reshape(F, 1, 1, []) ...
                   + F' .* G .* reshape(F, 1, 1, []) ...
                   + F' .* F .* reshape(G, 1, 1, []);

kernel_bending_2d = H' .* F + F' .* H + 2 * G' .* G;
kernel_bending_3d = H' .* F .* reshape(F, 1, 1, []) ...
                  + F' .* H .* reshape(F, 1, 1, []) ...
                  + F' .* F .* reshape(H, 1, 1, []) ...
                  + F' .* G .* reshape(G, 1, 1, []) * 2 ...
                  + G' .* F .* reshape(G, 1, 1, []) * 2 ...
                  + G' .* G .* reshape(F, 1, 1, []) * 2;

% plot
redblue = [linspace(0, 1, 128)' linspace(0, 1, 128)' linspace(1, 1, 128)'
           linspace(1, 1, 128)' linspace(1, 0, 128)' linspace(1, 0, 128)'];

figure
imagesc(kernel_membrane_2d)
clim([-max(abs(kernel_membrane_2d(:))), max(abs(kernel_membrane_2d(:)))])
colormap(redblue)
colorbar
title('membrane')

figure
imagesc(kernel_bending_2d)
clim([-max(abs(kernel_bending_2d(:))), max(abs(kernel_bending_2d(:)))])
colormap(redblue)
title('bending')
colorbar