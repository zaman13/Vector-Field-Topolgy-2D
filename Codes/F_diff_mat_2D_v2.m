function [D2x,D2y,DD2x,DD2y] = F_diff_mat_2D_v2(Nx,Ny)

ex = ones(Nx,1);
ey = ones(Ny,1);

Ix = spdiags(ex,0,Nx,Nx);
Iy = spdiags(ey,0,Ny,Ny);



% Constructing the 1D diffentiation matrices
%==============================================

% It is necessary to develop different matrixed for D1dx and D1dy as Nx and
% Ny may have different values. So, the dimensions for the two cases
% might be different. 

% The main diagonals are from 2nd order center difference formula
D1x = spdiags([-ex ex], [-1,1], Nx, Nx); % A division by (2*dx) is required later.
D1y = spdiags([-ey ey], [-1,1], Ny, Ny); % A division by (2*dy) is required later.
    
D1x(1,1:3) = [-3 4 -1];   % this is 2nd order forward difference
D1x(end,end-2:end) = [1 -4 3];  % this is 2nd order backward difference

D1y(1,1:3) = [-3 4 -1];   % this is 2nd order forward difference
D1y(end,end-2:end) = [1 -4 3];  % this is 2nd order backward difference


% Second derivatives
DDx =  spdiags([ex -2*ex ex], [-1,0,1], Nx, Nx);
DDx(1,1:4) = [2 -5 4 -1];   % this is 2nd order forward difference
DDx(end,end-3:end) = [-1 4 -5 2];  % this is 2nd order backward difference

DDy =  spdiags([ey -2*ey ey], [-1,0,1], Ny, Ny);
DDy(1,1:4) = [2 -5 4 -1];   % this is 2nd order forward difference
DDy(end,end-3:end) = [-1 4 -5 2];  % this is 2nd order backward difference




% Visualization of the sparse matrices
% Should see the same thing (but with different dimensions) 
% subplot(121),spy(D1x)
% subplot(122),spy(D1y)
% A = full(D1d)

% Constructing the 2D differentiation matrices
% =============================================
% D1x = ones(Nx,Nx);
% D1y = ones(Ny,Ny);
D2x = kron(D1x,Iy);
D2y = kron(Ix,D1y);

DD2x = kron(DDx,Iy);
DD2y = kron(Ix,DDy);

