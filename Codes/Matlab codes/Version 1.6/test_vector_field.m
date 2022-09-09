% Mohammad Asif Zaman
% April 24, 2019

%  Function to generate a sample 2D vector field


function [X,Y, ux, uy] = test_vector_field();

Nt = 40;
xL = -3;
xU = 3;
yL = -2;
yU = 2;
x = linspace(xL,xU,Nt);
y = linspace(yL,yU,Nt);

[X,Y] = meshgrid(x,y);

dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);
[Ny,Nx] = size(X); 
[D2x,D2y,DD2x,DD2y] = F_diff_mat_2D_v2(Ny,Nx);  % Loading the finite difference matrices

% 
% ux = exp(-X2.^2 - Y2.^2).*(-2*X2.^2 + 1);
% uy = exp(-X2.^2 - Y2.^2).*(-2*Y2.*X2);
% 

% Test field 1
ux = Y.^2 - X.^2;
uy = X.^2 + Y.^2 - 2;

% Test field 2: Gradient of peaks function
% V = peaks(X,Y);
% ux = reshape(D2x*V(:)./(2*dx), size(X)) + 4;
% uy = reshape(D2y*V(:)./(2*dy), size(X)) - 4;




% Test field 3: Gradient field
% V1 = exp(-3*(.3*Y-.1).^2 - 3*(.3*X-.3).^2);
% V2 = exp(-4*(.3*Y-.2).^2 - 4*(.3*X+1.2).^2);
% V3 = exp(-2*(.5*Y+.2).^2 - 2*(.5*X+1.2).^2);
% V4 = -exp(-2*(.5*X-.4).^2-2*(.5*Y-.4).^2);
% V = V1 + V2 + V3 + V4;
% 
% % contourf(X,Y,V,40,'linestyle','none'); colorbar; title('V'); xlabel('x'); ylabel('y');
% 
% ux = reshape(D2x*V(:)./(2*dx), size(X));
% uy = reshape(D2y*V(:)./(2*dy), size(X));

% Test field 2
% ux = -5*X.^2 - 6*(Y.^2-2).^2;
% uy = exp(-2*X.^2 - 3*(Y.^2-2).^2) -1;



% Test field 2
% ux = X.^2 - (Y-2).^2;
% uy = sin(.1*X) + cos(.5*(Y-2)+pi/2);


% size(ux)
% size(X2)

% Test field 3
% ux = Y.^2 - sin(0.5*pi*X.^2);
% uy = X.^2 + Y.^2 - 2;




