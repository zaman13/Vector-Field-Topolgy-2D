% Mohammad Asif Zaman
% April 24, 2019

%  Functio to generate a sample 2D vector field


function [ux, uy] = test_vector_field();

Nt = 40;
xL = -3;
xU = 3;
yL = -2;
yU = 2;
x = linspace(xL,xU,Nt);
y = linspace(yL,yU,Nt);

[X2,Y2] = meshgrid(x,y);

% ux = exp(-X.^2 - Y.^2).*(-2*X.^2 + 1);
% uy = exp(-X.^2 - Y.^2).*(-2*Y.*X);
% 
ux = Y2.^2 - X2.^2;
uy = X2.^2 + Y2.^2 - 2;

% size(ux)
% size(X2)
% ux = Y.^2 - sin(0.5*pi*X.^2);
% uy = X.^2 + Y.^2 - 2;

% ux = -Y;
% uy = 3*X;

