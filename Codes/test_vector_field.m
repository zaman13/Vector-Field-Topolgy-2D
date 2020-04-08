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
% 
% ux = exp(-X2.^2 - Y2.^2).*(-2*X2.^2 + 1);
% uy = exp(-X2.^2 - Y2.^2).*(-2*Y2.*X2);
% 

% Test field 1
ux = Y.^2 - X.^2;
uy = X.^2 + Y.^2 - 2;

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




