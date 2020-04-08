% Mohammad Asif Zaman
% April 7, 2020
% Testing the Newton-Raphson root finding funciton

clear all;
close all;
clc; clf;


x = linspace(-10,10,31);
y = linspace(-5,5,21);

[X,Y] = meshgrid(x,y);

dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);




ux = X.^2 - (Y-2).^2;
uy = sin(.1*X) + cos(.5*(Y-2)+pi/2);

contourf(x,y,ux,40,'linestyle','none');
figure,
contourf(x,y,uy,40,'linestyle','none');


[Ny,Nx] = size(ux);     % Size of the vector field

% D2x = d/dx, DD2x = d^2/dx^2
[D2x,D2y,DD2x,DD2y] = F_diff_mat_2D_v2(Ny,Nx);  % Loading the finite difference matrices

J11 = reshape(D2x*ux(:)./(2*dx), size(X));  % dux/dx
J12 = reshape(D2y*ux(:)./(2*dy), size(X));  % dux/dy
J21 = reshape(D2x*uy(:)./(2*dx), size(X));  % duy/dx
J22 = reshape(D2y*uy(:)./(2*dy), size(X));  % dut/dy

x0 = -1;
y0 = 1;

[xc,yc,err_flag] = NewtonRaphson2D(X,Y,ux,uy,J11,J12,J21,J22,x0,y0);
fprintf('Obtained root, xc = %1.2f, yc = %1.2f\n',xc,yc);
