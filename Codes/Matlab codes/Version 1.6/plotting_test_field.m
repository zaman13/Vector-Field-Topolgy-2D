clear all;
close all;
clc; clf;

%  Loading the vector field and defining the spatial grid
% ======================================================================>>>

[X,Y,ux,uy] = test_vector_field();     % Load a test vector field


% Lower and upper bound of x and y
xL = min(X(:)); xU = max(X(:)); 
yL = min(Y(:)); yU = max(Y(:)); 

[Ny,Nx] = size(ux);     % Size of the vector field


u2 = ux.^2 + uy.^2;     % Squared norm function

% figure,quiver(X,Y,ux,uy);
% hold on;
% startx = -3:.5:3;
% starty = -3*ones(size(startx));
% streamline(X,Y,ux,uy,startx,starty);
% startx = -3:.5:3;
% starty = 1*ones(size(startx));
% streamline(X,Y,ux,uy,startx,starty);

figure,
contourf(X,Y,ux,40,'linestyle','none'); colorbar; title('ux'); xlabel('x'); ylabel('y');
caxis([-max(abs(ux(:))) max(abs(ux(:)))]);
colormap(jet);
figure,
contourf(X,Y,uy,40,'linestyle','none'); colorbar; title('uy'); xlabel('x'); ylabel('y');
caxis([-max(abs(uy(:))) max(abs(uy(:)))]);
colormap(jet);
figure,
contourf(X,Y,ux.^2 + uy.^2,40,'linestyle','none'); colorbar; title('|u|^2'); xlabel('x'); ylabel('y');
colormap(jet);
% <<<======================================================================

