% Mohammad Asif Zaman
% April 26, 2019

% 2D Newton-Raphson to find the exact position of the critical points from an initial guess


function [xc,yc,err_flag] = NewtonRaphson2D(X,Y,ux,uy,J11,J12,J21,J22,x0,y0)


xc = x0;
yc = y0;
err_flag = 0;

xL = min(X(:)); xU = max(X(:)); yL = min(Y(:)); yU = max(Y(:));

dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);

max_iter = 30;      % Maximum number of iterations for the Newton-Raphson method
xy_tol = dx./100;   % Desired accuracy. If delta_x and delta_y are smaller than xy_tol, then Newton-Raphson stops.



for iter = 1:max_iter
    ux0 = interp2(X,Y,ux,x0,y0);
    uy0 = interp2(X,Y,uy,x0,y0);
    b = [ux0; uy0];
    J0 = JacobianInterp(J11,J12,J21,J22,X,Y,x0,y0);
    
    if sum(isnan(J0(:))) ~= 0
        err_flag = 2;
        fprintf('Singular matrix found while performing iterations in Newton-Raphson.\n');
        break;
    end
    delta = -J0\b;
    
    if ~isnan(delta)
        x0 = x0 + delta(1);
        y0 = y0 + delta(2);
    end
    
    % Break if the values of delta(1) and delta(2) are smaller than the
    % xy_tol.
    if abs(delta(1)) < xy_tol & abs(delta(2)) < xy_tol
        break;
    end
end


if x0 < xL & x0 > xU & y0 < yL & y0 > yU
    err_flag = 1;
elseif err_flag == 2
    err_flag = 2;
else
    err_flag = 0;
    xc = x0;
    yc = y0;
end

