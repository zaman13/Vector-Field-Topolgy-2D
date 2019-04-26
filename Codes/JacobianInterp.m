% Mohammad Asif Zaman
% April 23, 2019


% Calculating the Jacobi matrix at an arbitrary location from interpolation

% Input: X,Y =             2D meshgrid space variables
%
%        J11,J12,J21,J22 = Elements of the Jacobian matrix at X,Y positions.
%                          Each of them are matrices of the same size as X.
%
%        xc,yc =           coordinates of the point where the Jacobian matrix must be
%                          calculated. Note that xc and yc aren't arrays.


function Jout = JacobianInterp(J11,J12,J21,J22, X,Y, xc,yc)



temp_J11 = interp2(X,Y,J11,xc,yc);
temp_J12 = interp2(X,Y,J12,xc,yc);
temp_J21 = interp2(X,Y,J21,xc,yc);
temp_J22 = interp2(X,Y,J22,xc,yc);

Jout = [temp_J11 temp_J12; temp_J21 temp_J22];





