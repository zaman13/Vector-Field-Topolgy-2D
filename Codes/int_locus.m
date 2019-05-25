% Mohammad Asif Zaman
% April 23, 2019


% Calculating the the locus of the integral lines

% Input: X,Y   =           2D meshgrid space variables
%
%        ux,uy =           x and y components of the field
%        x0,y0 =           Starging point (critical point)
%        xslope,yslope =   Initial slope of the integral lines (parallel to
%                          the eigen vectors at critical point)
%        xsink,ysink =     Coordinates of the sink points where the lines
%                          terminate
%        drct =            Direction pointer. +1 means the integral lines
%                          follow the field lines. -1 means that the
%                          integral lines follow the opposite direction of
%                          the field lines.


function xylocus = int_locus(X,Y,ux,uy,x0,y0,xslope,yslope,xsink,ysink,drct)

xU = max(X(:));
xL = min(X(:));
yU = max(Y(:));
yL = min(Y(:));


dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);

red_fct = 0.2*drct;   % red_fct * dx = length of locus line-segment
vic_tol = 0.01;       % vic_tol * dx = maximum radial distance from a sink point where the locus terminates
max_length = 220;    % Maximum no. of locus points

xlocus(1) = x0;
ylocus(1) = y0;


n = 0;
flag = 1;

% Calculating locus points. The loops ends when the locus hits the domain
% boundary or a sink point.
while flag & n < max_length
    n = n + 1;
    xlocus(n+1) = xlocus(n) + xslope*dx*red_fct;
    ylocus(n+1) = ylocus(n) + yslope*dy*red_fct;
    
    
    
    xslope = interp2(X,Y,ux,xlocus(n+1),ylocus(n+1));
    yslope = interp2(X,Y,uy,xlocus(n+1),ylocus(n+1));
    
    
    if xlocus(n+1) >= xU | xlocus(n+1) <= xL | ylocus(n+1) >= yU | ylocus(n+1) <= yL
        flag = 0;
    end
    
    for snk = 1:length(xsink)
        dist = sqrt((xlocus(n+1) - xsink(snk))^2 + (ylocus(n+1) - ysink(snk))^2);
        if dist < dx*vic_tol
            flag = 0;
        end
    end
    
end

xylocus(:,1) = xlocus';
xylocus(:,2) = ylocus';

