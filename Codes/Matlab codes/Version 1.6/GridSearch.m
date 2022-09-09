% Mohammad Asif Zaman
% April 26, 2019

% Grid search to find the exact position of the critical points from an initial guess


function [xc,yc] = GridSearch(X,Y,ux,uy,x0,y0)

dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);


xc = x0;
yc = y0;

max_iter = 10;



% for m = 1:max_iter
%     
%     
%     
%     
%     xL = xc - 2*dx;
%     xU = xc + 2*dx;
%     yL = yc - 2*dy;
%     yU = yc + 2*dy;
%     
%     
%     
%     
%     xfine = linspace(xL,xU,30);
%     yfine = linspace(yL,yU,30);
%     
%     
%     [Xfine,Yfine] = meshgrid(xfine,yfine);
%     
%     Ufine = interp2(X,Y,u2,Xfine,Yfine);
%     ind = find(Ufine == min(Ufine(:)));
%     
%     xc = Xfine(ind);
%     yc = Yfine(ind);
%     
%     
%     dx = xfine(2) - xfine(1);
%     dy = yfine(2) - yfine(1);
%     
% end



for m = 1:max_iter
    
    
    
    
    xL = xc - 2*dx;
    xU = xc + 2*dx;
    yL = yc - 2*dy;
    yU = yc + 2*dy;
    
    
    
    
    xfine = linspace(xL,xU,30);
    yfine = linspace(yL,yU,30);
    
    
    [Xfine,Yfine] = meshgrid(xfine,yfine);
    
    Uxfine = interp2(X,Y,ux,Xfine,Yfine);
    Uyfine = interp2(X,Y,uy,Xfine,Yfine);
    
    [ind_ux_r,ind_ux_c] = find(Uxfine == min(Uxfine(:)));
    [ind_uy_r,ind_uy_c] = find(Uyfine == min(Uyfine(:)));
     
    
    xc = 0.5*(Xfine(ind_ux_r,ind_ux_c) + Xfine(ind_uy_r,ind_uy_c));
    yc = 0.5*(Yfine(ind_ux_r,ind_ux_c) + Yfine(ind_uy_r,ind_uy_c));
    
    
    dx = xfine(2) - xfine(1);
    dy = yfine(2) - yfine(1);
    
end