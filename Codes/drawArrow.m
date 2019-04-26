% Mohammad Asif Zaman
% April 25, 2019

% Function to plot arrows along a curved line.
% Input: x,y = arrays that describe the curve.

function yout = drawArrow(x,y,X,Y,ux,uy,clr)

arrow_pos_fct = 0.3;
arrow_len_fct = 0.04;

% Calculate the length of the curve
% arc length = sum(dx^2 + dy^2);

s = 0;   % s = arc length

for m = 1:length(x)-1
    dx = x(m+1) - x(m);
    dy = y(m+1) - y(m);
     s = s + sqrt(dx^2 + dy^2);
end



% Find the the points (x1,y1) and (x2,y2) between which the arrow will be
% drawn.

temp = 0;
flag_x1 = 0;

% Find the location where the arc length = desired fraction of total arc
% length
for m = 1:length(x)-1
    dx = x(m+1) - x(m);
    dy = y(m+1) - y(m);
    
    temp = temp + sqrt(dx^2 + dy^2);
    if temp > arrow_pos_fct*s & flag_x1 == 0
       flag_x1 = 1;
       x1 = x(m);
       y1 = y(m);
       ind = m;
    end
    
    if flag_x1 == 1 &  m > ind & temp > (arrow_pos_fct - arrow_len_fct)*s
        x2 = x(m);
        y2 = y(m);
        break;
    end
    
end




% Find the direction of the arrows from the field vectors

ux1 = interp2(X,Y,ux,x1,y1);
uy1 = interp2(X,Y,uy,x1,y1);

dot_prd = ux1*(x2-x1) + uy1*(y2-y1);  % calculate the dot product
drct = sign(dot_prd);    % Sign of dot product says if the direction of the arrow.



xp1 = x1;
yp1 = y1;
xp2 = drct*(x2-x1);    % Change direction depending on the dot product result
yp2 = drct*(y2-y1);

quiver(xp1,yp1,xp2,yp2,2,'MaxHeadSize',14,'linewidth',2,'color',clr,'autoscale','off');



