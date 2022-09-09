
% Mohammad Asif Zaman
% Vector Field Topology (2D)

% Version 1.3 April 25, 2019
% Version 1.4.4 April 29, 2019
% Version 1.5 May 24, 2019
% Version 1.6 April 7, 2020
%         Included the definition of the space variables within the
%         text_vector_field.m file.



% Functions required
% 1. tes_vector_field()   :  If a test vector field is required to be loaded
%                            v1.6 the space variables X,Y are also outputted now along with ux and uy.   
% 2. F_diff_mat_2D_v2(.)  :  Load the finite difference matrices
% 3. JacobianInterp(.)    :  A function that evaluates the Jacobian matrix at an arbitrary point
% 4. critical_class(.)    :  Classifies the critical points from corresponding Eigne values
% 5. int_locus(.)         :  Draws the integral lines
% 6. drawArrow(.)         :  Function to draw arrows on the integral lines


clear all;
close all;
clc; clf;



% ======================================================================>>>
% <<<======================================================================


% I/O
% ======================================================================>>>
path = 'Data/';
path_write = [path 'Data_out/'];
data_write = 'n';
locus_write = 'n';
% <<<======================================================================


% Control parameters
% ======================================================================>>>
ctrl_locus_test = 'n';

% <<<======================================================================



% Analysis parameters
% ======================================================================>>>
N_init_guess = 20;% Guess for number of critical points

min_tol = 2e-2;   % Defining a threshold value.
                  % If the difference between a minima and the minimum values is smaller than
                  % this, then the minima is considered a global minima.

rlim = 1000;      % allowed maximum radial extent of critical points 
                  % (critical points outside this radius will be ignored)
                  
min_dist_tol = 2; % If separation between 2 critical points is less than 
                  % min_dist_tol * dx, then the 2 points are considered
                  % identical.

locus_dist_tol = 1e-2;   % Minimum RMS separation between loci for them to be considered unique
% <<<======================================================================



%  Loading the vector field and defining the spatial grid
% ======================================================================>>>

[X,Y,ux,uy] = test_vector_field();     % Load a test vector field


% Lower and upper bound of x and y
xL = min(X(:)); xU = max(X(:)); 
yL = min(Y(:)); yU = max(Y(:)); 

[Ny,Nx] = size(ux);     % Size of the vector field


u2 = ux.^2 + uy.^2;     % Squared norm function

quiver(X,Y,ux,uy);
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











% Jacobian matrix calculation
% ======================================================================>>>

% The Jacobian matrix is calculated at all points in the grid. This will
% be used as a reference dataset. The Jacobian matrix at an arbitrary point
% in space will be calculated from this dataset using 2D interpolation.

dx = X(1,2) - X(1,1);
dy = Y(2,1) - Y(1,1);

% D2x = d/dx, DD2x = d^2/dx^2
[D2x,D2y,DD2x,DD2y] = F_diff_mat_2D_v2(Ny,Nx);  % Loading the finite difference matrices

J11 = reshape(D2x*ux(:)./(2*dx), size(X));  % dux/dx
J12 = reshape(D2y*ux(:)./(2*dy), size(X));  % dux/dy
J21 = reshape(D2x*uy(:)./(2*dx), size(X));  % duy/dx
J22 = reshape(D2y*uy(:)./(2*dy), size(X));  % dut/dy



% figure,
% contourf(X,Y,J11,40,'linestyle','none');
%
% figure,
% contourf(X,Y,J12,40,'linestyle','none');
%
% figure,
% contourf(X,Y,J21,40,'linestyle','none');
%
% figure,
% contourf(X,Y,J22,40,'linestyle','none');


% <<<======================================================================








% Calculating the critical points: April 23, 2019
% ======================================================================>>>
tic

% Initial guesses
% ======================================================================>>>

fprintf('Critical point search:\n');
fprintf('==============================================\n');


counter = 0;
% Find local minima of |u|^2. This will be a good approximate for the zeros
% of [ux uy]' = 0.
% Note: April 7, 2020: A better approach for finding initial guesses should
% needs to be developed.

% ======================================================
% April 8, 2020: Disabled
% ======================================================

% If the value of a function at a point is lower than its value at neighboring 4 points, then it is a local minima.



% for m = 2:Nx-1
%     for n = 2:Ny-1
%         if u2(n,m) < u2(n,m-1) & u2(n,m) < u2(n,m+1) & u2(n,m) < u2(n-1,m) & u2(n,m) < u2(n+1,m)
%             counter = counter + 1;
%             xc_guess(counter) = X(n,m);
%             yc_guess(counter) = Y(n,m);
%             xc_ind(counter) = m;
%             yc_ind(counter) = n;
%             
%         end
%     end
% end
% 
% fprintf('Initial guess of no. of critical points = %d\n',counter);
% Mtemp = [xc_guess' yc_guess'];
% fprintf('(%1.4f %1.4f) \n',Mtemp' );



% =========================================================================
% April 8, 2020: New approach for finding the guesses for the critical
% points
% =========================================================================

% Sorting the norm function to find the minimum values and their location
fprintf('Guesses of critical points (%i)\n', N_init_guess);
[u2s, ind_s] = sort(u2(:));

for m = 1:N_init_guess
    counter = counter + 1;
    xc_guess(counter) = X(ind_s(m));
    yc_guess(counter) = Y(ind_s(m));
    fprintf('(%1.4f %1.4f) \n',xc_guess(counter),yc_guess(counter) );
end
    
% <<<======================================================================


% Revising based on difference with global minimum
% ======================================================================>>>

% Verifying whether the minimas are zeros (global) or local vallies
% This step needs more work
u2min = min(u2(:));

fprintf('\nRevising\n');
fprintf('===================================================\n')

counter2 = 0;
for m = 1:length(xc_guess)
    temp = interp2(X,Y,u2,xc_guess(m),yc_guess(m));
    fprintf('(%1.4f,%1.4f) -> diff. from minimum = %1.4f', xc_guess(m),yc_guess(m),abs(temp-u2min));
    if abs(temp-u2min) <=  min_tol & xc_guess(m)^2 + yc_guess(m)^2 < rlim^2
        fprintf(': Accepted');
        counter2 = counter2 + 1;
        xc_guess_v(counter2) = xc_guess(m);
        yc_guess_v(counter2) = yc_guess(m);
    else
        fprintf(': Rejected');
    end
    fprintf('\n');
end

fprintf('\nRevised no. of critical points = %d\n' ,counter2);
Mtemp = [xc_guess_v' yc_guess_v'];
fprintf('(%1.4f %1.4f) \n',Mtemp' );
fprintf('===================================================\n')
% <<<======================================================================



% Refining roots using Newton-Raphson method
% ======================================================================>>>

% xc_guess_v = xc_guess;
% yc_guess_v = yc_guess;


% 2D Newton-Raphson method to refine the position of the critical points
% Finding zero of u = [ux uy]' = 0; J = [uxx uxy; uyx uyy]; delta = -J^-1*u
% x_new = x_old + delta



% xc_guess_v = 0;
% yc_guess_v = 0;
counter = 0;
for n_critical = 1:length(xc_guess_v)
  
    x0 = xc_guess_v(n_critical);
    y0 = yc_guess_v(n_critical);
    
    fprintf('\nN-R refinement of critical point %i (%1.3f,%1.3f) \n' ,n_critical, x0, y0);

    
    xc_temp = x0; yc_temp = y0;
%     [xc_temp,yc_temp] = GridSearch(X,Y,ux,uy,x0,y0);
    [xc_temp,yc_temp,err_flag] = NewtonRaphson2D(X,Y,ux,uy,J11,J12,J21,J22,xc_temp,yc_temp);  % current critical point
    
    
    % Check if the current critical point is very close to a previously
    % calculated critical point. In such a case, set flag = 1 and do not
    % consider the current point as a new critical point.
    min_dist = inf;
    
    % calculate distance from the current critical point to all other
    % critical points. Then find the minimum distance.
    for q = 1:counter                                                      % Loop over all previous critical points
        dist = sqrt( (xc_temp - xc(q)).^2 + (yc_temp - yc(q)).^2 );        % Distance 
        if dist < min_dist
            min_dist = dist;
        end
    end
    
    if min_dist < min_dist_tol*dx
        err_flag = 1;
        fprintf('Duplicate found. Separation = %1.4f \n',min_dist);
    end
    
    if err_flag ~= 1
        
        counter = counter + 1;
        xc(counter) = xc_temp;
        yc(counter) = yc_temp;
    end
    
    fprintf('Refined root: (%1.4f, %1.4f)\n\n', xc(counter),yc(counter));
end


% xc = 0;
% yc = 0;

tcritical = toc;
fprintf('Final no. of critical points found = %d \n',length(xc));
Mtemp = [xc' yc'];
fprintf('(%1.4f %1.4f) \n',Mtemp' );
fprintf('Time required for calculating the critical points = %1.2f sec\n',tcritical);
% <<<======================================================================

% <<<======================================================================





% Critical point classification
% ======================================================================>>>

% 1 = att_focus, 2 = rep_focus, 3 = saddle, 4 = center, 5 = att_node
% 6 = rep_node


% Color array. Associating a different color for each of the critical
% point types.
clr_arr = [0.0   0.4   0.8;
    1.0   0.0   0.0;
    0.0   0.6   0.0;
    1.0   1.0   0.0;
    0.5   0.5   0.5;
    0.8   0.0   0.8;
    ];


txt_arr = {'Att. focus', 'Rep. focus', 'Saddle', 'Center', 'Att. node','Rep. node'};


for m = 1:length(xc)
    
    temp_J = JacobianInterp(J11,J12,J21,J22,X,Y,xc(m),yc(m));  % Calculating the Jacobian matrix at critical point
    [temp_V,temp_D] = eig(temp_J);  % Calculating the Eigen values and Eigen vectors
    
    temp_lambda = diag(temp_D);     % Eigen values
    
    cr_clss(m) = critical_class(temp_lambda);  % Classifying the critical points
    
    figure(1); hold on;
    plot(xc(m),yc(m),'Marker','o','MarkerFaceColor',clr_arr(cr_clss(m),:),'markersize',8,'color','none');
    axis([xL xU yL yU]);
    
    daspect([1 1 1]);
    set(gca,'Position',[0.1 0.1 .65 .65]);
end

if data_write == 'y'
    Mdata = [xc' yc' cr_clss'];
    fname = [path_write 'critical_points' num2str(m) '.csv'];
    csvwrite(fname, Mdata);
end
    
% <<<======================================================================





% Calculating the integral lines: April 23, 20119
% ======================================================================>>>
tic
xsink = xc(cr_clss == 1 | cr_clss == 5);
ysink = yc(cr_clss == 1 | cr_clss == 5);

cnt = 0;
for m = 1:length(xc)
    
    % Check critical point class. For sinks, the direction has to be
    % reversed (the lines propagate in the opposite direction of the field
    % vector).
    if cr_clss(m) == 1 | cr_clss(m) == 4 | cr_clss(m) == 5
        drct = -1;
    else
        drct = 1;
    end
    
    
    x0 = xc(m);
    y0 = yc(m);
    J0 = JacobianInterp(J11,J12,J21,J22,X,Y,x0,y0);
    [temp_V,temp_D] = eig(temp_J);
    
    
    temp_V = real(temp_V);    % This is necessary as eigen vectors can be complex
    % Setting the intial direction of the integral lines to be parallel to
    % the eigen vectors.
    
    
    
    if rank(temp_V) == 2
        % For each vector, draw 2 lines in +/- direction. So, a total of 4
        % lines per critical point. This works for distinct vectors.
        xslope_set = [temp_V(1,1) -temp_V(1,1) temp_V(1,2) -temp_V(1,2)];
        yslope_set = [temp_V(2,1) -temp_V(2,1) temp_V(2,2) -temp_V(2,2)];
    else
        % If the matrix is not full rank, then both eigen vectors are the
        % same. Then you can just have 2 lines from them. To draw the
        % orthogonal lines to the vector [a,b] you need the vector [-b,a]
        % or [b, -a] which will create the other 2 lines.
        xslope_set = [temp_V(1,1) -temp_V(1,1) temp_V(2,1) -temp_V(2,1)];
        yslope_set = [temp_V(2,1) -temp_V(2,1) -temp_V(1,1) temp_V(1,1)];
    end
    
    
    
    for eg_vct = 1:length(xslope_set)
        xslope = xslope_set(eg_vct);
        yslope = yslope_set(eg_vct);
        
        
        
        xy_locus = int_locus(X,Y,ux,uy,x0,y0,xslope,yslope,xsink,ysink,drct);
        
        cnt = cnt + 1;
        xlocus{cnt} = xy_locus(:,1);
        ylocus{cnt} = xy_locus(:,2);
        cr_clss_save(cnt) = cr_clss(m);
        
        % For saddle points, draw in +drct and -drct directions. There are
        % redundancies though. For each direction, there are just 2 lines,
        % instead of 4.
        if cr_clss(m) == 3
            xy_locus = int_locus(X,Y,ux,uy,x0,y0,xslope,yslope,xsink,ysink,-drct);
            cnt = cnt + 1;
            xlocus{cnt} = xy_locus(:,1);
            ylocus{cnt} = xy_locus(:,2);
            cr_clss_save(cnt) = cr_clss(m);
        end
        
    end
    
    fprintf('Finished calculating integral lines for critical point %d \n',m);
    
end
tlocus = toc;
fprintf('Time required for calculating the integral lines = %1.2f sec\n',tlocus);

% <<<======================================================================




% Finding the unique locuses
% ======================================================================>>>

if ctrl_locus_test == 'y'
    [xlocus,ylocus] = unique_locus(xlocus,ylocus,locus_dist_tol);
end

% <<<======================================================================


% Plotting the integral lines
% ======================================================================>>>

figure(1);
hold on;
for m = 1:length(xlocus)
    plot(xlocus{m},ylocus{m},'k');
    
    drawArrow(xlocus{m},ylocus{m},X,Y,ux,uy,clr_arr(cr_clss_save(m),:));
    
    if data_write == 'y'
        Mdata = [xlocus{m} ylocus{m}];
        fname = [path_write 'locus' num2str(m) '.csv'];
        csvwrite(fname, Mdata);
    end
    
end

xlabel('x'); ylabel('y');
axis([xL xU yL yU]);

daspect([1 1 1]);
set(gca,'Position',[0.1 0.1 .65 .65]);


for m = 1:6
    annotation('ellipse',[.82 m/10+.05 .015 .02],'Facecolor',clr_arr(m,:),'color','none');
    annotation('textbox',[.83 m/10+.085 0 0],'string',txt_arr{m});
    
end
pbaspect([1 1 1]);

% <<<======================================================================


% Matrix write
% ======================================================================>>>
if locus_write == 'y'
    sz1 = length(xlocus);
    sz2 = length(xlocus{end});
    xM = zeros(sz2,sz1);
    yM = xM;
    
    for m = 1:length(xlocus)
        for n = 1:length(xlocus{m})
            xM(n,m) = xlocus{m}(n);
            yM(n,m) = xlocus{m}(n);
        end
    end
    
    csvwrite('xM.csv',xM);
    csvwrite('yM.csv',yM);
    
end
        

% <<<======================================================================






