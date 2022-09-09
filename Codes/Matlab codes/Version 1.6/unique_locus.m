% Mohammad Asif Zaman
% May 21, 2019
% Finding unique integral lines

% Input: x{}, y{}: A cell array of (x,y) loci, 
%        dist_tol = minimum separation between two loci for them to be
%        considered unique
% Output: x_out{}, y_out{}: Output cell arrays of unique loci

function [x_out,y_out] = unique_locus(x_in,y_in,dist_tol)


% Sort the cell-arrays according to length. The shortest one appears first
[~,I] = sort(cellfun(@length,x_in));
x = x_in(I);
y = y_in(I);


counter = 0;  % counter for unique loci
% dist_tol = 1e-2;



for m = 1:length(x)-1   % for each locus
%     clear xt;
%     clear yt;
    
    xt = x{m};
    yt = y{m};
    min_dist = inf;
    
    
    for n = m+1:length(x)  % compare with all other locuses
%         clear xtemp;
%         clear ytemp;
        xytemp = sortrows([x{n} y{n}],1);  % This sorting is necessary for the interpolation function

        xtemp = xytemp(:,1);
        ytemp = xytemp(:,2);
        
        
        % Check if the x values have some overlap.
        if max(xtemp) < min(xt) | max(xt) < min(xtemp)
            flg = 1;
%             fprintf('flag for m = %d and n = %d \n',m,n);
        else
            
            fy = griddedInterpolant(xtemp,ytemp,'nearest','nearest');   % Scattered interpolant doesn't work on 1D. 
            yt2 = fy(xt);
            
            dist = sum((yt-yt2).^2)./length(yt);
%             dist_store(m,n) = dist;
            
            if dist < min_dist
                min_dist = dist;
                ind_store = n;
            end
        end
        
        
    end
    
    
    
    
    %     min_dist
    if min_dist > dist_tol & min_dist < inf
        fprintf('Locus %d is unique (Minimum dist = %1.4f)\n',m,min_dist);
        counter = counter + 1;
        x_out{counter} = xt;
        y_out{counter} = yt;
    else
        fprintf('Locus %d is same as %d locus (dist = %1.4f)\n',m,ind_store,min_dist);
        
        
    end
    
    
end

% Include the last locus
x_out{counter + 1} = x{end};
y_out{counter + 1} = y{end};


fprintf('No. of unique loci = %d\n',length(x_out));

% min(dist(:))
% max(dist(:))
% dist_store(1,:)'