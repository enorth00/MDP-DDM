%%% Cases:
% 1) single x-value, vector of y-values
% 2) single y-vlaue, vector of x-values
% 3) two vectors of x&y-values generated from polar coords (single r-value
% and vector of theta-values
%
% The remedy is to figure out which coordinate (x or y) is the fixed value
% to be evaluated alongside the other vector. If case 3, then we are
% already given the two vectors of the appropriate pairs.

function F = radiative_source(X,Y,k,dx,dy)

    N1 = length(X); N2 = length(Y);
    if(N1==1)
        %%% Case 1
        X = X*ones(N2,1);        
    elseif(N2==1)
        %%% Case 2
        Y = Y*ones(N1,1);
    else
        %%% Case 3
        % Nothing because X and Y are input how we want them.
    end
    N = max(N1, N2);
    F = zeros(N,1);
    
    for i=1:N
        x = X(i); y = Y(i);
        R = 1/2;
        r = sqrt(x^2+y^2);    
        if(r >= R)
            f = 0;
        elseif(dx==0 && dy==0)
            f = radiative_solution(x,y,k,2,0) + radiative_solution(x,y,k,0,2) + ...
                k^2 * radiative_solution(x,y,k,0,0);
        elseif(dx==1 && dy==0)
            f = radiative_solution(x,y,k,3,0) + radiative_solution(x,y,k,1,2) + ...
                k^2 * radiative_solution(x,y,k,1,0);
        elseif(dx==0 && dy==1)
            f = radiative_solution(x,y,k,2,1) + radiative_solution(x,y,k,0,3) + ...
                k^2 * radiative_solution(x,y,k,0,1);
        elseif(dx==2 && dy==0)
            f = radiative_solution(x,y,k,4,0) + radiative_solution(x,y,k,2,2) + ...
                k^2 * radiative_solution(x,y,k,2,0);
        elseif(dx==0 && dy==2)
            f = radiative_solution(x,y,k,2,2) + radiative_solution(x,y,k,0,4) + ...
                k^2 * radiative_solution(x,y,k,0,2);
        elseif(dx==1 && dy==1)
            f = radiative_solution(x,y,k,3,1) + radiative_solution(x,y,k,1,3) + ...
                k^2 * radiative_solution(x,y,k,1,1);
        end
        F(i) = f;
    end
end