%-------------------------------------------------------------------------
%
%   Function:   N_plus_potential
%
%   Inputs:    xi - (N+1) x (N+1)
%              Mp - logical (N+1 x N+1)
%              Np - logical (N+1 x N+1)
%               x - (N+1) length vector
%               y - (N+1) length vector
%
%   Outputs:    p - potential
%
%   Description:
%
%-------------------------------------------------------------------------

function [p,solve_time] = N_plus_potential_Q(xi,BLOCK,SUBDOMAIN,REAL)

    N = BLOCK.N;
    h = BLOCK.grids.x(2) - BLOCK.grids.x(1);
    k = SUBDOMAIN.params.k(1);
    
    % Apply Lh and truncate to M+
    gr = apply_L(real(xi), k, N, h);
    gi = apply_L(imag(xi), k, N, h);
    g = complex(gr, gi);
%     g = apply_L(xi, k, N, h);
    g(~BLOCK.grids.Mp) = 0;

    % Solve the AP
    SOLVE = tic;
    if(REAL)
%         solved = REAL_solver(g,SUBDOMAIN.k,x,y,true);
    else
        solved = solver(g, k, N, h,true);
    end
    solve_time = toc(SOLVE);
    
    % Last steps of computing difference potential
    p = (xi - solved);
    p(~BLOCK.grids.Np) = 0;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Applies the discrete differential operator (stencil) L.
%   This is not a solver, and in fact performs the task that
%   a solver traditionally inverts. Created for the helmholtz
%   equation L = \laplace{u} + k^2 u
%
%   Input:  u - (M+1)x(N+1) matrix
%      params - struct containing at least: N, h, k
%
%   Output: g - (M+1)x(N+1) matrix of stenciled values.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function g = apply_L(u, k, N, h)
    % N = [N1, N2] is dimensions of x and y discretizations
    % h is grid spacing
    % k is wavenumber

    % initialize the output
    g = complex(zeros(N));
    
    % convenience def
    h2 = h^2;

    % workhorse loop
    for j = 2:N-1
        for i = 2:N-1
            g(i,j) = (u(i+1,j+1) + u(i-1,j+1) + u(i+1,j-1) + u(i-1,j-1))*(1/(6*h2)) + ...
                     (u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1))*(2/(3*h2) + (k^2)/12) + ...
                     u(i,j) * (-10/(3*h2) + (8*k^2)/12);
        end
    end
end