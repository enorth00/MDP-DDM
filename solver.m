function U = solver(f1, k, NN, h, P)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Solves the 2d constant coefficient Helmholtz equation with homogeneous
%   DBC on the y-boundaries and Sommerfeld radiation boundary conditions
%   on the x-boundaries.
%
%   Requires a tridiagonal solver "tri_diag.m" as well as PDE Toolbox
%   for (dst, idst)
%
%   Inputs -
%       f1: M+1xN+1 matrix of precomputed source values
%       k: wavenumber for the helmholtz equation
%       x,y: (1xM+1 and 1xN+1) 1D grids for evaluation
%       P: Flag for right-hand stencil. True (1) for pre-stenciled,
%           False (0) for needs stenciled.
%   Output -
%       U: (M+1xN+1) matrix. The finite difference solution
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--- Grid variables initialization ---%

    % h is given
    M = NN(1)-1; N = NN(2)-1;
    F = zeros(M+1,N-1);
    F = complex(F,0);
    U = zeros(M+1,N-1);

%--- Right-hand side stencil ---%

% Determines whether or not the right-hand side is pre-stenciled. If not,
% we stencil it here. If it is, just take the given input.
    if(~P)
        for m = 2:M
            for n = 1:N-1
                F(m,n) = 8*f1(m,n+1) + f1(m,n+2) + f1(m+1,n+1) ...
                    + f1(m,n) + f1(m-1,n+1);
            end
        end
        for n = 1:N-1
            F(1,n) = 8*f1(1,n+1) + f1(1,n+2) + f1(1,n) + f1(2,n+1);
            F(M+1,n) = 8*f1(M+1,n+1) + f1(M+1,n+2) + f1(M+1,n) + f1(M,n+1);
        end
        F = F/12;
    else
        F = f1(:,2:N);
    end

    F_hat = dst(F.'); F_hat = F_hat.';

%--- Precompute scheme coefficients ---%
    N = N-1; % there's a lot of N-1 that would float around, so we rename
    j_vec = 1:N;
    boundary_coeff = zeros(4,N);
    sin_terms = (sin((j_vec * pi)/(2*(N+1)))).^2;
    center_coeff = (5*k^2 / 6)*ones(1,N) - (2/h^2) * ones(1,N) ...
                    - (8/(3*h^2) + (k^2 / 3)) * sin_terms;
    edge_coeff = ((1/h^2) + (k^2 / 12))*ones(1,N) - (2/(3*h^2))*sin_terms;
    boundary_coeff(1,:) = ((-4/(3*h) + (k^2 * h)/24) - 1i*h^2*k^3/16)*ones(1,N) + ...
                          (1/(3*h) - 1i*k/2)*(ones(1,N) - 2*sin_terms);
    boundary_coeff(2,:) = -conj(boundary_coeff(1,:));
    boundary_coeff(3,:) = conj(boundary_coeff(1,:));
    boundary_coeff(4,:) = conj(boundary_coeff(2,:));

    N = N+1; % rename it back for the rest of the program

 %--- Run the loop and tridiagonal solve every time ---%
 o = complex(ones(M+1,1), zeros(M+1,1));
%  tic;
    for j = 1:N-1

        a = center_coeff(j)*o;
        b = edge_coeff(j)*o;
        c = b;

        a(1) = boundary_coeff(1,j);
        c(1) = boundary_coeff(2,j);
        a(M+1) = boundary_coeff(3,j);
        b(M+1) = boundary_coeff(4,j);

        U(:,j) = tridiag(a, b, c, F_hat(:,j));
    end

%--- Untransform and post-process the solution with its DBC ---%
    U = (idst(U.'));
    U = ([zeros(1,M+1); U; zeros(1,M+1)]).';

end