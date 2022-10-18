function c = recover_coefficient_handler(SUBDOMAIN, j, c_bar, OPTIONS)

    if(OPTIONS.RADIATION>0)
        c = exterior_problem(SUBDOMAIN, j, c_bar, OPTIONS);
    else
        c = interior_problem(SUBDOMAIN, j, c_bar);
    end
end


function c = exterior_problem(SUBDOMAIN, j, c_bar, OPTIONS)
    M = SUBDOMAIN.params.M(j);
    k = SUBDOMAIN.params.k(1);
    c = zeros(M,2);

    if(OPTIONS.RADIATION == 1) % Sommerfeld condition
        % alpha = ik, beta = 1
        c(:,1) = c_bar;
        c(:,2) = -1i*k*c_bar;
        
    elseif(OPTIONS.RADIATION == 2)
        I = eye(M); A = zeros(M);
        for j=1:M
            A(j,:) = transpose(cheb_derivatives(I(:,j),2));
        end
        B = (k^2 * I + (1/2)*transpose(A)) / (1i*k);        
        
        c(:,1) = c_bar;
        c(:,2) = B * c_bar;
        
    elseif(OPTIONS.RADIATION == 3)
        I = eye(M); A = zeros(M);
        for j=1:M
            A(j,:) = transpose(cheb_derivatives(I(:,j),2));
        end
        B = (k^2 * I + (1/4)*transpose(A)) \ (-1i*k^3 * I - (3/4)*1i*k*transpose(A));
        
        c(:,1) = c_bar;
        c(:,2) = B*c_bar;
        
    elseif(OPTIONS.RADIATION == 4)
        I = eye(M); A2 = zeros(M); A4 = zeros(M);
        for j=1:M
            A2(j,:) = transpose(cheb_derivatives(I(:,j),2));
            A4(j,:) = transpose(cheb_derivatives(I(:,j),4));
        end
        B = (1i*k^3 * I + (1/2)*1i*k*transpose(A2))\(k^4 * I + ...
                k^2 * transpose(A2) + (1/8)*transpose(A4));
        c(:,1) = c_bar;
        c(:,2) = B*c_bar;
        
    else
        N = (M-1)/2;

        [a,b] = mode_coefficients(N, SUBDOMAIN.params.k(1), OPTIONS);

        c(:,1) = c_bar;
        c(:,2) = (-(a./b)) .* c_bar;
    end
end

function c = interior_problem(SUBDOMAIN, j, c_bar)
    a = SUBDOMAIN.boundary(j).alpha;
    b = SUBDOMAIN.boundary(j).beta;
    M = SUBDOMAIN.params.M(j);
    
    c = zeros(M,2);
    d = basiscoeffs(SUBDOMAIN.boundary(j).fun, M, SUBDOMAIN.basis(j));
    if(a~=0)
        % Recover c0 via the formula, either Dirichlet or
        % proper Robin case.
        c(:,1) = (1/a)*d - (b/a)*c_bar;
        c(:,2) = c_bar;
    else
        % This is only the proper Neumann case.
        c(:,1) = c_bar;
        c(:,2) = (1/b)*d;
    end

end