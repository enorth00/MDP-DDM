function [Q_star, Q_bar] = boundary_condition_handler(SUBDOMAIN, Q0, Q1, adj, OPTIONS)
    
    if(OPTIONS.RADIATION>0)
        [Q_star, Q_bar] = exterior_problem(SUBDOMAIN, Q0, Q1, adj, OPTIONS);
    else
        [Q_star, Q_bar] = interior_problem(SUBDOMAIN, Q0, Q1, adj);
    end


end

function [Q_star, Q_bar] = interior_problem(SUBDOMAIN, Q0, Q1, adj)
    Q_bar = cell(5,1);
    [NN, ~] = size(Q0);
    Q_star = zeros(NN,1);
    M = SUBDOMAIN.params.M(1:4);
    MM = [0, M(1), sum(M(1:2)), sum(M(1:3)), sum(M(1:4))];
    
    for i=1:4
        if(adj(i) == 0)
            % Exterior side: Resolve boundary conditions 
            a = SUBDOMAIN.boundary(i).alpha; 
            b = SUBDOMAIN.boundary(i).beta;
            
            % extract coefficients from boundary data in the appropriate
            % basis
            d = basiscoeffs(SUBDOMAIN.boundary(i).fun, M(i), SUBDOMAIN.basis(i));

            if(a~=0)
                Q_star = Q_star + (1/a)*Q0(:,MM(i)+1:MM(i+1))*d;
                Q_bar(i) = {Q1(:,MM(i)+1:MM(i+1)) - (b/a)*Q0(:,MM(i)+1:MM(i+1))};
            else
                Q_star = Q_star + (1/b)*Q1(:,MM(i)+1:MM(i+1))*d;
                Q_bar(i) = {Q0(:,MM(i)+1:MM(i+1))};
            end
        else % adj(i) ~= 0
            % Interface case: keep both parts to sort out later.
            Q_bar(i) = {[Q0(:,MM(i)+1:MM(i+1)), Q1(:,MM(i)+1:MM(i+1))]};
        end
    end

end

function [Q_star, Q_bar] = exterior_problem(SUBDOMAIN, Q0, Q1, adj, OPTIONS)
    Q_bar = cell(5,1);
    [NN, ~] = size(Q0);
    Q_star = zeros(NN,1);
    M = SUBDOMAIN.params.M(1:4);
    MM = [0, M(1), sum(M(1:2)), sum(M(1:3)), sum(M(1:4))];
    k = SUBDOMAIN.params.k(1);
    for i=1:4
        if(adj(i)==0)
            
            if(OPTIONS.RADIATION == 1) % Sommerfeld condition
                % alpha = ik, beta = 1
                Q_bar(i) = {Q0(:,MM(i)+1:MM(i+1)) - 1i*k*Q1(:,MM(i)+1:MM(i+1))};

            elseif(OPTIONS.RADIATION == 2) % EM 2nd order
                I = eye(M(i)); A = zeros(M(i));
                for j=1:M(i)
                    A(j,:) = transpose(cheb_derivatives(I(:,j),2));
                end
                B = (k^2 * I + (1/2) * transpose(A)) / (1i*k);
                
                Q_bar(i) = {Q0(:,MM(i)+1:MM(i+1)) + Q1(:,MM(i)+1:MM(i+1))*B};

            elseif(OPTIONS.RADIATION == 3) % EM 3rd order
                I = eye(M(i)); A = zeros(M(i));
                for j=1:M(i)
                    A(j,:) = transpose(cheb_derivatives(I(:,j),2));
                end
                B = (k^2 * I + (1/4)*transpose(A))\(-1i*k^3 * I - (3/4)*1i*k*transpose(A));
                
                Q_bar(i) = {Q0(:,MM(i)+1:MM(i+1)) + Q1(:,MM(i)+1:MM(i+1))*B};
                
            elseif(OPTIONS.RADIATION == 4) % EM 4th order
                I = eye(M(i)); A2 = zeros(M(i)); A4 = zeros(M(i));
                for j=1:M(i)
                    A2(j,:) = transpose(cheb_derivatives(I(:,j),2));
                    A4(j,:) = transpose(cheb_derivatives(I(:,j),4));
                end
                B = (1i*k^3 * I + (1/2)*1i*k*transpose(A2))\(k^4 * I + ...
                        k^2 * transpose(A2) + (1/8)*transpose(A4));
                Q_bar(i) = {Q0(:,MM(i)+1:MM(i+1)) + Q1(:,MM(i)+1:MM(i+1))*B};

            else
                N = (M(i) - 1)/2; % number of modes present in the expansion
                % This subroutine constructs the radiation coefficients
                % when fourier basis functions are being used.
                [a,b] = mode_coefficients(N, k, OPTIONS);
                Q_bar(i) = {Q0(:,MM(i)+1:MM(i+1)) - ((a./b).').*Q1(:,MM(i)+1:MM(i+1))};
            end

        else % adj(i) ~= 0, interface
            Q_bar(i) = {[Q0(:,MM(i)+1:MM(i+1)), Q1(:,MM(i)+1:MM(i+1))]};
        end
        % END adj==0 IF
    end
    % END sides FOR
%     disp(cond(B));
end