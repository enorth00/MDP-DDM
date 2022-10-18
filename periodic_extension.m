
% The goal is for this function to accept, as input, a function and its
% relevant derivative (in one particular dimension). This function will
% return a function handle that can be used in place of the original
% function to evaluate the periodically extended version.

% This creates tails that are degree n-1.


function u = periodic_extension(U, X)

    x1 = X(1); x2 = X(2); xt1 = X(3); xt2 = X(4);
    % n is number of terms in polynomial extensions (degree n-1).
    n = 12;
    A = zeros(2*n); % n terms on each tail requires 2n conditions.
    
    m = 8; % number of conditions on the domain endpoints.
           % Then 2n - 2m is the remaining conditions to put on periodicity
           % continuity

    % Initialize the cascading numbers used for the powers in the
    % polynomial representation.
    P = n-1:-1:0;
    powers = zeros(n);
    powers(1,:) = ones(1,n);
    for i=1:n-1
        powers(i+1,:) = [P(i:end), zeros(1,i-1)];
    end

    % Fill in the main coefficient matrix for the physical domain endpoints
    for i=1:m
        A(i,:) = [prod(powers(1:i,:)).*x1.^(powers(i+1,:)), zeros(1,n)];
        A(m+i,:) = [zeros(1,n), prod(powers(1:i,:)).*x2.^(powers(i+1,:))];
    end
    
    % Fill in the main coefficient matrix for periodicity conditions.
    for i=1:2*(n-m)
        A(2*m+i,:) = [prod(powers(1:i,:)).*xt1.^powers(i+1,:), ...
                    -(prod(powers(1:i,:)).*xt2.^powers(i+1,:))];
    end

    % Construct the right-hand side. U is a chebfun so diff is efficient,
    % accurate, and stable.
    b = zeros(2*n,1); Udiff = U;
    for i=1:m
        b(i) = Udiff(x1);
        b(m+i) = Udiff(x2);
        
        Udiff = diff(Udiff);
    end

    % Solve
    c = A\b;

    % Return an in-line, piecewise defined function handle.
    u = @(x) ((x.^(n-1:-1:0)) * c(1:n)) .* (x<x1) + ...
              U(x).*((x>=x1) & (x<=x2)) + ...
             ((x.^(n-1:-1:0)) * c(n+1:2*n)) .* (x>x2);
end