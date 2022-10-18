function b = cheb_derivatives(a, q)
    switch(q)
        case 1
            b = first_derivative(a);
        case 2
            b = second_derivative(a);
        case 3
            b = third_derivative(a);
        case 4
            b = fourth_derivative(a);
        otherwise
            disp('Invalid derivative request in cheb_derivatives');
    end
end

function b = first_derivative(a)
    M = length(a);
    b = zeros(M,1);

end

function b = second_derivative(a)
    M = length(a);
    b = zeros(M,1);

    for n=0:M-1
        CAP = floor((M-(n+1))/2);
        for j=1:CAP
            b(n+1) = b(n+1) + (j*(n+j)*(n+2*j))*a((n+1) + 2*j);
        end

        % In the formula, there is an extra factor of 16 outside the sum.
        b(n+1) = 4 * b(n+1);
        % The formula is built on the half-interval and we are taking 2
        % derivatives, so we need to divide out a factor of 2, twice. This
        % leaves an ultimate factor of only 4.
    end
    
    % The first coefficient needs halved in the formula
    b(1) = b(1)/2;
end

function b = third_derivative(a)
    M = length(a);
    b = zeros(M,1);

end

function b = fourth_derivative(a)
    M = length(a);
    b = zeros(M,1);

    for n=0:M-1
        CAP = floor((M-(n+1)-2)/2);
        for j=1:CAP
            c = (j+2)*(j+1)*j*(n+j+2)*(n+j+1)*(n+j)*(n+2*j+2);
            b(n+1) = b(n+1) + c*a((n+1) + 2*j + 2);
        end

        % In the formula, there is an extra factor of 128/3 outside the sum.
        b(n+1) = (8/3) * b(n+1);
        % The formula is built on the half-interval and we are taking 4
        % derivatives, so we need to divide out a factor of 2, four times.
        % This leaves an ultimate factor of only 8/3.
    end
    
    % The first coefficient needs halved in the formula
    b(1) = b(1)/2;

end