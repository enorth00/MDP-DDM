function y = evaluate(f, x)
    % intentionally de-vectorizing function evaluation
    
    N = length(x);
    y = zeros(size(x));
    for i=1:N
        y(i) = f(x(i));
    end

end