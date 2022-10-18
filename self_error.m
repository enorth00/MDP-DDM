function e = self_error(SOLUTION, LAST, OPTIONS)
    Ns = length(SOLUTION);
    error = zeros(Ns, 1);
    if(OPTIONS.PARALLEL)
        parfor i=1:Ns
            [N1, N2] = size(SOLUTION(i).U);
            ERR = LAST(i).U - SOLUTION(i).U(1:2:N1, 1:2:N2);
            error(i) = max(max(abs(real(ERR))));
        end
    else
    	for i=1:Ns
            [N1, N2] = size(SOLUTION(i).U);
            ERR = LAST(i).U - SOLUTION(i).U(1:2:N1,1:2:N2);
            error(i) = max(max(abs(real(ERR))));
        end
    end
    e = max(error);
end
