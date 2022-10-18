% DRIVER_exact.m script

% assume all relevant parameters, subdomains, etc.. have been defined.

if(OPTIONS.PARALLEL)
    parfor i = 1:OPTIONS.Ns
        SUBDOMAIN(i) = SET_SubdomainInfo(SUBDOMAIN(i), ...
                                        [OPTIONS.BOUNDARY.alpha(i,:); ...
                                         OPTIONS.BOUNDARY.beta(i,:)], ...
                                         OPTIONS.ADJACENCY(i,:), OPTIONS);
    end
else   
    for i = 1:OPTIONS.Ns
        SUBDOMAIN(i) = SET_SubdomainInfo(SUBDOMAIN(i), ...
                                        [OPTIONS.BOUNDARY.alpha(i,:); ...
                                         OPTIONS.BOUNDARY.beta(i,:)], ...
                                         OPTIONS.ADJACENCY(i,:), OPTIONS);
    end
end

temp = false;
errors = zeros(length(OPTIONS.GRIDS),1);
for n = 1:length(OPTIONS.GRIDS)
    N = 2^(OPTIONS.GRIDS(n)) + 1;
    
    % Prep the grid-size-dependent parts of reference blocks
    for i=[1,2]
        BLOCK(i).N = [N,N];
        BLOCK(i).grids = SET_grids(BLOCK(i));
    end
    
    % Run the solving algorithm
    [SOLUTION, QR_time, solve_time] = problem_solver(SUBDOMAIN, BLOCK, OPTIONS);
    if(OPTIONS.FUNC == 14)
        temp = true;
        OPTIONS.FUNC = 13;
    end
    % Compute errors
    E = zeros(1,OPTIONS.Ns);
    u_ex = cell(OPTIONS.Ns,1);
    for i=1:OPTIONS.Ns
        T = SUBDOMAIN(i).type;
        
        u_ex{i} = get_exact(BLOCK(T), SUBDOMAIN(i), BLOCK(T).grids.Mp, OPTIONS);
        if(T==2)
            u_ex{i} = u_ex{i} + get_exact(BLOCK(T), SUBDOMAIN(i), BLOCK(T).grids.rodMp, OPTIONS);
        end
        E(i) = max(max(abs(real(u_ex{i} - SOLUTION(i).U))));
%         E(i) = max(max(real((u_ex{i} - SOLUTION(i).U))));
    end
    errors(n) = max(E);
    if(temp)
        OPTIONS.FUNC = 14;
        temp = false;
    end
end
if(n>1)
    ratios = errors(1:n-1)./errors(2:n);
    disp([errors, [0; log2(ratios)]]);
else
    disp([errors]);
end

