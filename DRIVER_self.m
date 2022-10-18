% DRIVER_self.m

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

% T = length(OPTIONS.types);
errors = zeros(length(OPTIONS.GRIDS), 1);

for n = 1:length(OPTIONS.GRIDS)
    N = 2^(OPTIONS.GRIDS(n)) + 1;
    
    % Prep the grid-size-dependent parts of reference blocks
    for i=1:2
        BLOCK(i).N = [N,N];
        BLOCK(i).grids = SET_grids(BLOCK(i));
    end
    
    % Run the solving algorithm
    [SOLUTION, QR_time, solve_time] = problem_solver(SUBDOMAIN, BLOCK, OPTIONS);
    
    % Compute errors (in self-norm)
    if(n>1)
        errors(n) = self_error(SOLUTION, LAST, OPTIONS);
    end
    LAST = SOLUTION;
%     run('PLOT_solution');
%     run('PLOT_scattered');
end

if(n>2)
    ratios = errors(2:end-1)./errors(3:end);
    disp([errors, [0; 0; log2(ratios)]]);
else
    disp(errors);
end

