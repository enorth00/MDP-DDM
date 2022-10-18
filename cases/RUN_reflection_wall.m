% 3x3 empty subdomains, k=40, varying values of M.
clear all;

%% Set everything for the first one.

OPTIONS.T = 1;

Mc = 35;

k1 = 40;
k2 = 0;

Mf=31; Mr=35; M = [Mc, Mf, Mr];
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
DOMAIN = 1*ones(3); 
start = [-3, -1, 1, 3];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);


OPTIONS.GRIDS = [6];
OPTIONS.FUNC = 14;
OPTIONS.PARALLEL = false;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 2;
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

run('DRIVER_exact');

errors1 = errors;

% Change Mc and OPTIONS.GRIDS. Run again.
Mc = 45;
M = [Mc, Mf, Mr];
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);

OPTIONS.GRIDS = [7, 8];

run('DRIVER_exact');

errors23 = errors;

% Change Mc (and relevant followup pieces) and OPTIONS.GRIDS. Run again.
Mc = 60;
M = [Mc, Mf, Mr];
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);

OPTIONS.GRIDS = [9, 10];

run('DRIVER_exact');

errors45 = errors;

% Process final results.
errors = [errors1; errors23; errors45];
ratios = errors(1:4) ./ errors(2:5);

% If you have access to plotting, also run:
% PLOT_reflection_case;

