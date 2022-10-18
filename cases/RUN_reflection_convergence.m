% Reflection errors and their relative percents at the exterior wall
clear all;

%% Set everything for the first one. 1x1 case

OPTIONS.T = 1;

Mc = 30;

k1 = 10;
k2 = 0;

Mf=31; Mr=35; M = [Mc, Mf, Mr];
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
DOMAIN = 1*ones(1); 
start = [-1, 1, -1, 1];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);


OPTIONS.GRIDS = [6:9];
OPTIONS.FUNC = 14;
OPTIONS.PARALLEL = true;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 2;
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

run('DRIVER_self');

errors1 = errors;
ratios1 = errors1(1:end-1) ./ errors1(2:end);





%% 3x3 case

DOMAIN = 1*ones(3);
start = [-3, -1, 1, 3];
OPTIONS.Ns = numel(DOMAIN);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);

run('DRIVER_self');

errors2 = errors;
ratios2 = errors2(1:end-1) ./ errors2(2:end);

%% 5x5 case

DOMAIN = 1*ones(5);
start = [-5, -3, 3, 5];
OPTIONS.Ns = numel(DOMAIN);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);

run('DRIVER_self');

errors3 = errors;
ratios3 = errors3(1:end-1) ./ errors3(2:end);



%% output a latex table
n = (2 .^ [6:9])';
output = [n, errors1, [0; log2(ratios1)], errors2, [0; log2(ratios2)], errors3, [0; log2(ratios3)]];

formatspec = "%i & %.2e & %.2f & %.2e & %.2f & %.2e & %.2f \\\\\n";

fprintf(formatspec, output');
