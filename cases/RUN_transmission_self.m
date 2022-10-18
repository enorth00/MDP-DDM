%% Case 1. Empty and rod, centered on empty.
clear all;

% k = 1 to 3

OPTIONS.T = 2;

Mc = 25;
Mr = 31;

k1 = 1;
k2 = 3;

Mf=31;  M = [Mc, Mf, Mr];
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];
DOMAIN = [1, 2];
start = [-1, 1, -1, 1];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);


OPTIONS.GRIDS = [6:10];
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


%% Case 2: k = 3 to 10
OPTIONS.T = 2;

Mc = 25;
Mr = 31;

k1 = 3;
k2 = 10;

Mf=31;  M = [Mc, Mf, Mr];
TYPES = cell(1,3);
TYPES(1,:) = [{1}, {'cccc'}, {[k1, 0]}];
TYPES(2,:) = [{2}, {'cccc'}, {[k1, k2]}];
DOMAIN = [1, 2];
start = [-1, 1, -1, 1];

OPTIONS.Ns = numel(DOMAIN);
OPTIONS.TYPES = TYPES;
BLOCK = allocate_blocks([1, 2]);
SUBDOMAIN = allocate_subdomains(OPTIONS.Ns);
[SUBDOMAIN, OPTIONS.ADJACENCY] = subdomain_construction(SUBDOMAIN, DOMAIN, TYPES, M, start);


OPTIONS.GRIDS = [6:10];
OPTIONS.FUNC = 14;
OPTIONS.PARALLEL = true;
OPTIONS.REAL = false;
OPTIONS.RADIATION = 2;
OPTIONS.BOUNDARY.alpha = ones(OPTIONS.Ns,4);
OPTIONS.BOUNDARY.beta = zeros(OPTIONS.Ns,4);
OPTIONS.FOURIER_AUX = [-3*pi/8, 3*pi/8];

run('DRIVER_self');
errors2 = errors;
ratios2 = errors2(1:end-1) ./ errors2(2:end);


%% output a latex table
n = (2 .^ (OPTIONS.GRIDS - ones(1, length(OPTIONS.GRIDS))))';
output = [n, errors1, [0; log2(ratios1)], errors2, [0; log2(ratios2)]];

formatspec = "%i & %.2e & %.2f & %.2e & %.2f \\\\\n";

fprintf(formatspec, output');

